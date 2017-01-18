import numpy as np
import math
from fractions import gcd
from PIL import Image
from scipy.signal import convolve2d
import imtools

from cv2 import VideoCapture

# uses an image to create a grid of background values
def find_bg(images, out=None, conv_len=2, bg_img=50, l1cal=False, clear_hotpix=False):

    vid = False
    
    def divisorGen(n):
        large_divisors = []
        for i in xrange(1, int(math.sqrt(n) + 1)):
            if n % i == 0:
                yield i
                if i*i != n:
                    large_divisors.append(n / i)
        for divisor in reversed(large_divisors):
            yield divisor

    n_img_bg = len(images)

    # establish grid dimensions
    if imtools.is_video(images[0]):
        vid = True
        print "Processing as video..."
        cap = VideoCapture(images[0])
        retval, frame = cap.read()
        bands = ['R','G','B']
        h,w,n_bands = frame.shape
        if bg_cutoff:
            cutoff = imtools.outlier_cutoff(frame.transpose(2,0,1), 2)
        cap.release()
        
        max_grid = np.zeros((h,w,n_bands), dtype=int)
        s_grid = np.zeros((h,w,n_bands), dtype=float)
        
    else:
        if bg_img and (len(images) > bg_img):
            l1test = images[bg_img:]
	    images = images[:bg_img]
        elif l1cal:
            l1cal = None
            print "Unable to perform L1 Calibration: too few images"
        im_grid = imtools.ImGrid(images[0])
        w,h = im_grid.width, im_grid.height
        bands = im_grid.bands
        n_bands = im_grid.n_bands
 
        max_grid = np.zeros((n_bands,h,w), dtype=int)
        s_grid = np.zeros((n_bands,h,w), dtype=int)
    

    # determine sampling resolution
    full_block_len = gcd(h,w)
    divisor_list = list(divisorGen(full_block_len))
    aspect_ratio = (w/full_block_len, h/full_block_len)
    print "Possible sample dimensions:"
    for i,d in enumerate(divisor_list):
        print " [%d] %d x %d" % (i+1, w/divisor_list[i], h/divisor_list[i])

    sample_block = divisor_list[int(raw_input("Select [1]-[%d]: " % len(divisor_list)))-1]
            

    print
    print "Calibrating S thresholds..."
    print "Processing levels..." 
    
    
    # find second largest ADC count of each pixel
    if vid:
        cap = VideoCapture(images[0])
        iframe = 0
        ret, frame = cap.read()
        while ret and cap.isOpened(): 
            s_grid = np.median([max_grid, s_grid, frame], axis=0).astype(int)
            max_grid = np.amax([max_grid, s_grid, frame], axis=0).astype(int)
            if (iframe+1) % 100 == 0:
                print " %d" % (iframe+1)
            iframe += 1
            ret, frame = cap.read()
            if bg_image and iframe >= bg_img: break
        s_grid = s_grid.transpose(2,0,1)
    else:
        print " 0/%d" % n_img_bg
        for i,im in enumerate(images):
            if (i+1) % 10 == 0:
                print " %d/%d" % (i+1,n_img_bg)
            im_grid = imtools.ImGrid(im)
            s_grid = np.median([max_grid, s_grid, im_grid], axis=0).astype(int)
            max_grid = np.amax([max_grid, s_grid, im_grid], axis=0).astype(int)

    print "Downsampling image..."

    s_grid = np.mean([s_grid[:,x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
    
    if conv_len:
        print "Applying convolution kernel..."

        kernel_side = 2*conv_len+1
        s_kernel = np.repeat(1, kernel_side**2).reshape((kernel_side,kernel_side))/float(kernel_side)**2
        convolved_grid = np.array([convolve2d(s_grid[cval], s_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)])
	if clear_hotpix:
            hot_pix_cutoff = 1.5
            hot_pix = (s_grid - convolved_grid >= hot_pix_cutoff)
            print "%d hot pixels found" % np.sum(hot_pix)
            s_grid = np.where(hot_pix, 256, convolved_grid)
            
        else:
            s_grid = convolved_grid

    # resize
    s_grid = np.repeat(np.repeat(s_grid, sample_block, axis=1), sample_block, axis=2)
    
    # calibrate L1
    if l1cal:
	print "Calibrating L1 threshold..."
	l1thresh = 256
	n_bins = 20
	failed = False
	
	while l1thresh > 3:
            if failed:
                s_grid *= 1.1
            max_array = np.zeros(n_bins)
            if vid:
                test_frames = 10/l1cal
                frames_passed = test_frames
                ret, frame = cap.read()
                iframe = 0
                while ret and cap.isOpened():
                    iframe += 1
                    max_minus_bg = np.amax(frame-s_grid)
                    if max_minus_bg >= n_bins-1:
                        max_array[n_bins-1] += 1
                    elif max_minus_bg >= 0:
                        max_array[int(max_minus_bg)] += 1
                    if iframe >= test_frames: break
                    ret, frame = cap.read()
            else:
                frames_passed = len(l1test)
                for i,im in enumerate(l1test):
                    if i%100==0:
                        print " %d/%d" % (i, len(l1test))
                    im_grid = imtools.ImGrid(im)
                    max_minus_bg = np.amax(im_grid-s_grid)
                    if max_minus_bg >= n_bins-1:
                        max_array[n_bins-1] += 1
                    elif max_minus_bg >= 0:
                        max_array[int(max_minus_bg)] += 1
            l1thresh = 0
            target_passed = frames_passed*l1cal
            frames_passed = np.sum(max_array)
            while frames_passed > target_passed:
                l1thresh += 1
                frames_passed -= max_array[l1thresh-1]
            print "L1 threshold = %d" % l1thresh
            failed = True
	s_grid += l1thresh
	

    
    if out:
	save_grid = np.where(s_grid > 255, 255, np.ceil(s_grid).astype(int))
        if n_bands == 1:
            s_img = Image.fromarray(save_grid[0].astype(np.uint8), mode='L')
            
        else:
            img_mode = ''.join(bands)
            s_img = Image.fromarray(save_grid.transpose(1,2,0).astype(np.uint8), mode=img_mode)

        # save as png
        print "Saving background as %s" % out
        s_img.save(out)

    if vid:
	cap.release();
    
    return s_grid

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Calculates the background levels for a set of images')
    parser.add_argument('--in', required=True, dest='infiles', nargs='+', help='Images used to set thresholds')
    parser.add_argument('--out', default='bg.png', help='Output file name')
    parser.add_argument("--conv_len", type=int, default=0, help='Distance to which pixels are included in averaging')
    parser.add_argument('--show', action='store_true', help='Display resulting threshold image')
    parser.add_argument('--bg_img', type=int, default=0, help='Limits number of images to be processed')
    parser.add_argument('--l1cal', type=float, help='Target rate of passing L1 threshold')
    parser.add_argument('--clear_hotpix', action='store_true', help='If convolved, raise hot pixel thresholds to 256')
    
    args = parser.parse_args()
    
    if args.l1cal and args.bg_img==0:
        args.bg_img = 50
    bg = find_bg(args.infiles, args.out, args.conv_len, args.bg_img, args.l1cal, args.clear_hotpix)
    if args.show:
       import imshow
       imshow.imshow(args.out)
