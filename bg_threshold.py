import numpy as np
import matplotlib.pyplot as plt
import math
from fractions import gcd
from PIL import Image
from scipy.signal import convolve2d
import imtools
import imshow
from cv2 import VideoCapture

# uses an image to create a grid of background values
def find_bg(images, out, conv_len=5, bg_cutoff=False, max_img=50):

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
        s_grid = np.zeros((h,w,n_bands), dtype=int)
        
    else:
        if max_img and (len(images) > max_img):
            images = images[:max_img]
        im_grid = imtools.ImGrid(images[0])
        w,h = im_grid.width, im_grid.height
        bands = im_grid.bands
        n_bands = im_grid.n_bands
        if bg_cutoff:
            cutoff = imtools.outlier_cutoff(im_grid, 2)
 
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
            print "Max pixel = %d" % np.amax(frame)
            s_grid = np.median([max_grid, s_grid, frame], axis=0).astype(int)
            max_grid = np.amax([max_grid, s_grid, frame], axis=0).astype(int)
            if (iframe+1) % 10 == 0:
                print " %d" % iframe+1
            iframe += 1
            ret, frame = cap.read()
            if iframe >= max_img: break
        cap.release()
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

    s_grid = np.amax([s_grid[:,x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
    
    # remove hot pixels and tracks
    if bg_cutoff:
        print "Removing thresholds above ",
        for cval in xrange(n_bands-1):
            print "%d," % cutoff[cval],
        print "%d" % cutoff[n_bands-1]
          
        mask_kernel = np.array([[1,1,1,1,1,1,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1],\
                                [1,0,0,0,0,0,1],[1,1,1,1,1,1,1]],dtype=float)/24.
        masked_grid = 0*s_grid
        cutoff = cutoff.reshape(n_bands,1,1)
        while np.any(np.amax(np.amax(s_grid, axis=1), axis=1) > cutoff):
            for cval in xrange(n_bands):
                masked_grid[cval] = convolve2d(s_grid[cval], mask_kernel, mode='same', boundary='symm')
            s_grid = np.where(s_grid <= cutoff, s_grid, cutoff)
    
    if conv_len:
        print "Applying convolution kernel..."

        kernel_side = 2*conv_len+1
        s_kernel = np.repeat(1, kernel_side**2).reshape((kernel_side,kernel_side))/float(kernel_side)**2
        convolved_grid = np.array([convolve2d(s_grid[cval], s_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)])
        s_grid = np.maximum(s_grid, convolved_grid)
        s_grid = np.ceil(s_grid+0.9).astype(int)

    # resize
    s_grid = np.repeat(np.repeat(s_grid, sample_block, axis=1), sample_block, axis=2)
    if n_bands == 1:
        img_mode = 'L'
    else:
        img_mode = ''.join(bands)
    s_img = Image.fromarray(s_grid.transpose(1,2,0).astype(np.uint8), mode=img_mode)

    # save as png
    img_name = out
    print "Saving background as %s" % img_name
    s_img.save(img_name)
    
    return s_grid

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Calculates the background levels for a set of images')
    parser.add_argument('--in', required=True, dest='infiles', nargs='+', help='Images used to set thresholds')
    parser.add_argument('--out', default='bg.png', help='Output file name')
    parser.add_argument("--conv_len", type=int, default=0, help='Distance to which pixels are included in averaging')
    parser.add_argument("--bg_cutoff", action='store_true', help='Removes tracks during sauto processing.')
    parser.add_argument('--show', action='store_true', help='Display resulting threshold image')
    parser.add_argument('--max_img', type=int, default=0, help='Limits number of images to be processed')
    
    
    args = parser.parse_args()
    
    bg = find_bg(args.infiles, args.out, args.conv_len, args.bg_cutoff, args.max_img)
    if args.show:
        imshow.imshow(args.out)      
 
