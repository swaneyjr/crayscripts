import numpy as np
import matplotlib.pyplot as plt
import math
from PIL import Image

# uses an image to create a grid of background values
def find_bg(images, out, conv_len=0, bg_cutoff=False):

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
    im = Image.open(images[0])
    w,h = im.width, im.height
    im_pix = w*h
    img_mode = im.mode
    n_bands = len(img_mode)
    #mean_grid = np.zeros((h, w, n_bands))
    #var_grid = np.zeros((h, w, n_bands))
    s_grid = np.zeros((h, w, n_bands),dtype=int)
    cutoff_vals = np.zeros(n_bands)

    # determine sampling resolution
    full_block_len = gcd(h,w)
    divisor_list = list(divisorGen(full_block_len))
    aspect_ratio = (w/full_block_len, h/full_block_len)
    print "Possible sample dimensions:"
    for i,d in enumerate(divisor_list):
        print " [%d] %d x %d" % (i+1, w/divisor_list[i], h/divisor_list[i])

    sample_block = divisor_list[int(raw_input("Select [1]-[%d]: " % len(divisor_list)))-1]

    # set cutoff for tracks
    if bg_cutoff:
        imarray = np.array(im)
        mean_vals = np.mean(np.mean(imarray, axis=0), axis=0)
        empty_vals = [np.argwhere(np.bincount(imarray[:,:,cval].flatten())==0) for cval in xrange(n_bands)]
        for cval,vals in enumerate(empty_vals):
            above_mean = vals[vals>mean_vals[cval]]
            if len(above_mean)>0:
                cutoff_vals[cval] = min(above_mean)
            else:
                cutoff_vals[cval] = np.amax(np.amax(imarray[:,:,cval], axis=0), axis=0) + 1

        max_cutoff = np.amax(cutoff_vals)

    im.close()
            

    print
    print "Calibrating S thresholds..."
    print "Processing levels..." 
    print " 0/%d" % n_img_bg

    """
    for i,im in enumerate(images):
        if (i+1) % 10 == 0:
            print " %d/%d" % (i+1,n_img_bg)
                        
        mean_grid += np.array(Image.open(im))

    mean_grid /= float(n_img_bg)

    print "Processing variances..."
    print " 0/%d" % n_img_bg
    
    for i,im in enumerate(images):
        if (i+1) % 10 == 0:
            print " %d/%d" % (i+1,n_img_bg)
                        
        var_grid += np.square(np.array(Image.open(im))-mean_grid)

    var_grid /= float(n_img_bg)

    s_grid = mean_grid + std_devs * np.sqrt(var_grid)

    """
    
    # find max of each pixel
    for i,im in enumerate(images):
        if (i+1) % 10 == 0:
            print " %d/%d" % (i+1,n_img_bg)

        s_grid = np.maximum(np.array(Image.open(im)), s_grid)

    # remove hot pixels and tracks
    if bg_cutoff:
        print "Removing thresholds above %d..." % bg_cutoff
          
        mask_kernel = np.array([[1,1,1,1,1],[1,0,0,0,1],[1,0,0,0,1],[1,0,0,0,1],[1,1,1,1,1]],dtype=float)/16.
        masked_grid = np.zeros((h, w, n_bands))
        for cval in xrange(n_bands):
            masked_grid[:,:,cval] = convolve2d(s_grid[:,:,cval], mask_kernel, mode='same', boundary='symm')
        s_grid = np.where(s_grid <= max_cutoff, s_grid, masked_grid)


    print "Downsampling image..."

    s_grid = np.amax([s_grid[x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
    
    print "Applying convolution kernel..."

    kernel_side = 2*conv_len+1
    s_kernel = np.repeat(1, kernel_side**2).reshape((kernel_side,kernel_side))/float(kernel_side)**2
    convolved_grid = np.array([convolve2d(s_grid[:,:,cval], s_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)]).transpose(1,2,0)
    s_grid = np.maximum(s_grid, convolved_grid)
    s_grid = np.ceil(s_grid+0.9).astype(int)

    # resize
    s_grid = np.repeat(np.repeat(s_grid, sample_block, axis=0), sample_block, axis=1)
                     
    s_img = Image.fromarray(s_grid.astype(np.uint8), mode=img_mode)

    # save as png
    img_name = out
    print "Saving background as %s" % img_name
    s_img.save(img_name)
    
    return s_grid.transpose(1,0,2)

if '__name__' == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='')
    parser.add_argument('--in', metavar='infiles', n_args='+', help='Images used to set thresholds')
    parser.add_argument('--out', default='bg.png', help='Output file name')
    parser.add_argument("--conv_len", type=int, default=0, help='Distance to which pixels are included in averaging')
    parser.add_argument("--bg_cutoff", action='store_true', help='Removes tracks during sauto processing.')
    parser.add_argument('--show', action='store_true', help='Display resulting threshold image')
    
    
    args = parser.parse_args()
    
    bg = find_bg(list(args.infiles), args.out, args.conv_len, args.bg_cutoff)
    if args.show:
        plt.figure(1)
        for cval in n_bands:
            plt.subplot(2,2,cval)
            plt.imshow(bg[:,:,cval], cmap='plasma', interpolation='nearest')
            plt.colorbar()
        plt.show()
      
        
    
 
