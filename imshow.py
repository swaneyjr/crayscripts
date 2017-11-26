#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from scipy.signal import convolve2d
import imtools

# enable X-forwarding
matplotlib.use('tkagg')

def plotchannel(imarray, cval, minus_bg=False, div_bg=False, colormap='plasma'):
    if minus_bg:
        mx = 10
    elif div_bg:
        mx = 5
    else:
        mx = imtools.outlier_cutoff(imarray, imarray.size/10000)[cval]
        mx += 5 - (mx%5)
    plt.imshow(imarray[cval], cmap=colormap, interpolation='nearest',vmin=0, vmax=mx)
    plt.colorbar()
    
def imshow(img, minus_bg=False, div_bg=False, sample=False, cmap='plasma'):

    imarray = imtools.ImGrid(img)
    sample_block = 4

    if minus_bg or div_bg:
        bg = np.median([imarray[x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
        bg = np.repeat(np.repeat(bg, sample_block, axis=0), sample_block, axis=1)
    if minus_bg:
        imarray -= bg.astype(int)
    elif div_bg:
        imarray = imarray.astype(float) / bg

    if sample:
        imarray = np.amax([imarray[:,x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
        imarray = np.repeat(np.repeat(imarray, sample_block, axis=1), sample_block, axis=2)


    for cval,b in enumerate(imarray.bands):
        plt.figure(cval+1)
        plotchannel(imarray, cval, minus_bg, div_bg, cmap)
        plt.title(b)

    plt.show()

    
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'Make a colormap out of an image')
    parser.add_argument('img', help='Image to be displayed')
    bg_option = parser.add_mutually_exclusive_group()
    bg_option.add_argument('-m','--minus_bg', action='store_true', help='Subtract avg7 values at each pixel')
    bg_option.add_argument('-d','--div_bg', action='store_true', help='Divide by avg7 values at each pixel.')
    parser.add_argument('-s',"--sample", action='store_true', help='Downsamples to max value in 4x4 block.')
    parser.add_argument("--cmap", default='plasma', help='Matplotlib color map to be used.')
    args = parser.parse_args()
    
    imshow(args.img, args.minus_bg, args.div_bg, args.sample, args.cmap)


