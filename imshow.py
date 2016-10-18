import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from scipy.signal import convolve2d
from bg_threshold import outlier_cutoff

# enable X-forwarding
matplotlib.use('tkagg')

from argparse import ArgumentParser
parser = ArgumentParser(description = 'Make a colormap out of an image')
parser.add_argument('img', help='Image to be displayed')
parser.add_argument('-r', "--red", action='store_true', help='Draw R channel')
parser.add_argument('-g', "--green", action='store_true', help='Draw G channel')
parser.add_argument('-b', "--blue", action='store_true', help='Draw B channel')
bg_option = parser.add_mutually_exclusive_group()
bg_option.add_argument('-m','--minus_bg', action='store_true', help='Subtract avg7 values at each pixel')
bg_option.add_argument('-d','--div_bg', action='store_true', help='Divide by avg7 values at each pixel.')
parser.add_argument('-s',"--sample", action='store_true', help='Downsamples to max value in 8x8 block.')
args = parser.parse_args()


def plotchannel(imarray, cval):
    if args.minus_bg:
        mx = 10
    elif args.div_bg:
        mx = 5
    else:
        mx = outlier_cutoff(imarray)
        mx += 5 - (mx%5)
    plt.imshow(imarray[:,:,cval], cmap='plasma', interpolation='nearest',vmin=0, vmax=mx)
    plt.colorbar()
    
imarray = np.array(Image.open(args.img)).astype(int)
sample_block = 8

if args.minus_bg or args.div_bg:
    bg = np.median([imarray[x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
    bg = np.repeat(np.repeat(bg, sample_block, axis=0), sample_block, axis=1)
if args.minus_bg:
    imarray -= bg.astype(int)
elif args.div_bg:
    imarray = imarray.astype(float) / bg

if args.sample:
    imarray = np.amax([imarray[x::sample_block,y::sample_block] for x,y in np.ndindex(sample_block,sample_block)], axis=0)
    imarray = np.repeat(np.repeat(imarray, sample_block, axis=0), sample_block, axis=1)

    
n_fig = 1

if args.red:
    plt.figure(n_fig)
    plotchannel(imarray, 0)
    n_fig += 1
if args.green:
    plt.figure(n_fig)
    plotchannel(imarray, 1)
    n_fig += 1
if args.blue:
    plt.figure(n_fig)
    plotchannel(imarray, 2)

plt.show()



