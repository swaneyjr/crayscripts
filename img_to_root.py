from PIL import Image
import ROOT as r
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import time
from scipy.signal import convolve2d
import imtools
from hotcell import vbranch
from bg_threshold import find_bg

# find appropriate L1 thresholds for a given target rate
def find_l1(imlist, l1_target_rate, bg_grid, dev_grid, b=0):
    
    # we don't need to survey the whole video
    if len(imlist) > 10/l1_target_rate:
        imlist = imlist[:int(10/l1_target_rate)]
    target_saved = int(l1_target_rate*len(imlist))
    max_vals = np.zeros((len(imlist), bg_grid.shape[0]))
    for i,im in enumerate(imlist):
        max_vals[i] = np.amax(np.amax((imtools.ImGrid(im)[:,b:-b,b:-b]-bg_grid)/dev_grid, axis=1), axis=1)
    l1array = np.sort(max_vals,axis=0)[-target_saved]
    
    return l1array

# given an ImGrid object and a threshold of pixels to keep, returns L2 values
def set_L2_thresh(imarray, thresh):

    thresh_array = np.repeat(-1, imarray.n_bands)
    target_pix = thresh*imarray.size
    histarray = [np.bincount(imarray[cval].flatten()) for cval in xrange(imarray.n_bands)]

    for cval in xrange(imarray.n_bands):
        thresh_pix = target_pix
        while target_pix <= thresh_pix:
            thresh_array[cval] += 1
            thresh_pix = np.sum(histarray[cval][thresh_array[cval]:])

    return thresh_array

def convert_to_root(infiles, l1_target_rate=None, l2auto=0, l2manual=0, s_thresh=True, border=0, max_img=0, rawcam_format=False):
 
    avg3_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]])/8.0
    avg5_kernel = np.array([[1,1,1,1,1],[1,0,0,0,1],[1,0,0,0,1],[1,0,0,0,1],[1,1,1,1,1]])/16.0
    
    n_images_auto = 50
    
    saved_pix = 0
    total_pix = 0

    # handle S threshold settings
    if s_thresh:
        bg_grid = find_bg(infiles[:n_images_auto])
        if border: bg_grid = bg_grid[:,border:-border,border:-border]
        infiles = infiles[n_images_auto:]
        dev_grid = np.sqrt(bg_grid)
    else:
        bg_grid = np.zeros(imtools.ImGrid(infiles[0]).shape)
        if border: bg_grid = bg_grid[:,border:-border,border:-border]
        dev_grid = np.ones(bg_grid.shape)
        


    # create TTree
    print
    print "Creating TTree..."
    t = r.TTree("events", 'TTree with image data')

    n_img = np.array([0], dtype=int)
    pix_n = np.array([0], dtype=int)
    bg_lvl = np.array([0], dtype=float)
    name = np.array('', dtype=str)
    color = np.array('', dtype=str)
   
    t.Branch('n_img', n_img, 'n_img/i')
    t.Branch('pix_n', pix_n, 'pix_n/i')
    t.Branch('bg_lvl', bg_lvl, 'bg_lvl/D')
    t.Branch('col', color, 'col/C')
    t.Branch('name', name, 'name/C') 
    
    if rawcam_format:
        time = np.array([0], dtype=long)
        t.Branch('t', time, 't/L')
        im_split = (infiles[0].split('.')[0]).split('_')
        if len(im_split) > 2:
            brancharray = np.zeros((len(im_split)-2,1)).astype(int)
            for i,b in enumerate(im_split[2:]):
                t.Branch(b[0], brancharray[i], b[0]+'/i')
        

    vbranch(t, 'pix_x', btype=int)
    vbranch(t, 'pix_y', btype=int)
    vbranch(t, 'pix_val', btype=int)
    vbranch(t, 'pix_avg3', btype=float)
    vbranch(t, 'pix_avg5', btype=float)
    vbranch(t, 'pix_bg', btype=float)

    print "Finding L1 Threshold..."

    if l1_target_rate:
        l1array = find_l1(infiles, l1_target_rate, bg_grid, dev_grid)
    else:
        l1array = np.zeros(bg_grid.shape[0])
    
    l1_grid = l1array.reshape(l1array.size,1,1) * dev_grid + bg_grid

    # fill TTree
    print "Starting loop..."
    prev_name = ''
    for i,im_name in enumerate(infiles):

        imarray = imtools.ImGrid(im_name)
        if border:
            imarray = imarray[:,border:-border,border:-border]
        im_base = im_name.split('.')

        # enforce L1S
        if np.count_nonzero(imarray >= l1_grid) == 0: continue
        if im_base[0] != prev_name:
            n_img[0] += 1

        print
        print "File %d/%d:" % (i+1,len(infiles))

        # set L2 threshold
        if l2auto:
            if s_thresh:
                cx, cy = imarray.width/2, imarray.height/2
                c_side_half = int(imarray.height/2/math.sqrt(2))
                l2array = set_L2_thresh(imarray[cx-c_side_half:cx+c_side_half,cy-c_side_half:cy+c_side_half], l2auto)
            else:
                l2array = set_L2_thresh(imarray, l2auto)
            
                
        elif l2manual:
            l2array = np.repeat(max(l2manual),imarray.n_bands)
            for i,v in enumerate(l2manual):
                l2array[i] = v

        elif s_thresh:
            l2array = 0.9*l1array

        else:
            l2array = l1array-1

        l2_grid = l2array.reshape(imarray.n_bands,1,1)*dev_grid + bg_grid    
        
            
        print "L2 threshold: \t",
        for cval, c in enumerate(imarray.bands):
            print "%s: %1.1f" % (c, l2array[cval]),
        print

        avg3_array = [convolve2d(imarray[cval], avg3_kernel, mode='same', boundary='symm') for cval in xrange(imarray.n_bands)]
        avg5_array = [convolve2d(imarray[cval], avg5_kernel, mode='same', boundary='symm') for cval in xrange(imarray.n_bands)]
            
        # fill TTree with image data for each band
        for cval, c in enumerate(imarray.bands):
            
            t.pix_x.clear()
            t.pix_y.clear()
            t.pix_val.clear()
            t.pix_avg3.clear()
            t.pix_avg5.clear()
            t.pix_bg.clear()
            
            for y,x in np.argwhere(imarray[cval] >= l2_grid[cval]):
                t.pix_x.push_back(x)
                t.pix_y.push_back(y)
                t.pix_val.push_back(imarray[cval][y,x])
                t.pix_avg3.push_back(avg3_array[cval][y,x])
                t.pix_avg5.push_back(avg5_array[cval][y,x])
                t.pix_bg.push_back(bg_grid[cval][y,x])
            
            color = np.array(c+'\0')
            name = np.array(im_base[0]+'\0')
            t.SetBranchAddress('col', color)
            t.SetBranchAddress('name', name)
            
            if rawcam_format:
                im_split = im_base[0].split('_')
                time[0] = int(im_split[1])
                if len(im_split) > 2:
                    for i,b in enumerate(im_split[2:]):
                        brancharray[i][0] = int(b[1:])
                        t.SetBranchAddress(b[0], brancharray[i])
            
            pix_n[0] = t.pix_x.size()
            bg_lvl[0] = np.mean(imarray[cval])
            print "%s: pix_n = %d" % (c, pix_n[0])
            saved_pix += pix_n[0]
            
            t.Fill()
        
        total_pix += imarray.size
        prev_name = im_base[0]
        if max_img and n_img >= max_img: break

    print "Done!"
    print
    print "Images saved:    %d\t(%.2f%%)" % (n_img[0], 100.*n_img[0]/float(len(infiles)))
    print "Total pixels:    %d\t(%.2f%%)" % (saved_pix, 100.*saved_pix/total_pix)           
    
    return t
    
    

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'Converts JPEG files to a ROOT file')
    parser.add_argument("--in", required=True, dest='infiles', nargs='+', help='Images to be converted to a ROOT file')
    parser.add_argument("--out", default='images.root', help='Output file name')
    parser.add_argument("--l1_rate", type=float, default=None, help='L1 threshold for highest band')
    
    l2option = parser.add_mutually_exclusive_group()
    l2option.add_argument("--l2auto", type=float, help='Target fraction of pixels kept after configuring L2 threshold')
    l2option.add_argument("--l2manual", nargs='+', type=int, help='L2 threshold for each band')

    parser.add_argument("--border", type=int, help="Number of pixels to discard around edges")
    parser.add_argument("--max_img", type=int, help='Maximum number of images to convert')
    parser.add_argument('-r', "--rawcam_format", action='store_true', help='Creates branches from image name: name_t_(bname val)_(bname val)...')
    parser.add_argument('-s', "--show", action='store_true', help='Generate graphs of background thresholds and saved pixels')
    parser.add_argument('-S',"--s_thresh", action='store_true', help='Add a spatially dependent threshold gradient.')

    args = parser.parse_args()

    outfile = r.TFile(args.out, "recreate")        
        
    ti = time.clock()
    t = convert_to_root(args.infiles, args.l1_rate, args.l2auto, args.l2manual,
                        args.s_thresh, args.border, args.max_img, args.rawcam_format)

    tf = time.clock()
      
   
    outfile.Write()
    
    if args.show:

        matplotlib.use('tkagg')                           
        
        print "Drawing saved pixels..."
        im = imtools.ImGrid(args.infiles[0])
        c1 = r.TCanvas('c1','Saved Pixels',300,250*im.n_bands)
        c1.Divide(1,im.n_bands,0,0)
        for cval,c in enumerate(im.bands):
            c1.cd(cval+1)
            t.Draw('pix_y:pix_x','col=="%s"' % c, 'colz')

    m,s = divmod(tf-ti,60)
    h,m = divmod(m,60)
    
    print "Done! Wrote to %s." % args.out
    print "Total time: ",
    if tf-ti > 3600:
        print "%d h %02d m %02d s" % (h,m,s)
    elif tf-ti > 60:
        print "%d m %02d s" % (m,s)
    else:
        print "%f s" % (tf-ti)
        
    outfile.Close()
