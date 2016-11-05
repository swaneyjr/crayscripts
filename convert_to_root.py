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

matplotlib.use('tkagg')

# given an ImGrid object and a threshold of pixels to keep, returns L2 values
def set_thresh(imarray, thresh):

    thresh_array = np.repeat(-1,n_bands)
    target_pix = thresh*imarray.size
    histarray = [np.bincount(imarray[cval].flatten()) for cval in xrange(imarray.n_bands)]

    for cval in xrange(n_bands):
        thresh_pix = target_pix
        while target_pix <= thresh_pix:
            thresh_array[cval] += 1
            thresh_pix = np.sum(histarray[cval][thresh_array[cval]:])

    return thresh_array

def convert_to_root(infiles, out, l1thresh=0, l2auto=True, l2manual=0, l2plus=0, sauto=True, smanual=False, max_img=0):
 
    avg3_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]])/8.0
    avg5_kernel = np.array([[1,1,1,1,1],[1,0,0,0,1],[1,0,0,0,1],[1,0,0,0,1],[1,1,1,1,1]])/16.0
    
    n_images_auto = 50
    
    saved_pix = 0
    total_pix = 0

    # handle S threshold settings
    if sauto:
        s_grid = find_bg(infiles[:n_images_auto], out.split('.')[0]+"_bg.png")
    elif smanual:
        s_grid = np.array(Image.open(smanual)).astype(int)
    else:
        s_grid = 0*imtools.ImGrid(infiles[0])
        
    if l2plus:
        s_grid += l2plus


    # create TTree
    print
    print "Creating TTree..."
    t = r.TTree("events", 'TTree with image data')

    n_img = np.array([0], dtype=int)
    pix_n = np.array([0], dtype=int)
    name = np.array('', dtype=str)
    color = np.array('', dtype=str)
   
    t.Branch('n_img', n_img, 'n_img/i')
    t.Branch('pix_n', pix_n, 'pix_n/i')
    t.Branch('col', color, 'col/C')
    t.Branch('name', name, 'name/C')   

    vbranch(t, 'pix_x', btype=int)
    vbranch(t, 'pix_y', btype=int)
    vbranch(t, 'pix_val', btype=int)
    vbranch(t, 'pix_avg3', btype=float)
    vbranch(t, 'pix_avg5', btype=float)
    vbranch(t, 'l2s', btype=int)


    # fill TTree
    print "Starting loop..."
    prev_name = ''
    for i,im_name in enumerate(infiles):
        
        print
        print "Image %d/%d:" % (i+1,len(infiles))

        imarray = imtools.ImGrid(im_name)
        n_bands = len(imarray.bands)
        im_split = im_name.split('.')

        im_pix = imarray.width*imarray.height

        # set L2 threshold
        if l2auto:
            if sauto or smanual:
                cx, cy = imarray.width/2, imarray.height/2
                c_side_half = int(imarray.height/2/math.sqrt(2))
                l2array = set_thresh(imarray[cx-c_side_half:cx+c_side_half,cy-c_side_half:cy+c_side_half], l2auto)
            else:
                l2array = set_thresh(imarray, l2auto)
            
                
        elif l2manual:
            l2array = np.repeat(max(l2manual),n_bands)
            for i,v in enumerate(l2manual):
                l2array[i] = v

        else:
            l2array = np.amin(np.amin(s_grid, axis=1), axis=1).astype(int)                   
            
        minl2thresh = np.amin(l2array)
        
            
        print "L2 threshold: \t",
        for cval, c in enumerate(imarray.bands):
            print "%s: %d" % (c, l2array[cval]),
        print

        # determine relevant S thresholds    
        grid_min = np.amin(np.amin(s_grid,axis=1),axis=1)
        l1diff = np.repeat(l1thresh,3)-grid_min
        l2diff = l2array-grid_min

        # enforce L1S
        if np.count_nonzero(imarray>=s_grid+l1diff) == 0: continue
        if im_split[0] != prev_name:
            n_img[0] += 1

        avg3_array = [convolve2d(imarray[cval], avg3_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)]
        avg5_array = [convolve2d(imarray[cval], avg5_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)]
            
        # fill TTree with image data for each band
        for cval, c in enumerate(imarray.bands):
            
            t.pix_x.clear()
            t.pix_y.clear()
            t.pix_val.clear()
            t.pix_avg3.clear()
            t.pix_avg5.clear()
            t.l2s.clear()
            
            for y,x in np.argwhere(imarray[cval] >= s_grid[cval]+l2diff[cval]):
                t.pix_x.push_back(x)
                t.pix_y.push_back(y)
                t.pix_val.push_back(imarray[cval][y,x])
                t.pix_avg3.push_back(avg3_array[cval][y,x])
                t.pix_avg5.push_back(avg5_array[cval][y,x])
                t.l2s.push_back(s_grid[cval][y,x]+l2diff[cval])
            
            color = np.array(c+'\0')
            name = np.array(im_split[0]+'\0')
            t.SetBranchAddress('col', color)
            t.SetBranchAddress('name', name)
            
            pix_n[0] = t.pix_x.size()
            print "%s: pix_n = %d" % (c, pix_n[0])
            saved_pix += pix_n[0]
            total_pix += im_pix
            
            t.Fill()
        
        prev_name = im_split[0]
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
    parser.add_argument("--l1", type=int, default=0, help='L1 threshold for highest band')
    
    l2option = parser.add_mutually_exclusive_group()
    l2option.add_argument("--l2auto", type=float, help='Target fraction of pixels kept after configuring L2 threshold')
    l2option.add_argument("--l2manual", nargs='+', type=int, help='L2 threshold for each band')
    l2option.add_argument("--l2plus", nargs='+', type=int, help='Increase the S thresholds by fixed amounts for each band.') 

    soption = parser.add_mutually_exclusive_group()
    soption.add_argument("--sauto", action='store_true', help='Add a spatially dependent threshold gradient.')
    soption.add_argument("--smanual", help='Manually add a threshold function. Input is an image of the same mode as those processed.')

    parser.add_argument("--max_img", type=int, help='Maximum number of images to convert')
    parser.add_argument("--show", action='store_true', help='Generate graphs of background thresholds and saved pixels')
    
    args = parser.parse_args()

    outfile = r.TFile(args.out, "recreate")        
        
    ti = time.clock()
    t = convert_to_root(args.infiles, args.out, args.l1, args.l2auto, args.l2manual, args.l2plus, args.sauto, args.smanual, args.max_img)

    tf = time.clock()
      
   
    outfile.Write()
    
    if args.show:
        
        print "Drawing saved pixels..."
        im = imtools.ImGrid(infiles[0])
        bands = im.bands
        n_bands = len(bands)
        c1 = r.TCanvas('c1','Saved Pixels',300,250*n_bands)
        c1.Divide(1,n_bands,0,0)
        for cval,c in enumerate(bands):
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
