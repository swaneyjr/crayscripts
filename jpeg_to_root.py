from PIL import Image
import ROOT as r
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import rawpy
import math
from fractions import gcd
import time
from scipy.signal import convolve2d
from hotcell import vbranch

matplotlib.use('tkagg')

# given an Image object and a threshold of pixels to keep, returns L2 values
def set_thresh(im, thresh):

    n_bands = len(im.mode)
    thresh_array = np.repeat(-1,n_bands)
    im_pix = im.height*im.width
    target_pix = thresh*im_pix
    histarray = np.array(im.histogram()).reshape((n_bands,256))

    for cval in xrange(n_bands):
        thresh_pix = target_pix
        while target_pix <= thresh_pix:
            thresh_array[cval] += 1
            thresh_pix = np.sum(histarray[cval][thresh_array[cval]:])

    return thresh_array

def convert_to_root(images, l1thresh, l2auto, l2manual, l2plus, sauto, smanual, border, max_img):

    saved_pix = 0
    total_pix = 0

    # handle S threshold settings
    if sauto:
        n_images_auto = 100
        s_grid = find_bg(images[:n_images_auto])
        images = images[n_images_auto:]
    elif smanual:
        s_grid = np.array(Image.open(smanual)).transpose(1,0,2).astype(int)
    else:
        im = Image.open(images[0])
        s_grid = np.zeros((im.width, im.height, len(im.mode)),dtype=int)


    # create TTree
    print
    print "Creating TTree..."
    t = r.TTree("events", 'TTree with JPEG data')

    n_img = np.array([0], dtype=int)
    phone_time = np.array([0], dtype=long)
    pix_n = np.array([0], dtype=int)
    color = np.array('', dtype=str)
    source = np.array([0], dtype=bool)
    l2 = np.array([0], dtype=int)
   
    t.Branch('n_img', n_img, 'n_img/i')
    t.Branch('t', phone_time, 't/i')
    t.Branch('pix_n', pix_n, 'pix_n/i')
    t.Branch('col', color, 'col/C')
    t.Branch('source', source, 'source/O')
    t.Branch('l2thresh', l2, 'l2thresh/i')
    

    vbranch(t, 'pix_x', btype=int)
    vbranch(t, 'pix_y', btype=int)
    vbranch(t, 'pix_val', btype=int)
    vbranch(t, 'pix_avg3', btype=float)
    vbranch(t, 'pix_avg5', btype=float)
    vbranch(t, 'l2s', btype=int)


    avg3_pix = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
    avg5_pix = [(-2,-2),(-2,-1),(-2,0),(-2,1),(-2,2),(-1,-2),(-1,2),(0,-2), \
                (0,2),(1,-2),(1,2),(2,-2),(2,-1),(2,0),(2,1),(2,2)]

    avg3_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]])/8.0
    avg5_kernel = np.array([[1,1,1,1,1,1,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1], \
                            [1,0,0,0,0,0,1],[1,0,0,0,0,0,1],[1,1,1,1,1,1,1]])/24.0


    # fill TTree
    print "Starting loop..."
    for i,jpg in enumerate(images):
        
        phone_time[0] = long(jpg[-20:-7])
        
        print
        print "Image %d/%d: t = %d" % (i+1,len(images),phone_time[0])

        im = Image.open(jpg)
        mode = im.mode
        n_bands = len(mode)
        
        # enforce L1 threshold
        min_band, max_band = zip(*list(im.getextrema()))
        if n_bands > 1 and max(max_band) < l1thresh: continue
        elif max_band < l1thresh: continue
        n_img[0] = i+1

        im_pix = im.width*im.height
        source[0] = bool(int(jpg[-5]))

        # set L2 threshold
        if l2auto:
            if sauto or smanual:
                cx, cy = im.width/2, im.height/2
                c_side_half = int(im.height/2/math.sqrt(2))
                cropped = im.crop((cx-c_side_half, cy-c_side_half, cx+c_side_half, cy+c_side_half))
                l2array = set_thresh(cropped, l2auto)
            else:
                l2array = set_thresh(im, l2auto)
            
                
        elif l2manual:
            l2array = np.repeat(max(l2manual),n_bands)
            for i,v in enumerate(l2manual):
                l2array[i] = v

        else:
            l2array = np.amin(np.amin(s_grid, axis=0), axis=0).astype(int)                   
            
        minl2thresh = min(l2array)
        
            
        print "L2 threshold: \t",
        for cval, c in enumerate(im.mode):
            print "%s: %d" % (c, l2array[cval]),
        print

        # determine relevant S thresholds
        l1s_satisfied = not (sauto or smanual)
        if sauto or smanual:
            grid_min = np.amin(np.amin(s_grid,axis=0),axis=0)
            l1diff = np.repeat(l1thresh,3)-grid_min
            l2diff = l2array-grid_min
                
                                
        x_vectors = [r.vector('int')() for i in xrange(n_bands)]
        y_vectors = [r.vector('int')() for i in xrange(n_bands)]
        val_vectors = [r.vector('int')() for i in xrange(n_bands)]
        avg3_vectors = [r.vector('double')() for i in xrange(n_bands)]
        avg5_vectors = [r.vector('double')() for i in xrange(n_bands)]
        l2s_vectors = [r.vector('int')() for i in xrange(n_bands)]
        
        # find val, avg3, avg5

        if not args.numpy:
        
            for j,val in enumerate(im.getdata()):

                # enforce L2 threshold
                if n_bands > 1 and max(val) < minl2thresh: continue
                elif n_bands == 1 and val < minl2thresh: continue
                elif n_bands > 1 and max(val-l2array) < 0: continue

                x,y = j%im.width, j/im.width

                if x < border or x >= im.width-border or y < border or \
                   y >= im.height-border: continue

                #if args.sample16 and (x%4 != 0 or y%4 != 0): continue 

                # determine L1S and L2S
                if sauto or smanual:
                    val_minus_l2s = val-s_grid[x,y]-l2diff
                    if max(val_minus_l2s) < 0: continue
                    if max(val-s_grid[x,y]-l1diff) >= 0:
                        l1s_satisfied = True


                # find avg3 and avg5
                        
                sum3 = np.zeros(n_bands)
                sum5 = np.zeros(n_bands)

                n3 = 0
                n5 = 0
                    
                for dx,dy in avg3_pix:
                    xprime = x+dx
                    yprime = y+dy
                    if xprime >= 0 and xprime < im.width and yprime >= 0 and yprime < im.height:
                        n3 += 1
                        sum3 += im.getpixel((xprime,yprime))
                    
                for dx,dy in avg5_pix:
                    xprime = x+dx
                    yprime = y+dy
                    if xprime >= 0 and xprime < im.width and yprime >= 0 and yprime < im.height:
                        n5 += 1
                        sum5 += im.getpixel((xprime,yprime))
                                
                
                avg3 = np.divide(sum3,n3)
                avg5 = np.divide(sum5,n5)

                for cval in xrange(n_bands):
                    if val_minus_l2s[cval] >= 0:
                        x_vectors[cval].push_back(x)
                        y_vectors[cval].push_back(y)
                        val_vectors[cval].push_back(val[cval])
                        avg3_vectors[cval].push_back(avg3[cval])
                        avg5_vectors[cval].push_back(avg5[cval])
                        l2s_vectors[cval].push_back(s_grid[x,y,cval]+l2diff[cval])
               
            im.close()       
            if not l1s_satisfied: continue

        # using numpy methods    
        else:
            imarray = np.array(im, dtype=int).transpose(1,0,2)
            im.close()

            # enforce L1S
            if np.count_nonzero(imarray>=s_grid+l1diff) == 0: continue

            avg3_array = [convolve2d(imarray[:,:,cval], avg3_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)]
            avg5_array = [convolve2d(imarray[:,:,cval], avg5_kernel, mode='same', boundary='symm') for cval in xrange(n_bands)]

            
            for x,y,cval in np.argwhere(imarray>=s_grid+l2diff):
                x_vectors[cval].push_back(x)
                y_vectors[cval].push_back(y)
                val_vectors[cval].push_back(imarray[x,y,cval])
                avg3_vectors[cval].push_back(avg3_array[cval][x,y])
                avg5_vectors[cval].push_back(avg5_array[cval][x,y])
                l2s_vectors[cval].push_back(s_grid[x,y,cval]+l2diff[cval])
                

        # fill TTree with image data for each band
        for cval, c in enumerate(im.mode):
            
            color = np.array(c+'\0')
            t.SetBranchAddress('col', color)
            l2[0] = l2array[cval]

            setattr(t,'pix_x', x_vectors[cval])
            setattr(t,'pix_y', y_vectors[cval])
            setattr(t,'pix_val', val_vectors[cval])
            setattr(t,'pix_avg3', avg3_vectors[cval])
            setattr(t,'pix_avg5', avg5_vectors[cval])
            setattr(t, 'l2s', l2s_vectors[cval])

            t.SetBranchAddress('pix_x', x_vectors[cval])
            t.SetBranchAddress('pix_y', y_vectors[cval])
            t.SetBranchAddress('pix_val', val_vectors[cval])
            t.SetBranchAddress('pix_avg3', avg3_vectors[cval])
            t.SetBranchAddress('pix_avg5', avg5_vectors[cval])
            t.SetBranchAddress('l2s', l2s_vectors[cval])
            
            pix_n[0] = t.pix_x.size()
            print "%s: pix_n = %d" % (c, pix_n[0])
            saved_pix += pix_n[0]
            total_pix += im_pix
            
            t.Fill()
            
        if max_img and n_img >= max_img: break

    print "Done!"
    print
    print "Images saved:    %d\t(%.2f%%)" % (n_img[0], 100.*n_img[0]/float(len(images)))
    print "Total pixels:    %d\t(%.2f%%)" % (saved_pix, 100.*saved_pix/total_pix)

    if args.show:
        
        print "Drawing saved pixels..."
        c1 = r.TCanvas('c1','Saved Pixels',300,250*n_bands)
        c1.Divide(1,n_bands,0,0)
        for cval,c in enumerate(mode):
            c1.cd(cval+1)
            t.Draw('pix_y:pix_x','col=="%s"' % c, 'colz')

        print "Drawing background..."
        
        plt.figure(1)
        d = math.ceil(math.sqrt(n_bands))
        for b in xrange(s_grid.shape[2]):
            plt.subplot(d,d,b+1)
            bg = plt.imshow(s_grid[:,:,b].transpose(1,0), cmap='plasma',interpolation='nearest', vmin=0, vmax=20)
            plt.colorbar()

        plt.show()            
    
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
    #parser.add_argument("--devs", type=float, default=0, help='Number of standard deviations above mean to set sauto')
    parser.add_argument("--conv_len", type=int, default=0, help='Distance to which pixels are included in sauto averaging')
    parser.add_argument("--bg_cutoff", action='store_true', help='Removes tracks during sauto processing.')

    parser.add_argument("--border", type=int, default=0, help='Number of pixels around edge to crop out')
    parser.add_argument("--max_img", type=int, help='Maximum number of images to convert')
    parser.add_argument("--show", action='store_true', help='Generate graphs of background thresholds and saved pixels')
    parser.add_argument("--numpy", action='store_true', help='Apply thresholds with numpy arrays.  Useful for large threshold gradients')

    args = parser.parse_args()

    images = list(args.infiles)

    outfile = r.TFile(args.out, "recreate")        
        
    ti = time.clock()
    t = convert_to_root(images, args.l1, args.l2auto, args.l2manual, args.l2plus, args.sauto, args.smanual, \
                        args.border, args.max_img)

    tf = time.clock()
      
   
    outfile.Write()
    outfile.Close()

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



