import rawpy
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from time import clock
import gzip
from hotcell import vbranch

matplotlib.use('tkagg')

def convert_to_root(images, l1thresh, l2thresh):
    
    saved_pix = 0
    total_pix = 0
    adc_hist = np.zeros(1024)

    avg3_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]])/8.0
    avg5_kernel = np.array([[1,1,1,1,1,1,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1],[1,0,0,0,0,0,1], \
                            [1,0,0,0,0,0,1],[1,0,0,0,0,0,1],[1,1,1,1,1,1,1]])/24.0


    # create TTree
    print
    print "Creating TTree..."
    t = r.TTree("events", 'TTree with JPEG data')

    n_img = np.array([0], dtype=int)
    pix_n = np.array([0], dtype=int)
    l2 = np.array([0], dtype=int)
   
    t.Branch('n_img', n_img, 'n_img/i')
    t.Branch('pix_n', pix_n, 'pix_n/i')
    t.Branch('l2thresh', l2, 'l2thresh/i')

    vbranch(t, 'pix_x', btype=int)
    vbranch(t, 'pix_y', btype=int)
    vbranch(t, 'pix_val', btype=int)
    vbranch(t, 'pix_avg3', btype=float)
    vbranch(t, 'pix_avg5', btype=float)     

    # fill TTree
    for i,z in enumerate(images):

        print
        print "Image %d/%d:" % (i+1,len(images))

        # open file
        with gzip.open(z) as r:
            with rawpy.imread(r) as raw:
                
                raw_im = raw.raw_image()
                adc_hist += np.bincount(raw_im, minlength=1024)

                if not args.out: continue
                
                if np.amax(raw_im) < l1thresh: continue
                n_img[0] += 1

                avg3_array = convolve2d(raw_im[:,:], avg3_kernel, mode='same', boundary='symm')
                avg5_array = convolve2d(raw_im[:,:], avg5_kernel, mode='same', boundary='symm')

                for x,y in np.argwhere(raw_im >= l2thresh):
                    t.pix_x.push_back(x)
                    t.pix_y.push_back(y)
                    t.pix_val.push_back(raw_im[x,y])
                    t.pix_avg3.push_back(avg3_array[x,y])
                    t.pix_avg5.push_back(avg5_array[x,y])

                pix_n[0] = t.pix_x.size()
                print "pix_n = %d" % pix_n[0]
                saved_pix += pix_n[0]
                total_pix += raw_im.size

        t.Fill()

        if max_img and n_img >= max_img: break

        print "Done!"
        print
        print "Images saved:    %d\t(%.2f%%)" % (n_img[0], 100.*n_img[0]/float(len(images)))
        print "Total pixels:    %d\t(%.2f%%)" % (saved_pix, 100.*saved_pix/total_pix)

        return t, adc_hist
                
          
    

if name == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument('--in', required=True, dest='infiles', nargs='+', help='Raw images to be processed')
    parser.add_argument('--out', help='Create a ROOT file with pixels above L2 threshold')
    parser.add_argument('--l1', type=int, default=0, help='L1 threshold')
    parser.add_argument('--l2', type=int, default=0, help='L2 threshold')
    parser.add_argument('--show', action='store_true', help='Draw all figures after processing')

    
    args = parser.parse_args()

    images = list(args.infiles)

    if args.out:
        outfile = r.TFile(args.out, "recreate")

    ti = clock()
    t, adc_hist = convert_to_root(images, args.l1, args.l2)
    tf = clock()

    if args.out:
        outfile.Write()
        outfile.Close()

    if args.show:

        print "Drawing saved pixels..."
        c1 = r.TCanvas('c1','Saved Pixels',300,250)
        t.Draw('pix_y:pix_x','', 'colz')

        print "Drawing ADC counts..."
        plt.figure(1)
        plt.hist(np.arange(1024), weights=adc_hist)
        plt.xlabel('ADC Count')
        plt.ylabel('Frequency')
        plt.show()
    
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
