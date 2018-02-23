#!/usr/bin/env python2

import ROOT as r
import numpy as np
from itertools import starmap, izip
from collections import namedtuple
from operator import itemgetter
import sys
import math
import cmath
from hotcell import vbranch

Pixel = namedtuple('Pixel',['x','y','val', 'adjusted_val', 'avg_3','avg_5'])
Clean_Pixel = namedtuple('Clean_Pixel',['x','y','val', 'adjusted_val', 'avg_3','avg_5','freq','masked'])

# finds the furthest two pixels in a list of pixels (or one such pair)
def find_endpoints(pixels):
    
    if len(pixels) == 1:
        return pixels[0],pixels[0]
    
    furthest_pix = {p: 0 for p in pixels}

    for i,p1 in enumerate(pixels):
        for p2 in pixels[i+1:]: 
            dist_sq = (p1.x-p2.x)**2 + (p1.y-p2.y)**2
            furthest_pix[p1] = max(furthest_pix[p1], dist_sq)
            furthest_pix[p2] = max(furthest_pix[p2], dist_sq)

    max_dist_sq = max(furthest_pix.itervalues())
    possible_endpoints = [p for p in pixels if furthest_pix[p] == max_dist_sq]
    start = possible_endpoints[0]
    for i in range(1,len(possible_endpoints)):
        end = possible_endpoints[i]
        if (start.x-end.x)**2+(start.y-end.y)**2 == max_dist_sq:
            break

    return start,end
    

def cluster(block, iso_thresh):
    
    def d(p1,p2):
        return math.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2)
    
    clusters = []
    while block != []:
        p = block.pop()
        match = None
        for cl in clusters:
            for pcl in cl:
                if d(p,pcl) < iso_thresh:
                    if match:
                        clusters.remove(match)
                        cl += match
                        match = cl
                        break
                    else:
                        cl.append(p)
                        match = cl
                        break
        if match == None:
            clusters.append([p])
                               
    return clusters
    
    
# returns a list of tracks separated by > iso_thresh and whose pixels are
# within iso_thresh of each other
def extract_tracks(pixels, min_pix_tr, iso_thresh, use_cluster=False):

    if len(pixels) == 1 and min_pix_tr <= 1:
        return [pixels]
    
    x_track_list = []

    # separate into groups by x
    pixels = sorted(pixels, key=itemgetter(0))
    pix_x = np.array(zip(*pixels)[0])
    x_diffs = np.diff(pix_x)
    x_track = [pixels[0]]
    for i,p in enumerate(pixels[1:]):
        if x_diffs[i] < iso_thresh:
            x_track.append(p)
        else:
            if len(x_track) >= min_pix_tr:
                x_track_list.append(x_track)
            x_track = [p]
    if len(x_track) >= min_pix_tr:
        x_track_list.append(x_track)

    xy_track_list = []
    # further separate each set of tracks by y
    for x_track in x_track_list:
        x_track = sorted(x_track, key=itemgetter(1))
        pix_y = np.array(zip(*x_track)[1])
        y_diffs = np.diff(pix_y)
        y_track = [x_track[0]]
        for i,p in enumerate(x_track[1:]):
            if y_diffs[i] < iso_thresh:
                y_track.append(p)
            else:
                if len(y_track) >= min_pix_tr:
                    xy_track_list.append(y_track)
                y_track = [p]
        if len(y_track) >= min_pix_tr:
            xy_track_list.append(y_track)
   
    tracks = []
    for block in xy_track_list:
        if len(block) < 3 or not use_cluster:
            tracks.append(block)
        else:
            tracks += cluster(block, iso_thresh)
        
    return tracks
    
    
# calculates the angle of a track via hough transform
def find_hough_theta(pixels):

    n = len(pixels)
    
    # orthogonal least squares fit
    sum_x = sum_y = sum_x_sq = sum_y_sq = sum_xy = sum_weights = 0
    for p in pixels:
        weight = p.val
        sum_x += weight*p.x
        sum_y += weight*p.y
        sum_x_sq += weight*p.x**2
        sum_y_sq += weight*p.y**2
        sum_xy += weight*p.x*p.y
        sum_weights += weight

    num = sum_weights*sum_xy-sum_x*sum_y
    den = sum_weights*sum_x_sq-sum_weights*sum_y_sq-sum_x**2+sum_y**2
    if den == 0:
        if num == 0:
            theta = float('nan')
        else:
            theta = math.copysign(math.pi/4., -num)
    else:
        theta = 0.5 * math.atan2(-2*num, -den)
        if theta == math.pi/2:
            theta *= -1
    rho = (sum_x*math.cos(theta)+sum_y*math.sin(theta))/sum_weights

    # calculate other stats based on theta fit

    rho_vals = [p.x*math.cos(theta)+p.y*math.sin(theta) for p in pixels]
    sigma_vals = [p.x*math.sin(theta)-p.y*math.cos(theta) for p in pixels]
    sig_min, sig_max = min(sigma_vals), max(sigma_vals)
    
    if n > 2:
        rho_std = np.std(rho_vals, ddof = 2)
        
    else:
        rho_std = float('nan')

    if n > 3:
        sig_plus = []
        sig_minus = []
        for i in xrange(n):
            if rho_vals[i] > rho:
                sig_plus.append(sigma_vals[i])
            else:
                sig_minus.append(sigma_vals[i])
        var_plus = np.var(sig_plus)
        var_minus = np.var(sig_minus)
        curv_ratio = min(var_plus, var_minus)/max(var_plus, var_minus)
    else:
        curv_ratio = None
    
    return theta, rho_std, rho, sig_min, sig_max, curv_ratio

def clean_tracks(t0, min_pix_tr=2, iso_thresh=3, fit=False, sel=None, cluster=False):
    saved_pix = total_pix = 0
    total_events = t0.GetEntries()

    print "Starting CopyTree..."
    if sel:
        sel += " && pix_n >= " + str(min_pix_tr)
    else:
        sel = "pix_n >= " + str(min_pix_tr)
    t1 = t0.CopyTree(sel)
    t1_events = t1.GetEntries()
    print "Pixels over min_pix_tr: %d/%d %.0f%%" % (t1_events, total_events, 100.*t1_events/total_events)
    print "Starting CloneTree..."
    t2 = t0.CloneTree(0)
    
    cleaned = 'pix_freq' in t2.GetListOfBranches()
             
    # add branches

    print "Adding TBranches..."
    pix_tr = np.array([0], dtype=int)
    
    max_val = np.array([0], dtype=int)
    max_x = np.array([0], dtype=int)
    max_y = np.array([0], dtype=int)
    tr_len = np.array([0], dtype=float)
    discr_len = np.array([0], dtype=int)
    tr_eff = np.array([0], dtype=float)
    if fit:
        hough_theta = np.array([0], dtype=float)
        rho_std = np.array([0], dtype=float)
        curv_ratio = np.array([0], dtype=float)

    t2.Branch('max_val', max_val, 'max_val/i')
    t2.Branch('max_x', max_x, 'max_x/i')
    t2.Branch('max_y', max_y, 'max_y/i')
    t2.Branch('tr_len', tr_len, 'tr_len/d')
    t2.Branch('discr_len', discr_len, 'discr_len/i')
    t2.Branch('tr_eff', tr_eff, 'tr_eff/d')
    t2.Branch('pix_tr', pix_tr, 'pix_tr/i')
    if fit:
        t2.Branch('theta', hough_theta, 'theta/d')
        t2.Branch('rho_std', rho_std, 'rho_std/d')
        t2.Branch('curv_ratio', curv_ratio, 'curve_ratio/d')
        vbranch(t2, 'd_rho')
        vbranch(t2, 'sigma')
    
    print 'Starting loop...'
    sys.stdout.flush()
    for ievt,evt in enumerate(t1):
        if ievt%(t1_events/10) == 0:
            print "%d/%d  %.1f%%" % (ievt, t1_events, 100.*ievt/t1_events)

        # read pixels from TTree
        if cleaned:
            raw_pixels = list(starmap(Clean_Pixel, izip(evt.pix_x, evt.pix_y, evt.pix_val, evt.pix_avg3, evt.pix_avg5, \
                                                        evt.pix_freq, evt.pix_masked)))
        else:    
            raw_pixels = list(starmap(Pixel, izip(evt.pix_x, evt.pix_y, evt.pix_val, evt.pix_avg3, evt.pix_avg5)))
        total_pix += len(raw_pixels)
        tracks = extract_tracks(raw_pixels, min_pix_tr, iso_thresh, use_cluster=cluster)

        
        for cleaned_pixels in tracks:
              
            # clear pixel data to rewrite in new TTree
            t2.pix_x.clear()
            t2.pix_y.clear()
            t2.pix_val.clear()
            t2.pix_avg3.clear()
            t2.pix_avg5.clear()

            if fit:
                t2.sigma.clear()
                t2.d_rho.clear()
                
            max_val[0] = 0
            n_max = 1

            # find new variables of interest
            if fit:
                hough_theta[0], rho_std[0], rho, sig_min, sig_max, curv_ratio[0] = find_hough_theta(cleaned_pixels)
            start, end = find_endpoints(cleaned_pixels)
            tr_len[0] = np.sqrt([(start.x-end.x)**2 + (start.y-end.y)**2])
            discr_len[0] = abs(start.x-end.x) + abs(start.y-end.y) + 1
            if discr_len > 2:
                tr_eff[0] = (len(cleaned_pixels)-2.)/(discr_len-2.)
            else:
                tr_eff[0] = None

            for p in cleaned_pixels:   

                # append pixels to new TTree data
                t2.pix_x.push_back(p.x)
                t2.pix_y.push_back(p.y)
                t2.pix_val.push_back(p.val)
                t2.pix_avg3.push_back(p.avg3)
                t2.pix_avg5.push_back(p.avg5)
                if cleaned:
                    t2.pix_freq.push_back(p.freq)
                    t2.pix_freq.push_back(p.masked)

                if fit:
                    t2.sigma.push_back((p.x*math.sin(hough_theta[0])-p.y*math.cos(hough_theta[0])-sig_min)/(sig_max-sig_min))
                    t2.d_rho.push_back(p.x*math.cos(hough_theta[0])+p.y*math.sin(hough_theta[0])-rho)
                

                if p.val > max_val[0]:
                    # replace with new maximum and reset counter
                    max_val[0] = p.val
                    max_x[0] = p.x
                    max_y[0] = p.y
                    n_max = 1
                    
                elif p.val == max_val[0]:
                    # average xy of maximum points
                    max_x[0] += p.x
                    max_y[0] += p.y
                    n_max += 1

            
            pix_tr[0] = t2.pix_x.size()
            saved_pix += t2.pix_x.size()
            
            max_x[0] /= n_max
            max_y[0] /= n_max

            t2.Fill()

    
    print "Done!"
    print
    print "Input events:    %d" % total_events
    print "Output events:   %d" % t2.GetEntries()
    print "Pixels in raw tracks:  %d" % total_pix
    print "Pixels saved:    %d\t(%.2f%%)" % (saved_pix, 100.*saved_pix/total_pix)

    return t2


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'Extract clean tracks from data')
    parser.add_argument("--in", required=True, dest='infiles', nargs='+', help='File to be cleaned')
    parser.add_argument("--out", default='tracks.root', help='Output file name')
    parser.add_argument("--min", type=int, default=2, help='Minimum pix_tr to be kept')
    parser.add_argument("--iso", type=int, default=5, help='Distance above which a pixel is classified as isolated')
    parser.add_argument('-f', "--fit", action='store_true', help='Compute linear fit parameters')
    parser.add_argument('-c', "--cluster", action='store_true', help='Use more advanced clustering algorithm')
    parser.add_argument("--sel", help='Selection to be applied to input tree')
    args = parser.parse_args()

    t0 = r.TChain("events")
    map(t0.Add, args.infiles)

    outfile = r.TFile(args.out, "recreate")

    print "Cleaning tracks..."
    t1 = clean_tracks(t0, min_pix_tr=args.min, iso_thresh=args.iso, fit=args.fit, sel=args.sel, cluster=args.cluster)
    outfile.Write()
    outfile.Close()
    print "Done! Wrote to %s." % args.out

    
