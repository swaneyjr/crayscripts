#!/usr/bin/env python

import ROOT as r
import numpy as np
from itertools import starmap, product, izip
from collections import namedtuple
import sys

Pixel = namedtuple('Pixel', ['x','y','val','avg3','avg5'])

def pix_mask(p):
    for x,y in BAD_PIX:
        if ((p.x-x)**2 + (p.y-y)**2) < (PIX_RADIUS**2): return False
    return True

def pix_stats(t0):
    rx,ry = args.res
    t0.Draw("pix_y:pix_x>>pstats({rx},0,{rx},{ry},0,{ry})".format(rx=rx, ry=ry), "", "goff")
    h = r.gROOT.FindObject("pstats")
    h.Scale(1./t0.GetEntries())

    arr = np.zeros(args.res, dtype=float)
    for x,y in product(xrange(rx), xrange(ry)):
        i = h.GetBin(x+1,y+1)
        arr[x,y] = h.GetBinContent(i)
    return arr

def draw_pix_stats(stats):
    import uuid
    hname = str(uuid.uuid1())
    xmin = np.floor(np.log10(np.min(stats[stats>0])))
    htemp = r.TH1F(hname, hname, 50, xmin, 0)
    map(htemp.Fill, np.log10(stats).flat)
    htemp.Draw()
    r.gROOT.FindObject("c1").Update()
    raw_input("Press enter")

def vbranch(t, bname, btype=float):
    btype_name = {float: 'double', int: 'int'}[btype]
    v = r.vector(btype_name)()
    setattr(t, '_%s'%bname, v)
    t.Branch(bname, v)

def clean_pix(t0, thresh, mask=None, stats=None, keep_empty=False, bad_regions=None, drop_source=False, l2thresh=None):
    if stats is None:
        print "Calculating pixel stats..."
        stats = pix_stats(t0)
        print "stats finished."

    def in_mask(x,y):
        if mask == None: return False
        for imask,(mx,my,mR) in enumerate(mask):
            if ((x-mx)**2 + (y-my)**2) < mR**2: return imask+1
        return 0
    
    def region_mask(p):
        if bad_regions == None: return False
        neval = 0
        for xlo,ylo,xhi,yhi in bad_regions:
            if (p.x>=xlo and p.x<=xhi and p.y>=ylo and p.y<=yhi): return True
        return False


    existing_branches = [b.GetName() for b in t0.GetListOfBranches()]
    print "calling clonetree()"
    t1 = t0.CloneTree(0)
    print "clonetree done"
    vbranch(t1, 'pix_freq')
    if 'pix_masked' in existing_branches:
        t1._pix_masked = r.vector('int')()
        t1.SetBranchAddress('pix_masked', t1._pix_masked)
    else:
        vbranch(t1, 'pix_masked', int)
    pix_n = np.array([0], dtype=int)
    t1.SetBranchAddress('pix_n', pix_n)

    total_pix = 0
    saved_pix = 0
    input_events = t0.GetEntries()
    print 'starting loop.'
    sys.stdout.flush()
    for ievt,evt in enumerate(t0):
        if ievt%(input_events/10) == 0:
            print "%d/%d  %.1f%%" % (ievt, input_events, 100.*ievt/input_events)
        t1._pix_freq.clear()
        t1._pix_masked.clear()

        # first, load the pixel data out of the TTree vectors
        pixels = list(starmap(Pixel, izip(evt.pix_x, evt.pix_y, evt.pix_val, evt.pix_avg3, evt.pix_avg5)))

        # now clear the vectors so we can re-write them
        t1.pix_x.clear()
        t1.pix_y.clear()
        t1.pix_val.clear()
        t1.pix_avg3.clear()
        t1.pix_avg5.clear()

        total_pix += len(pixels)

        if l2thresh:
            pixels = filter(lambda p: p.val > l2thresh, pixels)

        for p in pixels:
            freq = stats[p.x, p.y]
            masked = in_mask(p.x, p.y)
            if freq > thresh and not masked: continue
            if drop_source and masked: continue
            if region_mask(p): continue

            t1._pix_freq.push_back(stats[p.x, p.y])
            t1._pix_masked.push_back(masked)
            t1.pix_x.push_back(p.x)
            t1.pix_y.push_back(p.y)
            t1.pix_val.push_back(p.val)
            t1.pix_avg3.push_back(p.avg3)
            t1.pix_avg5.push_back(p.avg5)

        if t1.pix_x.size() == 0 and not keep_empty: continue
        pix_n[0] = t1.pix_x.size()
        saved_pix += t1.pix_x.size()

        t1.Fill()

    print "Done!"
    print
    print "Input events:    %d" % input_events
    print "Output events:   %d\t(%.2f%%)" % (t1.GetEntries(), 100.*t1.GetEntries()/input_events)
    print "Pixels removed:  %d\t(%.2f%%)" % ((total_pix-saved_pix), 100.*(total_pix-saved_pix)/total_pix)
    print "Pixels saved:    %d\t(%.2f%%)" % (saved_pix, 100.*saved_pix/total_pix)

    return t1

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Remove hotcells from ntuple.")
    parser.add_argument("--source-mask", metavar="X0,Y0,R", nargs="+", help="A comma-delimited list of (x0, y0, radius) to use as a source mask.")
    parser.add_argument("--drop-source", action='store_true', help="Don't save pixels in the source mask")
    parser.add_argument("--bad-region", metavar="X0,Y0,X1,Y1", nargs="+", help="Comma-delimited list of (xlo,ylo,xhi,yhi) regions to remove as bad pixels.")
    parser.add_argument("--thresh", default=1., type=float, help="The frequency threshold for keeping pixels")
    parser.add_argument("--L2", type=int, help="Apply an L2 threshold filter")
    parser.add_argument("--out", metavar="OUTFILE", default="cleaned.root", help="The output file")
    parser.add_argument("--sel", help="An event-level selection to apply to the input tree(s) before copying.")
    parser.add_argument("--in", dest="infiles", metavar="INFILE", required=True, nargs="+", help="The input file(s)")
    parser.add_argument("--res", metavar="RESX,RESY", default='1920,1080', required=True)
    args = parser.parse_args()

    if args.source_mask:
        args.source_mask = [map(float, sm.split(',')) for sm in args.source_mask]

    if args.bad_region:
        args.bad_region = [map(int, br.split(',')) for br in args.bad_region]

    args.res = map(int, args.res.split(','))
    
    t0 = r.TChain("events")
    map(t0.Add, args.infiles)

    if args.sel:
        t0_pre = t0
        t0 = t0_pre.CopyTree(args.sel)

    outfile = r.TFile(args.out, "recreate")

    print "Filtering pixels..."
    t1 = clean_pix(t0, args.thresh, args.source_mask, bad_regions=args.bad_region, drop_source=args.drop_source, l2thresh=args.L2)
    outfile.Write()
    outfile.Close()
    print "Done! Wrote to %s." % args.out
