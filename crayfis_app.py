#!/usr/bin/env python2

import os
import ROOT as r
from save_frames import save_frames
from img_to_root import convert_to_root

def crayfis_app(vid, l1rate=0.1, spatial_thresh=False):
    vid_name = vid.split('/')[-1].split('.')[0]
    os.mkdir(vid_name)
    os.chdir(vid_name)
    save_frames([os.getcwd()+'/../'+vid], l1thresh=0)
    imfiles = os.listdir(os.getcwd())
    f = r.TFile(vid_name + '.root','recreate')
    t=convert_to_root(imfiles, l1rate, s_thresh=spatial_thresh)
    f.Write()
    f.Close()
    os.rename(os.getcwd()+'/'+vid_name+'.root',os.pardir+'/'+vid_name+'.root')

    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument("vid", help='Video to be converted')
    parser.add_argument('-S',"--s_thresh", action='store_true', help='Use spatial thresholds')
    parser.add_argument("--l1_rate", type=float, default=0.1, help='Target rate for passing L1')
    args = parser.parse_args()

    crayfis_app(args.vid, args.l1_rate, args.s_thresh)
