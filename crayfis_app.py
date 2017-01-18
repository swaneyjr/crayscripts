import os
from save_frames import save_frames
from img_to_root import convert_to_root

def crayfis_app(vid, l1rate=0.1, spatial_thresh=False):
    vid_name = vid.split('/')[-1].split('.')[0]
    os.mkdir(vid_name)
    os.chdir(vid_name)
    save_frames([os.getcwd()+'/../'+vid], l1thresh=0)
    convert_to_root(os.listdir(os.getcwd()), vid_name+'.root', l1rate, s_thresh=spatial_thresh)
    
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument("vid", help='Video to be converted')
    parser.add_argument('-S',"--s_thresh", action='store_true', help='Use spatial thresholds')
    parser.add_argument("--l1_rate", default=0.1, help='Target rate for passing L1')
    args = parser.parse_args()

    crayfis_app(args.vid, args.l1_rate, args.s_thresh)
