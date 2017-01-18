import os
from save_frames import save_frames
from img_to_root import convert_to_root

def crayfis_app(vid):
    vid_name = vid.split('/')[-1].split('.')[0]
    os.mkdir(vid_name)
    os.chdir(vid_name)
    save_frames(vid, l1thresh=0)
    os.chdir(os.pardir)
    convert_to_root(listdir(vid_name), vid_name+'.root')
    
    
    
    
    

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument("vid", help='Video to be converted')
    args = parser.parse_args()

    crayfis_app(args.vid)
