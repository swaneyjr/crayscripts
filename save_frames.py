import cv2
import numpy as np
from bg_threshold import find_bg
from imtools import outlier_cutoff

def save_frames(vids, l1thresh=None):
 
  imlist = []
  autothresh = (l1thresh == None)
  
  for fname in vids:
    print fname
    
    if autothresh:
      l1thresh = np.amax(find_bg([fname]))
      print "L1 Threshold: %d" % l1thresh
  
    fbase = fname.split('.')[0]
    cap = cv2.VideoCapture(fname)
    iframe = 0
    ret, frame = cap.read()
    
    while ret and cap.isOpened():
    
      if np.amax(frame) >= l1thresh:
        imname = fbase + '_f' + str(iframe) + '.jpg'
        imlist.append(imname)
        print "Writing to %s" % imname
        cv2.imwrite(imname, frame, CV_IMWRITE_JPEG_QUALITY_100)
        
        
      iframe += 1
      ret, frame = cap.read()
    cap.release()
      
  return imlist
        
    

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Converts video to JPEG for frames above L1 threshold')
  parser.add_argument("--in", dest='infiles', nargs='+', help='Video(s) to convert')
  parser.add_argument("--l1thresh", type=int, help='L1 threshold for keeping frames')
  args = parser.parse_args()
  
  imlist = save_frames(args.infiles, args.l1thresh)
