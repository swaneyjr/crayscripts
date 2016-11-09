import cv2
import numpy as np

def save_frames(infiles, l1thresh=0):
  imlist = []
  
  for fname in infiles:
  
    fbase = fname.split('.')[0]
    cap = cv2.VideoCapture(fname)
    iframe = 0
    ret, frame = cap.read()
    
    while ret and cap.isOpened():
    
      if frame.amax() >= l1thresh:
        imname = fbase + '_f' + str(iframe) + '.jpg'
        imlist.append(imname)
        cv2.imwrite(imname, frame, CV_IMWRITE_JPEG_QUALITY=100)
        
        
      iframe += 1
      ret, frame = cap.read()
      
  return imlist
        
    

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Converts video to JPEG for frames above L1 threshold')
  parser.add_argument("--in", dest=infiles, n_args='+', help='Video(s) to convert')
  parser.add_argument("--l1thresh", type=int, help='L1 threshold for keeping frames')
  parser.parse_args()
  
  imlist = save_frames(args.infiles, args.l1thresh)
