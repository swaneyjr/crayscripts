import cv2
import numpy as np
from bg_threshold import find_bg
from imtools import outlier_cutoff
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('tkagg')

def save_frames(vids, l1thresh=None):
 
  imlist = []
  
  if l1thresh == None:
   
    print "Calculating video stats..."
    
    adc_array = np.zeros(256, dtype=int)
    iframe = 0
    
    print "Frames processed:"
    print " 0"
    
    for fname in vids:
      cap = cv2.VideoCapture(fname)
      ret, frame = cap.read()
      iframe += 1
      while ret and cap.isOpened():
        adc_array[np.amax(frame)] += 1
        ret, frame = cap.read()
        iframe += 1
        if iframe % 1000 == 0:
          print " %d" % iframe 
      cap.release()
        
    figure = plt.figure()
    ax = figure.add_subplot(111)
    plt.hist(np.arange(256), weights=adc_array, bins=256)
    ax.set_xlabel('ADC count')
    ax.set_ylabel('Frequency')
    ax.set_title('Max ADC count by frame')
    plt.show()
    
    l1thresh = int(raw_input('L1 Threshold: '))
  
  for fname in vids:
    print fname
      
    fbase = fname.split('.')[0]
    cap = cv2.VideoCapture(fname)
    iframe = 0
    ret, frame = cap.read()
    
    while ret and cap.isOpened():
    
      if np.amax(frame) >= l1thresh:
        imname = fbase + '_f' + str(iframe) + '.jpg'
        print "Writing to %s" % imname
        cv2.imwrite(imname, frame, [int(cv2.IMWRITE_JPEG_QUALITY), 100])
        
        
      iframe += 1
      ret, frame = cap.read()
    cap.release()
        
    

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Converts video to JPEG for frames above L1 threshold')
  parser.add_argument("--in", dest='infiles', nargs='+', help='Video(s) to convert')
  parser.add_argument("--l1thresh", type=int, help='L1 threshold for keeping frames')
  args = parser.parse_args()
  
  save_frames(args.infiles, args.l1thresh)
