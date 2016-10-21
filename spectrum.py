import imtools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from time import clock

matplotlib.use('tkagg')

# make plots for each color band of distribution of ADC 
def plot_spectrum(images, rg, normalize):

  adc_counts = 0
  n_images = len(images)
  print "Calculating ADC Counts..."
  for i, im in enumerate(images):
  
    if i%(n_images/10) == 0:
      print "%d/%d (%.1f%%)" % (i, n_images, 100.*i/n_images)
    
    imarray = imtools.ImGrid(im)
    adc_counts += imtools.spectrum(imarray, region=rg)
    
  if normalize:
    adc_counts /= float(adc_counts.sum())
    
  if np.count_nonzero(adc_counts[:,256:]) == 0:
    adc_counts = adc_counts[:,:256]
    
  max_count = adc_counts.shape[1]
    
  for cval in xrange(adc_counts.shape[0]):
    plt.figure(cval+1)
    plt.hist(np.arange(max_count), bins = max_count, weights=adc_counts[cval], log=True)
    plt.xlabel('ADC Counts')
    plt.ylabel('Frequency')
  plt.show()
    

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Generates histograms of ADC counts of images by color band')
  parser.add_argument("--in", dest='images', required=True, nargs='+', help='Images from which spectrum is to be constructed')
  parser.add_argument("--region", metavar="X0,Y0,X1,Y1", help='A comma delimited list of the region to compute spectrum across images')
  parser.add_argument("--normalize", action='store_true', help='Normalize data')
  args = parser.parse_args()
  
  if args.region:
    region = tuple(map(int, args.region.split(',')))
  else:
    region = None
  
  ti = clock()
  plot_spectrum(args.images, region, args.normalize)
  tf = clock()
  
  m,s = divmod(tf-ti,60)
  h,m = divmod(m,60)
    
  print "Done!"
  print "Total time: ",
  if tf-ti > 3600:
    print "%d h %02d m %02d s" % (h,m,s)
  elif tf-ti > 60:
    print "%d m %02d s" % (m,s)
  else:
    print "%f s" % (tf-ti)
