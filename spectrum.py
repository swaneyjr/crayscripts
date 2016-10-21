import imtools
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('tkagg')
plt.ion()

# make plots for each color band of distribution of ADC 
def plot_spectrum(images, rg, normalize):

  adc_counts = 0
  
  for i, im in enumerate(images):
  
   
    imarray = imtools.ImGrid(im)
    adc_counts += imtools.find_spectrum(imarray, region=rg)
    
  if normalize:
    adc_counts /= float(adc_counts.sum())
    
  max_count = max(adc_counts[0].nonzero())
  if max_count < 256:
    max_count = 256
    adc_counts = adc_counts[:,256]
  else:
    max_count = 1024
    
  for cval in xrange(adc_counts.shape[0]):
    plt.figure(cval)
    plt.hist(np.arange(max_count), weights=adc_counts)
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
  
  plot_spectrum(args.images, region, args.normalize)
  
  
