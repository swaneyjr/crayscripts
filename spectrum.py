import imtools
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('tkagg')
plt.ion()

def plot_spectrum(images, rg, average):

  adc_counts = 0
  
  for i, im in enumerate(images):
  
    print
    
   
    imarray = imtools.ImGrid(im)
    adc_counts += imtools.find_spectrum(imarray, region=rg)
    
  plt.figure(1)
  plt.hist(np.arange(counts), weights=adc_counts)
  plt.xlabel('ADC Counts')
  plt.ylabel('Frequency')
  plt.show()
    

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='')
  parser.add_argument("--in", dest='images', required=True, nargs='+', help='')
  parser.add_argument("--region, help='')
  parser.add_argument("--average", help='')
  args = parser.parse_args()
                      
  
