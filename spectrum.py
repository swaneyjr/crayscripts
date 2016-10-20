import imtools
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('tkagg')

def plot_spectrum(images, region, average):

  adc_counts = 0
  
  for i, im in enumerate(images):
  
    print
    
   
    imarray = imtools.ImGrid(im)
    adc_counts += imtools.find_spectrum(imarray)
    
  plt.figure(1)
  plt.hist(adc_counts, )
    

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser(description='')
  parser.add_argument("--in", dest='images', required=True, nargs='+', help='')
  parser.add_argument("--region,
  parser.add_argument("--average", help='')
  args = parser.parse_args()
  
