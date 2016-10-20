import numpy as np
import rawpy
from PIL import Image
from scipy.signal import convolve2d
import gzip

class ImGrid(np.array):
  def __new__(cls, input_array, bands=None):
    obj = np.asarray(input_array).view(cls)
    obj.bands = ###
    return obj

  def __array_finalize__(self, obj):
    if obj is None: return
    self.info = getattr(obj, 'info', None)

# returns an array of counts of ADC values in a region
def spectrum(imarray, counts=1024, region=(0,0)+imarray.shape):
  return np.array([np.bincount(imarray[region[0]:region[2],region[1]:region[3], cval], minlength=counts) \
                   for cval in xrange(imarray.shape[2])])

# opens an image file and returns an array of dimensions (h,w,band)
def get_imarray(file_name):
  
  raw_types = ['dng']
  
  # PIL
  if not file_name.split('.')[1] in raw_types:
    return np.array(Image.open(file_name))
  
  # rawpy
  elif file_name.endswith('.gz'):
    f = gzip.open(file_name)
    file_name = file_name[:-3]
  else:
    f = open(file_name)
  raw = rawpy.imread(f)
  imarray = raw.raw_image().copy()
  raw.close()
  return np.array([imarray])
  
