import numpy as np
import rawpy
from PIL import Image
from scipy.signal import convolve2d
import gzip

class ImGrid(np.ndarray):
  def __new__(cls, file_name, bands=None):
    
    raw_types = ['dng']
    img_types = ['jpg']
    zip_types = ['gz']
    
    extensions = file_name.split('.')
  
    # open file
    if extensions[-1] in zip_types:
      f = gzip.open(file_name)
      extensions = extensions[:-1]
    else:
      f = open(file_name)
    
    if extensions[-1] in raw_types:
      raw = rawpy.imread(f)
      imarray = np.array([raw.raw_image.copy()]).transpose(1,2,0).astype(int)
      bands = ['RAW']
      raw.close()
    
    # PIL
    else: 
      imarray = np.array(im).astype(int)
      bands = list(im.mode)
      if len(bands) == 1:
        imarray = np.array([imarray]).transpose(1,2,0)
          
    f.close()
    
    obj = np.asarray(imarray).view(cls)
    obj.bands = bands
    obj.height = imarray.shape[0]
    obj.width = imarray.shape[1]
    return obj

  def __array_finalize__(self, obj):
    if obj is None: return
    self.bands = getattr(obj, 'bands', None)
    self.height = getattr(obj, 'height', None)
    self.width = getattr(obj, 'width', None)

# returns an array of counts of ADC values in a region
def spectrum(imarray, counts=1024, region=None):
  if not region:
    region = (0,0,imarray.width, imarray.height)
  return np.array([np.bincount(imarray[region[1]:region[3],region[0]:region[2], cval].flatten(), minlength=counts) \
                   for cval in xrange(imarray.shape[2])])
  
