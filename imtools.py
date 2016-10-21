import numpy as np
import rawpy
from PIL import Image
from scipy.signal import convolve2d
import gzip

class ImGrid(np.ndarray):
  def __new__(cls, file_name, bands=None):
    
    raw_types = ['dng']
  
    # PIL
    if not file_name.split('.')[1] in raw_types:
      with Image.open(file_name) as im:
        imarray = np.array(im).astype(int)
        bands = list(im.mode)
        if len(bands) == 1:
          imarray = np.array([imarray]).transpose(1,2,0)
  
    # rawpy
    else:
      if file_name.endswith('.gz'):
        f = gzip.open(file_name)
      else:
        f = open(file_name)
      raw = rawpy.imread(f)
      imarray = np.array([raw.raw_image().copy()]).transose(1,2,0)
      bands = ['RAW']
      raw.close()
      f.close()
  
    obj = np.asarray(imarray).view(cls)
    obj.bands = bands
    return obj

  def __array_finalize__(self, obj):
    if obj is None: return
    self.bands = getattr(obj, 'bands', None)

# returns an array of counts of ADC values in a region
def spectrum(imarray, counts=1024, region=(0,0)+imarray.shape):
  return np.array([np.bincount(imarray[region[0]:region[2],region[1]:region[3], cval], minlength=counts) \
                   for cval in xrange(imarray.shape[2])])
  
