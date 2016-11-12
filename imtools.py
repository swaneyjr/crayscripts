import numpy as np
import rawpy
from PIL import Image
from scipy.signal import convolve2d
import gzip

class ImGrid(np.ndarray):
  
  def __new__(cls, file_name, bands=None):
    
    compressed_types = ['jpg','png','gif']
    zip_types = ['gz']
    
    extensions = file_name.split('.')[1:]
  
    # open file
    if extensions[-1] in zip_types:
      f = gzip.open(file_name)
      extensions = extensions[:-1]
    else:
      f = open(file_name)
    
    if extensions[-1] not in compressed_types:
      raw = rawpy.imread(f)
      imarray = np.array([raw.raw_image.copy()]).astype(int)
      bands = ['RAW']
      raw.close()
    
    # PIL
    else:
      with Image.open(f) as im:
        imarray = np.array(im).astype(int).transpose(2,0,1)
        bands = list(im.mode)
        if len(bands) == 1:
          imarray = np.array([imarray])
          
    f.close()
    
    obj = np.asarray(imarray).view(cls)
    obj.bands = bands
    obj.n_bands = len(bands)
    obj.height = imarray.shape[1]
    obj.width = imarray.shape[2]
    return obj


  def __array_finalize__(self, obj):
    if obj is None: return
    self.bands = getattr(obj, 'bands', None)
    self.height = getattr(obj, 'height', None)
    self.width = getattr(obj, 'width', None)
    self.n_bands = getattr(obj, 'n_bands', None)


# returns an array of counts of ADC values in a region
def spectrum(imarray, counts=1024, region=None):
  if not region:
    region = (0,0,imarray.width, imarray.height)
  return np.array([np.bincount(imarray[region[1]:region[3],region[0]:region[2], cval].flatten(), minlength=counts) \
                   for cval in xrange(imarray.shape[2])])
  
def is_type(file, ext):
  zip_types = ['gz','tar']
  file_ext = file.split('.')[1:]
  if file_ext[-1] in zip_types:
    file_ext = file_ext[:-1]
  return file[-1] == ext

def is_raw(file):
  raw_types = ['dng']
  for t in raw_types:
    if is_type(file, t):
      return True
  return False
  
def is_img(file):
  im_types = ['jpg']
  for t in im_types:
    if is_type(file, t):
      return True
  return False
  
def is_video(file):
  vid_types = ['mp4']
  for t in raw_types:
    if is_type(file, t):
      return True
  return False
