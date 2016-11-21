import numpy as np
import rawpy
from PIL import Image
from scipy.signal import convolve2d
import gzip

class ImGrid(np.ndarray):
  
  def __new__(cls, file_name, bands=None):
    
    zip_types = ['gz']
    
    extensions = file_name.split('.')[1:]
  
    # open file
    if extensions[-1] in zip_types:
      f = gzip.open(file_name)
      extensions = extensions[:-1]
    else:
      f = open(file_name)
    
    if is_raw(file_name):
      raw = rawpy.imread(f)
      imarray = np.array([raw.raw_image.copy()]).astype(int)
      bands = ['RAW']
      raw.close()
    
    # PIL
    elif is_img(file_name):
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

def unzipped_ext(file):
  zip_types = ['gz','tar']
  file_ext = file.split('.')
  while file_ext[-1] in zip_types:
    file_ext = file_ext[:-1]
  return file_ext[-1]

def is_type(file, ext):
  return unzipped_ext(file) == ext

def is_raw(file):
  raw_types = ['dng']
  for t in raw_types:
    if is_type(file, t):
      return True
  return False
  
def is_img(file):
  im_types = ['jpg','png','gif']
  for t in im_types:
    if is_type(file, t):
      return True
  return False
  
def is_video(file):
  vid_types = ['mp4']
  for t in vid_types:
    if is_type(file, t):
      return True
  return False


def outlier_cutoff(imarray, thresh=1):
    n_bands = imarray.shape[0]
    cutoff_vals = np.zeros(n_bands)
    median_vals = np.mean(np.median(imarray, axis=1), axis=1)
    empty_vals = [np.argwhere(np.bincount(imarray[cval].flatten())<thresh) for cval in xrange(n_bands)]
    for cval,vals in enumerate(empty_vals):
        above_median = vals[vals>median_vals[cval]]
        if len(above_median)>0:
            cutoff_vals[cval] = min(above_median)
        else:
            cutoff_vals[cval] = np.amax(np.amax(imarray[cval], axis=0), axis=0) + 1

    return cutoff_vals
