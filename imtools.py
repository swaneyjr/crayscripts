import numpy as np
import rawpy
from PIL import Image
from scipy.signal import convolve2d
import gzip

class Imarray(np.array):

# adds a vector TBranch to a TTree
def vbranch(t, bname, btype=float):
    btype_name = {float: 'double', int: 'int'}[btype]
    v = r.vector(btype_name)()
    setattr(t, '_%s'%bname, v)
    t.Branch(bname, v)

# returns an array of counts of ADC values in a region
def find_spectrum(imarray, counts=1024, region=(0,0)+imarray.shape):
  return np.bincount(imarray[region[0]:region[2],region[1]:region[3]], minlength=counts)

# opens an image file and returns an array of dimensions (h,w,band)
def get_imarray(file_name, israw=False):
  
  if not raw:
    return np.array(Image.open(file_name))
  elif file_name.endswith('.gz'):
    f = gzip.open(file_name)
    file_name = file_name[:-3]
  else:
    f = open(file_name)
  raw = rawpy.imread(f)
  imarray = raw.raw_image().copy()
  raw.close()
  return np.array([imarray])
  
