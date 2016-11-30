import cPickle as pickle
from cctbx.array_family import flex
import os,sys
import fnmatch
from libtbx.easy_mp import pool_map
import numpy as np

def zero_mproc(frame_no, frame_files, dstdir):
  fullpath = frame_files[frame_no]
  G = pickle.load(open(fullpath,"rb"))
  fullpath_arr = fullpath.split('/')
  item = fullpath_arr[len(fullpath_arr)-1]
  data = G["DATA"].as_numpy_array()

  data_len = data.size
  len_data_0 = data_len - np.count_nonzero(data)
  #define mask

  #L650-hydrognase mask
  #x_start, y_start = (1829, 1786)
  #mask_w, mask_h = (183, 197)

  #SFC mask
  #x_start, y_start = (1519,1531)
  #mask_w, mask_h = (504, 464)

  #ACT short mask
  #x_start, y_start = (898, 567)
  #mask_w, mask_h = (1351, 202)
  #ACT beam center mask
  #x_start, y_start = (898, 898)
  #mask_w, mask_h = (135, 135)
  #Mathews HEWL mask
  #x_start, y_start = (1125, 1140)
  #mask_w, mask_h = (120, 120)
  #SFC II-20B
  #x_start, y_start = (1918,0)
  #mask_w, mask_h = (3600,272)
  #munc 81B-58
  #x_start, y_start = (1900,0)
  #mask_w, mask_h = (4096,286)
  #Lior 2016
  x_start, y_start = (2030,0)
  mask_w, mask_h = (4096,60)

  #start masking
  for i in range(x_start, x_start + mask_h):
    data[i][y_start:y_start+mask_w] = -2

  len_data_0_after = data_len - np.count_nonzero(data)

  G["DATA"] = flex.int(data)
  pickle.dump(G,open(os.path.join(dstdir,item),"wb"),pickle.HIGHEST_PROTOCOL)
  print frame_no, item, len_data_0, len_data_0_after

if (__name__ == "__main__"):

  data=sys.argv[1]
  destdir=sys.argv[2]
  n_proc = 1

  frame_files = []
  for root, dirnames, filenames in os.walk(data):
    for filename in fnmatch.filter(filenames, '*.pickle'):
      fullpath = os.path.join(root, filename)
      frame_files.append(fullpath)

  def zero_mproc_wrapper(arg):
    return zero_mproc(arg, frame_files, destdir)

  frames = range(len(frame_files))
  result = pool_map(
            iterable=frames,
            func=zero_mproc_wrapper,
            processes=n_proc)
