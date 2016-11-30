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

  x_start, y_start = (0, 0)
  mask_w, mask_h = (4096, 4096)

  txt_out = ''
  for i in range(x_start, x_start + mask_h):
    for j in range(y_start, y_start + mask_w):
      txt_out += str(data[i][j]) + ' '
    txt_out += '\n'
  f = open(dstdir+'/'+item+'.txt', 'w')
  f.write(txt_out)
  f.close()
  print frame_no, fullpath, data_len

if (__name__ == "__main__"):

  data=sys.argv[1]
  destdir=sys.argv[2]
  n_proc = 16

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
