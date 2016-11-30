import cPickle as pickle
import os,sys
import fnmatch
from libtbx.easy_mp import pool_map

def zero_mproc(frame_no, frame_files, dstdir):
  fullpath = frame_files[frame_no]
  G = pickle.load(open(fullpath,"rb"))
  fullpath_arr = fullpath.split('/')
  item = fullpath_arr[len(fullpath_arr)-1]

  data = G["DATA"]
  len_data_0 = (data==0).count(True)
  len_data_20 = (data<=0).count(True)

  zdata = data.set_selected(data<=0, 0)
  len_data_0_after = (zdata==0).count(True)
  G["DATA"] = zdata
  pickle.dump(G,open(os.path.join(dstdir,item),"wb"),pickle.HIGHEST_PROTOCOL)
  print frame_no, item, len_data_0, len_data_20, len_data_0_after

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
