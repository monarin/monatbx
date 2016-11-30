import fnmatch
import os, sys
import cPickle as pickle
from libtbx.easy_mp import pool_map

def img2pickle_mproc(frame_no, frame_files, file_done_dict, lambda_pickle):
  fullpath = frame_files[frame_no]
  #cmd='cxi.image2pickle -c -x 1886.65 -y 2099.13 '+fullpath
  #cmd='cxi.image2pickle -d 153 '+fullpath
  filename = fullpath.split('/')[-1].split('.')[0]
  lambda_new = 0
  if filename in lambda_pickle:
    lambda_new = lambda_pickle[filename]
  if lambda_new > 0:
    cmd='cxi.image2pickle -w '+str(lambda_new)+' '+fullpath
  else:
    cmd='cxi.image2pickle '+fullpath
  os.system(cmd)
  print filename, cmd
  flag_done = True
  return fullpath, flag_done

if (__name__ == "__main__"):
  data=sys.argv[1]
  n_proc = 32
  #read lambda_pickle
  lambda_pickle = {}
  if len(sys.argv) > 2:
    lambda_pickle = pickle.load(open(sys.argv[2], 'rb'))
  #in case done.lst is there - exclude all files in this list
  file_done_dict = {}
  if os.path.isfile('done.lst'):
    file_done = open('done.lst','r')
    file_done_list=file_done.read().split("\n")
    file_done_dict.update(dict([(f, 0) for f in file_done_list]))
  #start reading data
  frame_files = []
  for root, dirnames, filenames in os.walk(data):
    for filename in fnmatch.filter(filenames, '*.mccd'):
      fullpath = os.path.join(root, filename)
      frame_files.append(fullpath)
  #process mproc
  def img2pickle_mproc_wrapper(arg):
    return img2pickle_mproc(arg, frame_files, file_done_dict, lambda_pickle)
  frames = range(len(frame_files))
  result = pool_map(
            iterable=frames,
            func=img2pickle_mproc_wrapper,
            processes=n_proc)
  file_all_dict = {}
  cn_done = 0
  cn_new = 0
  for res in result:
    if res is not None:
      filenew_name, flag_done = res
      if flag_done:
        cn_done += 1
      else:
        cn_new += 1
        file_all_dict[filenew_name] = 0

  print len(frame_files), cn_done, cn_new, len(file_all_dict.keys())
