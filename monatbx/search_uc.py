import cPickle as pickle
import os,sys
import fnmatch
from libtbx.easy_mp import pool_map
from iotbx import reflection_file_reader
from subprocess import call
import numpy as np

def zero_mproc(frame_no, frame_files, a_range, b_range, c_range, data, model):
  fname = frame_files[frame_no]
  fname_arr = fname.split('_')
  a = a_range[int(fname_arr[0])]
  b = b_range[int(fname_arr[1])]
  c = c_range[int(fname_arr[2])]
  uc = str(a)+','+str(b)+','+str(c)+',90,90,90'
  cb = call(["phenix.refine", "inp_0416_008.eff","refinement.input.pdb.file_name="+model, \
  "refinement.input.xray_data.file_name="+data,"refinement.input.xray_data.r_free_flags.file_name="+data, \
  "refinement.output.prefix=SFC_"+fname,"refinement.crystal_symmetry.unit_cell="+uc])
  print frame_no, fname, a, b, c
  return fname


if (__name__ == "__main__"):
  #read input
  data=sys.argv[1]
  model=sys.argv[2]
  n_proc = 32
  #grab unitcell from reflection file
  reflection_file = reflection_file_reader.any_reflection_file(data)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array = miller_arrays[0]
  miller_array.show_summary()
  #a,b,c = miller_array.unit_cell().parameters()[:3]
  a,b,c = (69.65700378,170.5290039,290.8610046)
  a_range = np.arange(a-0.3,a+0.3,0.15)
  b_range = np.arange(b-0.6,b+0.6,0.15)
  c_range = np.arange(c-0.6,c+0.6,0.15)
  print 'Unit-cell search range'
  print a_range
  print b_range
  print c_range
  frame_files = []
  for i in range(len(a_range)):
    for j in range(len(b_range)):
      for k in range(len(c_range)):
        frame_files.append(str(i)+'_'+str(j)+'_'+str(k))
  def zero_mproc_wrapper(arg):
    return zero_mproc(arg, frame_files, a_range, b_range, c_range, data, model)
  frames = range(len(frame_files))
  result = pool_map(
            iterable=frames,
            func=zero_mproc_wrapper,
            processes=n_proc)
