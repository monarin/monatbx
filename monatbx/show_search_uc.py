import cPickle as pickle
import os,sys
import fnmatch
from libtbx.easy_mp import pool_map
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.show()

def zero_mproc(frame_no, frame_files):
  fname = os.path.basename(frame_files[frame_no])
  fname_arr = fname.split('_')
  i,j,k,a,b,c,alpha,beta,gamma,dummy,rwork,rfree,bonds,angles,b_min,b_max,b_ave,n_water,shifts = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  if len(fname_arr) == 5:
    i,j,k = fname_arr[1:4]
    pf = open(frame_files[frame_no],'r')
    data = pf.read().split('\n')
    for row in data:
      if row.startswith('  refinement.crystal_symmetry.unit_cell = '):
        col1 = row.split('  refinement.crystal_symmetry.unit_cell = ')
        col2 = col1[1].split(',')
        a,b,c,alpha,beta,gamma = col2
      if row.startswith('        end:'):
        col = row.split()
        dummy,rwork,rfree,bonds,angles,b_min,b_max,b_ave,n_water,shifts = col
    print i,j,k, a,b,c,alpha,beta,gamma,rwork,rfree,bonds,angles,b_min,b_max,b_ave,n_water,shifts
  if a==0:
    return None
  else:
    return i,j,k, a,b,c,alpha,beta,gamma,rwork,rfree,bonds,angles,b_min,b_max,b_ave,n_water,shifts

if (__name__ == "__main__"):

  data=sys.argv[1]
  n_proc = 1
  frame_files = []
  for root, dirnames, filenames in os.walk(data):
    for filename in fnmatch.filter(filenames, '*.log'):
      fullpath = os.path.join(root, filename)
      frame_files.append(fullpath)

  def zero_mproc_wrapper(arg):
    return zero_mproc(arg, frame_files)

  frames = range(len(frame_files))
  result = pool_map(
            iterable=frames,
            func=zero_mproc_wrapper,
            processes=n_proc)
  good_result = []
  for res in result:
    if res is not None:
      good_result.append(res)
  a_set = np.array([float(res[3]) for res in good_result])
  b_set = np.array([float(res[4]) for res in good_result])
  c_set = np.array([float(res[5]) for res in good_result])
  rfree_set = np.array([float(res[10]) for res in good_result])
  i_min = np.argmin(rfree_set)
  print 'Minimum Rfree =%6.2f a=%8.4f b=%8.4f c=%8.4f'%(rfree_set[i_min]*100,a_set[i_min],b_set[i_min],c_set[i_min])
  scatter3d(a_set,b_set,c_set,rfree_set)
