import fnmatch
import os, sys
from cctbx.array_family import flex
import cPickle as pickle
from libtbx.easy_mp import pool_map
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def calc_I_mproc(frame_no, frame_files):
  fullpath = frame_files[frame_no]
  #G = pickle.load(open(fullpath,"rb"))["DATA"].as_numpy_array()[:,range(550,950)][range(350,750), :]
  #img = G["DATA"].as_numpy_array()[:,range(550,950)][range(350,750), :]

  #define area to calculate intensity stats
  #col_img = data[:,range(550,950)]
  #img = col_img[range(350,750), :]
  e = float(np.max(pickle.load(open(fullpath,"rb"))["DATA"].as_numpy_array()[:,range(1400,1800)][range(1400,1800), :]))
  print '%6.0f %10.1f'%(frame_no, e), fullpath
  return 0, 0, 0, 0, e, fullpath

if (__name__ == "__main__"):

  data=sys.argv[1]
  n_proc = 16

  frame_files = []
  for root, dirnames, filenames in os.walk(data):
    for filename in fnmatch.filter(filenames, '*.pickle'):
      fullpath = os.path.join(root, filename)
      frame_files.append(fullpath)
  def calc_I_mproc_wrapper(arg):
    return calc_I_mproc(arg, frame_files)
  frames = range(len(frame_files))
  result = pool_map(
            iterable=frames,
            func=calc_I_mproc_wrapper,
            processes=n_proc)
  txt_out = ''
  mx_set = []
  mean_set = []
  for res in result:
    if res is not None:
      mean, med, std, mn, mx, fullpath = res
      mx_set.append(mx)
      mean_set.append(mean)
      if mx > 60:
        txt_out += fullpath + '\n'
  print np.mean(mean_set), np.std(mean_set)
  print np.mean(mx_set), np.std(mx_set)
  f = open('non_blank.lst', 'w')
  f.write(txt_out)
  f.close()
  plt.hist(mx_set, 20, normed=0)
  plt.show()
  """
  plt.subplot(321)
  x = mean_I_set.as_numpy_array()
  mu = np.mean(x)
  med = np.median(x)
  sigma = np.std(x)
  num_bins = 20
  n, bins, patches = plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
  y = mlab.normpdf(bins, mu, sigma)
  plt.plot(bins, y, 'r--')
  plt.ylabel('Frequencies')
  plt.title('Mean I distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

  plt.subplot(322)
  x = med_I_set.as_numpy_array()
  mu = np.mean(x)
  med = np.median(x)
  sigma = np.std(x)
  num_bins = 20
  n, bins, patches = plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
  y = mlab.normpdf(bins, mu, sigma)
  plt.plot(bins, y, 'r--')
  plt.ylabel('Frequencies')
  plt.title('Median I distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

  plt.subplot(323)
  x = std_I_set.as_numpy_array()
  mu = np.mean(x)
  med = np.median(x)
  sigma = np.std(x)
  num_bins = 20
  n, bins, patches = plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
  y = mlab.normpdf(bins, mu, sigma)
  plt.plot(bins, y, 'r--')
  plt.ylabel('Frequencies')
  plt.title('Std. I distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

  plt.subplot(324)
  x = min_I_set.as_numpy_array()
  mu = np.mean(x)
  med = np.median(x)
  sigma = np.std(x)
  num_bins = 20
  n, bins, patches = plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
  y = mlab.normpdf(bins, mu, sigma)
  plt.plot(bins, y, 'r--')
  plt.ylabel('Frequencies')
  plt.title('Min I distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

  plt.subplot(325)
  x = max_I_set.as_numpy_array()
  mu = np.mean(x)
  med = np.median(x)
  sigma = np.std(x)
  num_bins = 20
  n, bins, patches = plt.hist(x, num_bins, normed=0, facecolor='green', alpha=0.5)
  y = mlab.normpdf(bins, mu, sigma)
  plt.plot(bins, y, 'r--')
  plt.ylabel('Frequencies')
  plt.title('Max I distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
  plt.show()
  """
