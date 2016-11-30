import dxtbx, sys, os, math
import numpy as np
from sklearn import mixture
import matplotlib.pyplot as plt

if (__name__ == "__main__"):
  imgpath = sys.argv[1]

  img = dxtbx.load('mccd/E1_0_00006_106.mccd')
  raw_data = img.get_raw_data()
  image = raw_data.as_numpy_array()

  data_1 = image[100][:]
  obs = np.array([0]*int(round(data_1[0])))
  for i in range(1,len(data_1)):
    x = np.array([i]*int(round(data_1[i])))
    obs = np.concatenate([obs, x])

  obs = np.reshape(obs, (len(obs),1))
  plt.subplot(311)
  plt.plot(data_1)
  plt.subplot(312)
  plt.plot(obs)
  plt.subplot(313)
  num_bins = 100
  n, bins, patches = plt.hist(obs, num_bins, normed=0, facecolor='green', alpha=0.5)
  plt.show()

  g = mixture.GMM(n_components=2)
  g.fit(obs)
  print np.round(g.means_, 2)
  print np.round(g.covars_, 2)
