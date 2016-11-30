import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys
import numpy as np
from cctbx.array_family import flex

fname = sys.argv[1]

#read .paramhist files and display refinement results
pf = open(fname,'r')
data = pf.read().split('\n')
vol_set = flex.double()
n_col = 20
for data_row in data:
  dc = data_row.split()
  if len(dc)==n_col:
    vol_set.append(float(dc[12])*float(dc[13])*float(dc[14]))

print 'N_images=%8.0f'%(len(vol_set))
x = vol_set.as_numpy_array()
mu = np.mean(x)
med = np.median(x)
sigma = np.std(x)
num_bins = 50
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.ylabel('Frequencies')
plt.title('beta distribution\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))
plt.show()
