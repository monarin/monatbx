import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
x = np.array([35,38,53,57,45,47,52,54,15,48,35,38,53,55,43,48,56,52,45,47,33,39,53,57,45,47,54,50,35,44,35,38,53,57,45,50,51,54,40,49])
mu = np.mean(x)
med = np.median(x)
sigma = np.std(x)
num_bins = 10
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.ylabel('Frequencies')
plt.title('Age distribution')
plt.show()
