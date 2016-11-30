import numpy as np
from cctbx.array_family import flex

G = [5, 6, 0.1, 20, 16, 12, 11, 10, 11.5, 15]
B = [100, 80, 200, 60, 70, 80, 85, 90, 70, 40]
rotx = [0.01, 0.002, 0.001, 0.05, 0.1, 0.025, 0.008, 0.01, 0.002, 0.001]
X = np.array([G, B, rotx])
print X
COV = np.cov(X)
print COV
CORR = np.correlate(X)
print CORR
