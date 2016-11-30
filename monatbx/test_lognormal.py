import math
import numpy as np
from cctbx.array_family import flex
import matplotlib.pyplot as plt

def lognpdf(x, fwhm, zero):
  #find sig from root of this function
  sig_set = np.arange(50)/100
  t = sig_set * math.sqrt(math.log(4))
  fx = fwhm - (zero * (np.exp(t) - np.exp(-1*t)))
  sig = sig_set[np.argmin(np.abs(fx))]
  plt.plot(sig_set, fx)
  plt.title('sig_f0=%6.2f'%(sig))
  plt.show()
  #calc x0
  x0 = math.log(zero) + sig**2
  g = 1/( sig * math.sqrt(2*math.pi) * np.exp(x0-((sig**2)/2)) )
  #calc lognpdf
  X = zero - x
  f1 = 1/( X * sig * math.sqrt(2*math.pi) )
  f2 = np.exp( -1 * (np.log(X)-x0)**2 / (2*(sig**2)) )
  svx = flex.double(f1 * f2 / g)
  return svx

zero = 0.008
fwhm = 0.005
x = flex.double(range(-1000,1000,1))
x = x/100000
fn = lognpdf(x, fwhm, zero)
plt.plot(x, fn)
plt.show()
