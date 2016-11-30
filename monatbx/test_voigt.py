import math
from cctbx.array_family import flex
import matplotlib.pyplot as plt

def voigt(x,sig):
  if nu < 0:
    nu = 0
  elif nu > 1:
    nu = 1

  f1 = nu * math.sqrt(math.log(2)/math.pi) * flex.exp(-4*math.log(2)*((x/sig)**2)) * (1/abs(sig))
  f2 = (1-nu)/(math.pi*abs(sig)*(1+(4*((x/sig)**2))))
  f3 = ((nu * math.sqrt(math.log(2)/math.pi))/abs(sig)) + ((1-nu)/(math.pi*abs(sig)))
  svx = (f1 + f2)/f3
  return svx

sig = 0.003
nu = 0.5
x = flex.double(range(-1000,1000,1))
x = x/100000
fn = voigt(x, sig)
plt.plot(x, fn)
plt.show()
