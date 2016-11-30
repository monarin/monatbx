import dxtbx, sys, os, math
from skimage import io, color, measure, draw, img_as_bool
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

def cost(params):
  x0, y0, r = params
  coords = draw.circle(y0, x0, r, shape=image.shape)
  template = np.zeros_like(image)
  template[coords] = mean_i
  print x0, y0, r
  return np.sum((template[coords] - image[coords])**2)

if (__name__ == "__main__"):
  imgpath = sys.argv[1]

  img = dxtbx.load('mccd/E1_0_00006_106.mccd')
  raw_data = img.get_raw_data()
  image = raw_data.as_numpy_array()
  detector = img.get_detector()
  detector = detector[0]
  beam = img.get_beam()
  pixel_size = detector.get_pixel_size()[0]
  x0 = int(round(detector.get_beam_centre(beam.get_s0())[0] / pixel_size))
  y0 = int(round(detector.get_beam_centre(beam.get_s0())[1] / pixel_size))
  r = 1400

  coords = draw.circle(y0, x0, r, shape=image.shape)
  mean_i = np.mean(image[coords])
  print x0, y0, r, mean_i

  """
  img = io.imread(imgpath)
  img_gray = color.rgb2gray(img)
  image = img_as_bool(img_gray)
  regions = measure.regionprops(image)
  bubble = regions[0]

  y0, x0 = bubble.centroid
  r = bubble.major_axis_length / 2.


  """
  x0, y0, r = optimize.fmin(cost, (x0, y0, r))
  import matplotlib.pyplot as plt

  f, ax = plt.subplots()
  circle = plt.Circle((x0, y0), r)
  plt.title('x0=%6.2f y0=%6.2f r=%6.2f'%(x0,y0,r))
  ax.imshow(image, cmap='gray', interpolation='nearest')
  ax.add_artist(circle)
  plt.show()
