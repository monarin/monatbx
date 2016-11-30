'''
Author      : Uervirojnangkoorn, M.
Created     : 7/29/2015
Description : Partiality calculation

'''
from __future__ import division
from scipy import optimize
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation, basis_type
from cctbx import statistics

class partiality_handler(object):
  '''
  A wrapper class for partiality calculation
  '''

  def __init__(self, iparams):
    '''
    Intialitze parameters
    '''
    self.iparams = iparams

  def calc_partiality_anisotropy_set(self, observations_hdlr):

    G, b11, b12, b13, b22, b23, b33, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = observations_hdlr.postref_params

    #use III.4 in Winkler et al 1979 (A35; P901) for set of miller indices
    O = sqr(unit_cell((a,b,c,alpha,beta,gamma)).orthogonalization_matrix()).transpose()
    R = sqr(observations_hdlr.crystal_orientation.crystal_rotation_matrix()).transpose()
    CO = crystal_orientation(O*R, basis_type.direct)
    CO_rotate = CO.rotate_thru((1,0,0), rotx
                 ).rotate_thru((0,1,0), roty)
    A_star = sqr(CO_rotate.reciprocal_matrix())
    S0 = -1*col((0,0,1./observations_hdlr.wavelength))
    partiality_set = flex.double()
    delta_xy_set = flex.double()
    rs_set = flex.double()
    rh_set = flex.double()
    for miller_index, bragg_angle, alpha_angle, spot_pred_x_mm, spot_pred_y_mm in \
        zip(observations_hdlr.observations_original.indices(), \
        observations_hdlr.bragg_angle_set, observations_hdlr.alpha_angle_set, \
        observations_hdlr.spot_pred_x_mm_set, observations_hdlr.spot_pred_y_mm_set):
      if self.iparams.flag_beam_divergence:
        rs = math.sqrt((ry * math.cos(alpha_angle))**2 + (rz * math.sin(alpha_angle))**2) + \
          (r0 + (re*math.tan(bragg_angle)))
      else:
        rs = r0 + (re*math.tan(bragg_angle))
      h = col(miller_index)
      x = A_star * h
      S = x + S0
      rh = S.length() - (1/observations_hdlr.wavelength)

      spot_partiality = ((rs**2)/((2*(rh**2))+(rs**2)))
      partiality_set.append(spot_partiality)
      rs_set.append(rs)
      rh_set.append(rh)

      #finding coordinate x,y on the detector
      d_ratio = -observations_hdlr.detector_distance_mm/S[2]
      dx_mm = S[0]*d_ratio
      dy_mm = S[1]*d_ratio
      pred_xy = col((spot_pred_x_mm, spot_pred_y_mm))
      calc_xy = col((dx_mm, dy_mm))
      diff_xy = pred_xy - calc_xy
      delta_xy_set.append(diff_xy.length())

    return partiality_set, delta_xy_set, rs_set, rh_set
