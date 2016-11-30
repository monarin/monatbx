from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import cPickle as pickle
from sets import Set
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal
from mod_polarity import polarity_handler
from mod_postrefine import postrefine_handler

class cluster_handler(object):
  '''
  A wrapper class for organizer class
  '''

  def __init__(self, iparams):
    """
    Intialitze parameters
    """
    self.iparams = iparams

  def get_dict(self):
    """
    Generate a complete miller array for a given unit-cell and space group and
    return a dictionary with miller indices as key.
    """
    uc = self.iparams.target_unit_cell.parameters()
    crystal_symmetry = crystal.symmetry(
        unit_cell=(uc[0], uc[1], uc[2], uc[3], uc[4], uc[5]),
        space_group_symbol=self.iparams.target_space_group)
    miller_set = crystal_symmetry.build_miller_set(
        anomalous_flag=self.iparams.target_anomalous_flag,
        d_min=self.iparams.d_min, d_max=self.iparams.d_max)
    miller_array = miller_set.array()

    ma_complete_dict = {}
    for miller_index in miller_array.indices():
      ma_complete_dict[miller_index]=[]

    return ma_complete_dict
