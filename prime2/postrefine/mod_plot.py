from __future__ import division
from cctbx.array_family import flex
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


class plot_handler(object):
  """
  A wrapper class for plotting class
  """

  def __init__(self, iparams):
    """
    Intialitze parameters
    """
    self.iparams = iparams

  def plot_wilson(self, observations_hdlr_ref, observations_hdlr_set):
    obs_ref = observations_hdlr_ref.observations.deep_copy()
    binner = obs_ref.setup_binner(auto_binning=True)

    #collect mean intensity of the raw and scaled set
    I_mean_set = []
    I_full_mean_set = []
    for obs_hdlr in observations_hdlr_set:
      try:
        obs = obs_hdlr.observations.deep_copy()
        obs.use_binning_of(obs_ref)
        obs_mean = obs.mean(use_binning=True)
        I_mean = flex.double([0]*binner.n_bins_used())

        obs_full = obs_hdlr.calc_full_observations()
        obs_full.use_binning_of(obs_ref)
        obs_full_mean = obs_full.mean(use_binning=True)
        I_full_mean = flex.double([0]*binner.n_bins_used())
        one_dsqr = flex.double()
        for i_bin in range(binner.n_bins_used()):
          one_dsqr.append(1/binner.bin_d_range(i_bin+1)[1]**2)
          if obs_mean.data[i_bin+1] is not None:
            I_mean[i_bin] = obs_mean.data[i_bin+1]
          if obs_full_mean.data[i_bin+1] is not None:
            I_full_mean[i_bin] = obs_full_mean.data[i_bin+1]
        I_mean_set.append(I_mean)
        I_full_mean_set.append(I_full_mean)
      except Exception:
        pass

    I_mean_ref = flex.double(obs_ref.mean(
        use_binning=True,
        use_multiplicities=True).data[1:-1])

    plt.subplot(211)
    for I_mean in I_mean_set:
      plt.plot(one_dsqr, I_mean, linestyle='-', linewidth=2.0, c='b')
    plt.plot(one_dsqr, I_mean_ref, linestyle='-', linewidth=2.0, c='g')
    plt.title('Wilson plot')
    plt.xlabel('1/(d^2)')
    plt.ylabel('<I>')
    plt.grid()
    plt.subplot(212)
    for I_full_mean in I_full_mean_set:
      plt.plot(one_dsqr, I_full_mean, linestyle='-', linewidth=2.0, c='r')
    plt.plot(one_dsqr, I_mean_ref, linestyle='-', linewidth=2.0, c='g')
    plt.title('Wilson plot')
    plt.xlabel('1/(d^2)')
    plt.ylabel('<I>')
    plt.grid()
    plt.show()
