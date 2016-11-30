from __future__ import division
from cctbx.array_family import flex
from cctbx import miller
import numpy as np
import cPickle as pickle
from cctbx.crystal import symmetry
import math
from scitbx.matrix import sqr
from cctbx import statistics
from mod_observations import observations_handler
from mod_leastsqr import leastsqr_handler
from cctbx import crystal

class postrefine_handler(object):
  """
  handle post-refinement
  - read-in and store input in input_handler object
  - generate a mean-intensity-scaled mtz file as a reference set
  - perform post-refinement
  """
  def __init__(self, iparams):
    """
    Constructor
    """
    self.iparams = iparams
    self.scale_0_d_max = 5.0
    self.sigma_I_cutoff = 2.5

  def scale_0(self, filename):
    """
    Perform 0th-order scaling.
    """
    observations_hdlr = observations_handler(self.iparams)
    observations_hdlr.init_params_from_file(filename)
    """
    #get the miller_array_template
    uc = self.iparams.target_unit_cell.parameters()
    crystal_symmetry = crystal.symmetry(
        unit_cell=(uc[0], uc[1], uc[2], uc[3], uc[4], uc[5]),
        space_group_symbol=self.iparams.target_space_group
      )
    miller_array_complete = crystal_symmetry.build_miller_set(anomalous_flag=self.iparams.target_anomalous_flag,
                      d_min=self.iparams.d_min, d_max=self.iparams.d_max).array()
    miller_array_complete = miller_array_complete.set_observation_type_xray_intensity()
    miller_array_complete = miller_array_complete.customized_copy(data=flex.double([1e4]*len(miller_array_complete.indices())),
                      sigmas=flex.double([1e2]*len(miller_array_complete.indices())))
    asu_contents = {}
    asu_volume = miller_array_complete.unit_cell().volume()/float(miller_array_complete.space_group().order_z())
    number_carbons = asu_volume/18.0
    asu_contents.setdefault('C', number_carbons)
    ftmpl = miller_array_complete.as_amplitude_array()
    ftmpl.setup_binner(auto_binning=True)
    wp = statistics.wilson_plot(ftmpl, asu_contents, e_statistics=True)
    G, B = (wp.wilson_intensity_scale_factor*1e3,wp.wilson_b)

    stolsq = miller_array_complete.two_theta(observations_hdlr.wavelength).sin_theta_over_lambda_sq().data()
    I_wp = miller_array_complete.data()/(G * flex.exp(-2 * B * stolsq))
    sigI_wp = miller_array_complete.sigmas()/(G * flex.exp(-2 * B * stolsq))
    miller_array_complete = miller_array_complete.customized_copy(data=I_wp, sigmas=sigI_wp)
    miller_array_hires = miller_array_complete.resolution_filter(d_min=self.iparams.d_min, d_max=self.scale_0_d_max)
    miller_array_lowres = miller_array_complete.resolution_filter(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)

    obs_lowres = observations_hdlr.observations.resolution_filter(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)
    G_lowres = obs_lowres.mean()/miller_array_lowres.mean()
    observations_hdlr.set_params(G_lowres=G_lowres)
    """
    from mod_util import utility_handler
    util_hdlr = utility_handler()
    flag_hklisoin_found, miller_array_iso = util_hdlr.get_miller_array_from_mtz(self.iparams.hklisoin)
    miller_array_complete = miller_array_iso.as_intensity_array()
    miller_array_complete = miller_array_complete.customized_copy(data=miller_array_complete.data(), sigmas=miller_array_complete.sigmas())
    if self.iparams.target_anomalous_flag:
      miller_array_complete = miller_array_complete.generate_bijvoet_mates()
    miller_array_hires = miller_array_complete.resolution_filter(d_min=self.iparams.d_min, d_max=self.scale_0_d_max)
    miller_array_lowres = miller_array_complete.resolution_filter(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)

    observations_hdlr_ref = observations_handler(self.iparams)
    observations_hdlr_ref.init_params_from_merge(miller_array_hires)

    G,B=(1,0)
    #caculate G and B-tensor for high resolution reflections
    if True:
      observations_full = observations_hdlr.calc_full_observations()
      lsqr_hdlr = leastsqr_handler(self.iparams)
      lsqr_hdlr.optimize_scalefactors(observations_hdlr_ref, observations_hdlr,  use_binning=False)
    #except Exception:
    #  return None, ' {0:40} ==> SCALING FAILED'.format(observations_hdlr.pickle_filename_only)

    #calculate G for low-resolution
    obs_lowres = observations_hdlr.observations.resolution_filter(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)
    obs_lowres = obs_lowres.select(flex.bool(np.absolute((obs_lowres.data()-np.mean(obs_lowres.data()))/np.std(obs_lowres.data())) < self.sigma_I_cutoff))
    G_lowres = obs_lowres.mean()/miller_array_lowres.mean()
    observations_hdlr.set_params(G_lowres=G_lowres)

    txt_out = ' {0:40} ==> SCALED  RES:{1:5.2f} NREFL:{2:5d} SG:{12:8} Glow:{15:6.1f} G:{3:6.1f} B:{4:6.1f} {5:6.1f} {6:6.1f} {7:6.1f} {8:6.1f} {9:6.1f} Rfin:{11:5.2f} Gwp:{13:6.2f} Bwp{14:6.1f}'.format(observations_hdlr.pickle_filename_only, observations_hdlr.observations.d_min(), len(observations_hdlr.observations.data()), observations_hdlr.G, observations_hdlr.b11, observations_hdlr.b22, observations_hdlr.b33, observations_hdlr.b12, observations_hdlr.b13, observations_hdlr.b23, observations_hdlr.r_init, observations_hdlr.r_final, observations_hdlr.observations.space_group_info(), G, B, G_lowres)

    return observations_hdlr, txt_out

  def scale_with_reference(self, observations_hdlr_ref, filename=None, observations_hdlr=None, use_binning=False):
    """
    Scale each emage to the reference.
    """
    if True:
      if filename is not None:
        observations_hdlr = observations_handler(self.iparams)
        observations_hdlr.init_params_from_file(filename)

      #find G and B-tensor for high resolution reflections
      lsqr_hdlr = leastsqr_handler(self.iparams)
      lsqr_hdlr.optimize_scalefactors(observations_hdlr_ref, observations_hdlr,  use_binning=use_binning)

      #find G low resolution
      obs_ref_lowres = observations_hdlr_ref.observations.resolution_filter(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)
      obs_lowres = observations_hdlr.observations.resolution_filter(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)
      obs_lowres = obs_lowres.select(flex.bool(np.absolute((obs_lowres.data()-np.mean(obs_lowres.data()))/np.std(obs_lowres.data())) < self.sigma_I_cutoff))
      try:
        G_lowres = obs_lowres.mean()/obs_ref_lowres.mean()
      except Exception:
        G_lowres = observations_hdlr.G
      observations_hdlr.set_params(G_lowres=G_lowres)

      txt_out = ' {0:40} ==> SCALED  RES:{1:5.2f} NREFL:{2:5d} SG:{14:8} Glow:{15:6.1f} G:{3:6.1f} B:{4:6.1f} {5:6.1f} {6:6.1f} {7:6.1f} {8:6.1f} {9:6.1f} Rfin:{11:5.2f} CCini:{12:5.2f} CCfin:{13:5.2f}'.format(observations_hdlr.pickle_filename_only, observations_hdlr.observations.d_min(), len(observations_hdlr.observations.data()), observations_hdlr.G, observations_hdlr.b11, observations_hdlr.b22, observations_hdlr.b33, observations_hdlr.b12, observations_hdlr.b13, observations_hdlr.b23, observations_hdlr.r_init, observations_hdlr.r_final, observations_hdlr.cc_init, observations_hdlr.cc_final,observations_hdlr.observations.space_group_info(), G_lowres)
    #except Exception:
    #  return None, ' {0:40} ==> SCALING FAILED'.format(observations_hdlr.pickle_filename_only)

    return observations_hdlr, txt_out

  def postrefine_with_reference(self, observations_hdlr_ref, observations_hdlr):
    """
    Post-refine each image using a reference set.
    """
    msg = 'POSTREF'
    try:
      lsqr_hdlr = leastsqr_handler(self.iparams)
      lsqr_hdlr.optimize(observations_hdlr_ref, observations_hdlr)
    except Exception:
      msg = 'FAILED'

    txt_out = ' {0:40} ==> {1:7} RES:{2:5.2f} NREFL:{3:6d} r0{4:7.4f} re{5:7.4f} CELL:{6:5.1f} {7:5.1f} {8:5.1f} {9:5.1f} {10:5.1f} {11:5.1f} Rini:{12:5.2f} Rfin:{13:5.2f} CCini:{14:5.2f} CCfin:{15:5.2f}'.format(observations_hdlr.pickle_filename_only, msg, observations_hdlr.observations.d_min(), len(observations_hdlr.observations.data()), observations_hdlr.r0, observations_hdlr.re, observations_hdlr.uc_params[0],observations_hdlr.uc_params[1],observations_hdlr.uc_params[2],observations_hdlr.uc_params[3],observations_hdlr.uc_params[4],observations_hdlr.uc_params[5],observations_hdlr.r_init, observations_hdlr.r_final, observations_hdlr.cc_init, observations_hdlr.cc_final)

    return observations_hdlr, txt_out
