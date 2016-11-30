from __future__ import division
from cctbx.array_family import flex
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from cctbx import miller
from cctbx import crystal
import cPickle as pickle
from sets import Set
from cctbx import statistics
from mod_observations import observations_handler

class merge_handler(object):
  '''
  A wrapper class for organizer class
  '''

  def __init__(self, iparams):
    """
    Intialitze parameters
    """
    self.iparams = iparams
    self.sigma_I_cutoff = 2.5

  def calc_mean_unit_cell(self, obs_hdlr_set):
    a_all = flex.double()
    b_all = flex.double()
    c_all = flex.double()
    alpha_all = flex.double()
    beta_all = flex.double()
    gamma_all = flex.double()
    for obs_hdlr in obs_hdlr_set:
      a_all.append(obs_hdlr.uc_params[0])
      b_all.append(obs_hdlr.uc_params[1])
      c_all.append(obs_hdlr.uc_params[2])
      alpha_all.append(obs_hdlr.uc_params[3])
      beta_all.append(obs_hdlr.uc_params[4])
      gamma_all.append(obs_hdlr.uc_params[5])

    uc_mean = flex.double([np.mean(a_all), np.mean(b_all), np.mean(c_all), np.mean(alpha_all), np.mean(beta_all), np.mean(gamma_all)])
    uc_med = flex.double([np.median(a_all), np.median(b_all), np.median(c_all), np.median(alpha_all), np.median(beta_all), np.median(gamma_all)])
    uc_std = flex.double([np.std(a_all), np.std(b_all), np.std(c_all), np.std(alpha_all), np.std(beta_all), np.std(gamma_all)])

    return uc_mean, uc_med, uc_std

  def merge_observations_handler(self, obs_hdlr_ref, obs_hdlr, r_final_cutoff=1.0):
    """
    Take the reference (full observations) and merge with the given
    observations (need to correct for scale and partiality).
    """
    if obs_hdlr.r_final < r_final_cutoff:
      obs_ref = obs_hdlr_ref.observations.deep_copy()
      obs = obs_hdlr.calc_full_observations()
      miller_indices_all = obs_ref.indices()
      I_all = obs_ref.data()
      sigI_all = obs_ref.sigmas()
      miller_indices_all.extend(obs.indices())
      I_all.extend(obs.data())
      sigI_all.extend(obs.sigmas())

      uc_mean, uc_med, uc_std = self.calc_mean_unit_cell([obs_hdlr_ref, obs_hdlr])
      crystal_symmetry = crystal.symmetry(
          unit_cell=(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]),
          space_group_symbol=self.iparams.target_space_group)
      miller_set=miller.set(
                  crystal_symmetry=crystal_symmetry,
                  indices=miller_indices_all,
                  anomalous_flag=self.iparams.target_anomalous_flag)
      miller_array_all = miller_set.array(
                data=I_all,
                sigmas=sigI_all).set_observation_type_xray_intensity()

      me_obj = miller_array_all.merge_equivalents()

      #build newly merged observations_handler
      observations_hdlr_merge = observations_handler(self.iparams)
      observations_hdlr_merge.init_params_from_merge(me_obj.array())

      txt_out = ' {0:40} ==> MERGED  RES:{1:5.2f} NREFL:{2:5d} R_MEAS:{3:6.2f} R_MERGE:{4:6.2f}'.format(observations_hdlr_merge.pickle_filename_only, observations_hdlr_merge.observations.d_min(), len(observations_hdlr_merge.observations.data()), me_obj.r_meas(), me_obj.r_merge())
    else:
      observations_hdlr_merge = observations_handler(self.iparams)
      observations_hdlr_merge.copy_from(obs_hdlr_ref)
      txt_out = ' {0:40} ==> NOT MERGED (Rfin too low)'.format(observations_hdlr_merge.pickle_filename_only)
    return observations_hdlr_merge, txt_out


  def merge_multi_observations_handler(self, observations_hdlr_set, flag_show_summary=False, r_final_cutoff=1.0):
    """
    Merge all observations_hdlr in the set
    """
    miller_indices_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()

    cn_good_observations = 0
    observations_hdlr_selected_set = []
    for i_frame in range(len(observations_hdlr_set)):
      if observations_hdlr_set[i_frame].r_final < r_final_cutoff:
        obs_full = observations_hdlr_set[i_frame].calc_full_observations()
        if obs_full is not None:
          cn_good_observations += 1
          observations_hdlr_selected_set.append(observations_hdlr_set[i_frame])
          miller_indices_all.extend(obs_full.indices())
          I_all.extend(obs_full.data())
          sigI_all.extend(obs_full.sigmas())

    uc_mean, uc_med, uc_std = self.calc_mean_unit_cell(observations_hdlr_selected_set)
    try:
      crystal_symmetry = crystal.symmetry(
          unit_cell=(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]),
          space_group_symbol=self.iparams.target_space_group)
    except Exception:
      print uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]
      uc_mean = self.iparams.target_unit_cell.parameters()
      crystal_symmetry = crystal.symmetry(
          unit_cell=(uc_mean[0], uc_mean[1], uc_mean[2], uc_mean[3], uc_mean[4], uc_mean[5]),
          space_group_symbol=self.iparams.target_space_group)
    miller_set=miller.set(
                crystal_symmetry=crystal_symmetry,
                indices=miller_indices_all,
                anomalous_flag=self.iparams.target_anomalous_flag)
    miller_array_all = miller_set.array(
              data=I_all,
              sigmas=sigI_all).set_observation_type_xray_intensity()

    if miller_array_all.size() == 0:
      txt_out = ' FINAL MERGE FAILED'
      return None, None, txt_out

    #remove outlier reflections
    perm = miller_array_all.sort_permutation(by_value="packed_indices")
    miller_array_all_sort = miller_array_all.select(perm)
    current_miller_index = flex.miller_index([miller_array_all_sort.indices()[0]])
    group_id = 0
    group_id_set = flex.int([0]*len(miller_array_all_sort.indices()))
    for miller_index, i_seq in zip(miller_array_all_sort.indices(), range(len(miller_array_all_sort.indices()))):
      if current_miller_index[0][0]==miller_index[0] and current_miller_index[0][1]==miller_index[1] and current_miller_index[0][2]==miller_index[2]:
        group_id_set[i_seq] = group_id
      else:
        current_miller_index = flex.miller_index([miller_array_all_sort.indices()[i_seq]])
        group_id += 1
        group_id_set[i_seq] = group_id

    miller_indices_all = flex.miller_index()
    I_all = flex.double()
    sigI_all = flex.double()
    for group_id in range(max(group_id_set)+1):
      miller_array_group = miller_array_all_sort.select(group_id_set==group_id)
      mean_I_over_sigI = flex.mean(miller_array_group.data()/miller_array_group.sigmas())
      std_I_over_sigI = math.sqrt(mean_I_over_sigI)
      #miller_array_group = miller_array_group.select(flex.bool(np.absolute((miller_array_group.data()-np.mean(miller_array_group.data()))/np.std(miller_array_group.data())) < self.sigma_I_cutoff))
      miller_array_group = miller_array_group.select(flex.bool((miller_array_group.data()/miller_array_group.sigmas()) > mean_I_over_sigI + (2*std_I_over_sigI)))
      miller_indices_all.extend(miller_array_group.indices())
      I_all.extend(miller_array_group.data())
      sigI_all.extend(miller_array_group.sigmas())
    miller_set=miller.set(
                crystal_symmetry=crystal_symmetry,
                indices=miller_indices_all,
                anomalous_flag=self.iparams.target_anomalous_flag)
    miller_array_all = miller_set.array(
              data=I_all,
              sigmas=sigI_all).set_observation_type_xray_intensity()

    miller_array_all = miller_array_all.resolution_filter(d_min=self.iparams.d_min, d_max=self.iparams.d_max)
    #merge
    me_obj = miller_array_all.merge_equivalents()
    miller_array_merge = me_obj.array().as_intensity_array()
    if miller_array_merge.size() == 0:
      txt_out = ' FINAL MERGE FAILED'
      return None, None, txt_out

    miller_array_merge.setup_binner(n_bins=self.iparams.n_bins)

    #get asu_contents
    asu_contents = {}
    asu_volume = miller_array_merge.unit_cell().volume()/float(miller_array_merge.space_group().order_z())
    number_carbons = asu_volume/18.0
    asu_contents.setdefault('C', number_carbons)



    #build newly merged observations_handler
    observations_hdlr_merge = observations_handler(self.iparams)
    observations_hdlr_merge.init_params_from_merge(miller_array_merge)
    #write as mtz file
    mtz_dataset_merge = miller_array_merge.as_mtz_dataset(column_root_label="IOBS")
    mtz_dataset_merge.mtz_object().write(file_name=self.iparams.run_no+'/postref_merged.mtz')

    txt_out = ' FINAL MERGE RES:{0:5.2f} NREFL:{1:5d} NIMG:{2:6d} R_MEAS:{3:6.2f} R_MERGE:{4:6.2f}'.format(observations_hdlr_merge.observations.d_min(), len(observations_hdlr_merge.observations.data()), cn_good_observations, me_obj.r_meas(), me_obj.r_merge())
    if flag_show_summary:
      print 'SUMMARY: mean intensity (%6.0f images)'%(cn_good_observations)
      ma_mean_I = miller_array_merge.mean(use_binning=True)
      ma_mean_I.show()
      print ''

      from mod_util import utility_handler
      util_hdlr = utility_handler()
      flag_hklisoin_found, miller_array_iso = util_hdlr.get_miller_array_from_mtz(self.iparams.hklisoin)
      if flag_hklisoin_found:
        try:
          print 'SUMMARY: CCiso'
          ma_cc_iso = miller_array_merge.correlation(miller_array_iso, use_binning=True, assert_is_similar_symmetry=False)
          ma_cc_iso.show()
          print ''
        except Exception:
          pass
      print 'SUMMARY: Merging statistics'
      me_obj.show_summary()
      print ''

      if self.iparams.target_anomalous_flag:
        miller_array_all.setup_binner(n_bins=self.iparams.n_bins)
        ma_cc_anom = miller_array_all.cc_anom(use_binning=True)
        ma_cc_anom.show()

      from mod_util import utility_handler
      util_hdlr = utility_handler()
      flag_hklisoin_found, miller_array_iso = util_hdlr.get_miller_array_from_mtz(self.iparams.hklisoin)
      miller_array_iso = miller_array_iso.as_intensity_array()
      if self.iparams.target_anomalous_flag:
        miller_array_iso = miller_array_iso.generate_bijvoet_mates()
      observations_hdlr_iso = observations_handler(self.iparams)
      observations_hdlr_iso.init_params_from_merge(miller_array_iso)

      from mod_plot import plot_handler
      plot_hdlr = plot_handler(self.iparams)
      plot_hdlr.plot_wilson(observations_hdlr_iso, observations_hdlr_selected_set)

    return observations_hdlr_merge, observations_hdlr_selected_set, txt_out
