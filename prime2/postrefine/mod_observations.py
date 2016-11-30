from __future__ import division
from cctbx.array_family import flex
from cctbx.uctbx import unit_cell
import numpy as np
import cPickle as pickle
import math
from scitbx.matrix import sqr, col, row
from cctbx import miller
from cctbx import crystal
from mod_partiality import partiality_handler
from mod_polarity import polarity_handler
from cctbx import adptbx
from cctbx import sgtbx

class observations_handler(object):
  '''
  Author      : Uervirojnangkoorn, M.
  Created     : 8/2/2015
  A wrapper class to store refinement result
  '''

  def __init__(self, iparams):
    """
    initializing postref_params and stats
    """
    self.iparams = iparams
    self.postref_params = None
    self.stats = None
    self.detector_distance_mm = None
    self.alpha_angle_set = None
    self.bragg_angle_set = None
    self.spot_pred_x_mm_set = None
    self.spot_pred_y_mm_set = None
    self.rh_set = None
    self.partiality = None
    self.rs_set = None
    self.crystal_orientation = None
    self.G_lowres = 1
    self.scale_0_d_max = 5.0

  def refresh(self):
    """
    Calling refresh will update the intensity in the observations
    and reset the scale factors back to initial values.
    """
    if True:
      obs_scaled = self.calc_full_observations(flag_apply_partiality=False)

      self.observations = self.observations.customized_copy(data=obs_scaled.data(), sigmas=obs_scaled.sigmas())
      self.observations_original = self.observations_original.customized_copy(data=obs_scaled.data(), sigmas=obs_scaled.sigmas())
      self.postref_params[:7] = [1,0,0,0,0,0,0]
      self.G, self.b11, self.b12, self.b13, self.b22, self.b23, self.b33 = (1,0,0,0,0,0,0)
      self.G_lowres = 1.0

    #except Exception:
    #  self = None

  def calc_r0(self, miller_indices, wavelength, a_star_matrix):
    """
    calculate spot_radius based on rms delta_S for all spots
    """
    delta_S_all = flex.double()
    for miller_index in miller_indices:
      S0 = -1*col((0,0,1./wavelength))
      h = col(miller_index)
      x = a_star_matrix * h
      S = x + S0
      delta_S = S.length() - (1./wavelength)
      delta_S_all.append(delta_S)

    #spot_radius = math.sqrt(flex.mean(delta_S_all*delta_S_all))
    spot_radius = np.std(delta_S_all)
    return spot_radius

  def init_params_from_merge(self, miller_array_merge):
    observations = miller_array_merge.deep_copy()
    observations_original = miller_array_merge.deep_copy()
    a,b,c,alpha,beta,gamma = observations.unit_cell().parameters()
    postref_params = [1,0,0,0,0,0,0,0,0,0,0,0,self.iparams.gamma_e,a,b,c,alpha,beta,gamma]
    stats = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    #create observations handler object
    self.set_params(observations=observations,
      observations_original=observations_original,
      pickle_filename='new reference',
      postref_params=postref_params,
      stats=stats,
      G_lowres=1)
    self.init_type = 'from_merge'

  def init_params_from_file(self, filename):
    """
    From a giving observation file, build observations_handler object.
    """
    self.init_type = 'from_file'
    observations_pickle = pickle.load(open(filename,"rb"))
    observations_original = observations_pickle["observations"][0]

    try:
      #apply constrain using the crystal system
      uc = observations_original.unit_cell().parameters()
      crystal_symmetry = crystal.symmetry(
          unit_cell=(uc[0], uc[1], uc[2], uc[3], uc[4], uc[5]),
          space_group_symbol=self.iparams.target_space_group
        )
      observations_original = observations_original.customized_copy(anomalous_flag=self.iparams.target_anomalous_flag,
                      crystal_symmetry=crystal_symmetry)
      miller_array_complete = crystal_symmetry.build_miller_set(anomalous_flag=self.iparams.target_anomalous_flag,
                      d_min=self.iparams.d_min, d_max=self.iparams.d_max).array()
    except Exception:
      txt_exception = 'rejected (bad cell) CELL:{0:6.2f} {1:6.2f} {2:6.2f} {3:6.2f} {4:6.2f} {5:6.2f}'.format(uc[0], uc[1], uc[2], uc[3], uc[4], uc[5])
      return False, txt_exception

    if len(observations_original.indices()) != len(observations_pickle['mapped_predictions'][0]):
      txt_exception = 'number of reflections does not match with number of predicted locations.'
      return False, txt_exception

    #handle polarity problem (if exisited)
    polar_hdlr = polarity_handler(self.iparams)
    observations_non_polar = polar_hdlr.get_observations_non_polar(observations_original, filename)

    #initialize partiality (postref) parameters
    a_star_matrix = sqr(observations_pickle["current_orientation"][0].reciprocal_matrix())
    r0 = self.calc_r0(observations_original.indices(), observations_pickle["wavelength"], a_star_matrix)
    a,b,c,alpha,beta,gamma = observations_original.unit_cell().parameters()
    postref_params = [1,0,0,0,0,0,0,0,0,0,0,r0,self.iparams.gamma_e,a,b,c,alpha,beta,gamma]
    stats = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

    #initialize lab coordinate parameters
    detector_distance_mm = observations_pickle['distance']
    mm_predictions = self.iparams.pixel_size_mm*(observations_pickle['mapped_predictions'][0])
    xbeam = observations_pickle["xbeam"]
    ybeam = observations_pickle["ybeam"]
    alpha_angle_set = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
    spot_pred_x_mm_set = flex.double([pred[0]-xbeam for pred in mm_predictions])
    spot_pred_y_mm_set = flex.double([pred[1]-ybeam for pred in mm_predictions])
    bragg_angle_set = observations_original.two_theta(observations_pickle["wavelength"]).data()

    #filter resolution
    i_sel_res = observations_original.resolution_filter_selection(d_min=self.iparams.d_min, d_max=self.iparams.d_max)
    observations_original = observations_original.select(i_sel_res)
    observations_non_polar = observations_non_polar.select(i_sel_res)
    alpha_angle_set = alpha_angle_set.select(i_sel_res)
    bragg_angle_set = bragg_angle_set.select(i_sel_res)
    spot_pred_x_mm_set = spot_pred_x_mm_set.select(i_sel_res)
    spot_pred_y_mm_set = spot_pred_y_mm_set.select(i_sel_res)

    #Filter weak
    i_sel = (observations_original.data()/observations_original.sigmas()) > 2
    observations_original = observations_original.select(i_sel)
    observations_non_polar = observations_non_polar.select(i_sel)
    alpha_angle_set = alpha_angle_set.select(i_sel)
    bragg_angle_set = bragg_angle_set.select(i_sel)
    spot_pred_x_mm_set = spot_pred_x_mm_set.select(i_sel)
    spot_pred_y_mm_set = spot_pred_y_mm_set.select(i_sel)

    #merge equivalents
    observations_non_polar_merge = observations_non_polar.merge_equivalents().array()
    ma_seq = observations_non_polar.customized_copy(data=flex.int(range(len(observations_non_polar.indices()))))
    ma_seq_match = ma_seq.common_set(observations_non_polar_merge)
    i_sel = flex.bool([False]*len(observations_non_polar.indices()))
    for seq in ma_seq_match.data(): i_sel[seq] = True

    observations_non_polar = observations_non_polar.select(i_sel)
    observations_original = observations_original.select(i_sel)
    alpha_angle_set = alpha_angle_set.select(i_sel)
    bragg_angle_set = bragg_angle_set.select(i_sel)
    spot_pred_x_mm_set = spot_pred_x_mm_set.select(i_sel)
    spot_pred_y_mm_set = spot_pred_y_mm_set.select(i_sel)

    #create observations handler object
    self.set_params(observations=observations_non_polar,
      observations_original=observations_original,
      G_lowres=1,
      postref_params=postref_params,
      stats=stats,
      bragg_angle_set=bragg_angle_set,
      alpha_angle_set=alpha_angle_set,
      spot_pred_x_mm_set=spot_pred_x_mm_set,
      spot_pred_y_mm_set=spot_pred_y_mm_set,
      wavelength=observations_pickle["wavelength"],
      crystal_orientation=observations_pickle["current_orientation"][0],
      detector_distance_mm=detector_distance_mm,
      pickle_filename=filename)


    return True, 'ok'


  def set_params(self, observations=None,
      observations_original=None,
      G_lowres=None,
      postref_params=None,
      partiality=None,
      rh_set=None,
      rs_set=None,
      stats=None,
      bragg_angle_set=None,
      alpha_angle_set=None,
      spot_pred_x_mm_set=None,
      spot_pred_y_mm_set=None,
      wavelength=None,
      crystal_orientation=None,
      detector_distance_mm=None,
      pickle_filename=None):
    """
    reset observations and postrefinement parameters
    """
    if observations is not None:
      self.observations = observations.deep_copy()
    if observations_original is not None:
      self.observations_original = observations_original.deep_copy()
    if G_lowres is not None:
      self.G_lowres = G_lowres
    if postref_params is not None:
      self.postref_params = postref_params[:]
      self.G,self.b11,self.b22,self.b33,self.b12,self.b13,self.b23,self.rotx,self.roty,self.ry,self.rz,self.r0,self.re = postref_params[:13]
      self.uc_params = flex.double(postref_params[13:])
      self.unit_cell = unit_cell((self.uc_params[0], self.uc_params[1], self.uc_params[2],
            self.uc_params[3], self.uc_params[4], self.uc_params[5]))
    if partiality is not None:
      self.partiality = partiality[:]
    if rh_set is not None:
      self.rh_set = rh_set[:]
    if rs_set is not None:
      self.rs_set = rs_set[:]
    if stats is not None:
      self.stats = stats[:]
      self.cc_init = stats[0]
      self.cc_final = stats[1]
      self.r_init = stats[2]
      self.r_final = stats[3]
      self.r_xy_init = stats[4]
      self.r_xy_final = stats[5]
      self.cc_iso_init = stats[6]
      self.cc_iso_final = stats[7]

      if self.r_init == 0:
        self.r_change = 0
      else:
        self.r_change = ((self.r_final - self.r_init)/self.r_init)*100

      if self.cc_init == 0:
        self.cc_change = 0
      else:
        self.cc_change = ((self.cc_final - self.cc_init)/self.cc_init)*100
    if bragg_angle_set is not None:
      self.bragg_angle_set = bragg_angle_set[:]
    if alpha_angle_set is not None:
      self.alpha_angle_set = alpha_angle_set[:]
    if spot_pred_x_mm_set is not None:
      self.spot_pred_x_mm_set = spot_pred_x_mm_set[:]
    if spot_pred_y_mm_set is not None:
      self.spot_pred_y_mm_set = spot_pred_y_mm_set[:]
    if wavelength is not None:
      self.wavelength = wavelength
    if crystal_orientation is not None:
      self.crystal_orientation = crystal_orientation
    if detector_distance_mm is not None:
      self.detector_distance_mm = detector_distance_mm
    if pickle_filename is not None:
      self.pickle_filename = pickle_filename
      pickle_filepaths = pickle_filename.split('/')
      if len(pickle_filepaths) > 0:
        self.pickle_filename_only = pickle_filepaths[len(pickle_filepaths)-1]
      else:
        self.pickle_filename_only = pickle_filename

    #set partiality
    if self.postref_params is not None and self.partiality is None and self.crystal_orientation is not None:
      partiality_hdlr = partiality_handler(self.iparams)
      self.partiality, dum0, self.rh_set, self.rs_set  = partiality_hdlr.calc_partiality_anisotropy_set(self)

    #get asu_contents
    self.asu_contents = {}
    asu_volume = self.observations.unit_cell().volume()/float(self.observations.space_group().order_z())
    number_carbons = asu_volume/18.0
    self.asu_contents.setdefault('C', number_carbons)


  def copy_from(self, observations_hdlr):
    """
    copy all the attribuites from the given observations_hdlr
    """
    if observations_hdlr.init_type == 'from_file':
      self.set_params(observations=observations_hdlr.observations.deep_copy(),
        observations_original=observations_hdlr.observations_original.deep_copy(),
        G_lowres = observations_hdlr.G_lowres,
        postref_params = observations_hdlr.postref_params[:],
        partiality = observations_hdlr.partiality[:],
        rh_set = observations_hdlr.rh_set[:],
        rs_set = observations_hdlr.rs_set[:],
        stats = observations_hdlr.stats[:],
        bragg_angle_set = observations_hdlr.bragg_angle_set[:],
        alpha_angle_set = observations_hdlr.alpha_angle_set[:],
        spot_pred_x_mm_set = observations_hdlr.spot_pred_x_mm_set[:],
        spot_pred_y_mm_set = observations_hdlr.spot_pred_y_mm_set[:],
        wavelength = observations_hdlr.wavelength,
        crystal_orientation = observations_hdlr.crystal_orientation,
        detector_distance_mm = observations_hdlr.detector_distance_mm,
        pickle_filename = observations_hdlr.pickle_filename)
      self.init_type = 'from_file'
    elif observations_hdlr.init_type == 'from_merge':
      self.set_params(observations=observations_hdlr.observations.deep_copy(),
        observations_original=observations_hdlr.observations_original.deep_copy(),
        pickle_filename = observations_hdlr.pickle_filename)
      a,b,c,alpha,beta,gamma = observations_hdlr.observations.unit_cell().parameters()
      self.uc_params = flex.double([a,b,c,alpha,beta,gamma])
      self.init_type = 'from_merge'

  def calc_debye_waller_factor(self, observations, b_cart):
    #apply space-group constraints on B-tensor
    space_group = self.observations.space_group()
    adp_constraints = sgtbx.tensor_rank_2_constraints(space_group=space_group,
              reciprocal_space=True)
    #convert B tensor direct-space to reciprocal space
    fmx = sqr(self.unit_cell.fractionalization_matrix())
    u_cart = b_cart / (8 * math.pi**2)
    u_star_mat = fmx * u_cart * fmx.transpose()
    u_star = u_star_mat.as_sym_mat3()
    u_star = space_group.average_u_star(u_star)
    u_indep = adp_constraints.independent_params(all_params=u_star)
    u_star = adp_constraints.all_params(independent_params=u_indep)
    dw_c = adptbx.debye_waller_factor_u_star(observations.indices(), u_star)

    return dw_c

  def calc_full_observations(self, observations_type='asymmetric_unit', flag_apply_partiality=False):
    """
    calculate partiality and full observations (miller array)
    """
    if True:
      if self.init_type == 'from_file':
        partiality_hdlr = partiality_handler(self.iparams)
        self.partiality, dum0, self.rh_set, self.rs_set  = partiality_hdlr.calc_partiality_anisotropy_set(self)

        b_cart = sqr((self.b11, self.b12, self.b13, self.b12, self.b22, self.b23, self.b13, self.b23, self.b33))
        dw = self.calc_debye_waller_factor(self.observations, b_cart)

        if flag_apply_partiality:
          partiality = self.partiality[:]
        else:
          partiality = flex.double([1]*len(self.observations.data()))

        G_set = flex.double([self.G]*len(self.observations.data()))
        i_lowres = self.observations.resolution_filter_selection(d_min=self.scale_0_d_max, d_max=self.iparams.d_max)
        G_set.set_selected(i_lowres, self.G_lowres)
        dw.set_selected(i_lowres, 1.0)

        I_full = flex.double(self.observations.data()/(G_set * dw * partiality))
        sigI_full = flex.double(self.observations.sigmas()/(G_set * dw * partiality))

      elif self.init_type == 'from_merge':
        I_full = self.observations.data()[:]
        sigI_full = self.observations.sigmas()[:]

      if observations_type=='asymmetric_unit':
        observations_full = self.observations.customized_copy(data=I_full, sigmas=sigI_full)
      elif observations_type=='original':
        observations_full = self.observations_original.customized_copy(data=I_full, sigmas=sigI_full)
      return observations_full
    #except Exception:
    #  return None

  def calc_full_intensity(self):
    """
    calculate partiality and full intensity (just I)
    """
    if self.init_type == 'from_file':
      return self.calc_full_observations.data()
    elif self.init_type == 'from_merge':
      return self.observations.data()[:]

  def resolution_filter(self, d_min=0, d_max=99):
    i_sel = self.observations_original.resolution_filter_selection(d_min=d_min, d_max=d_max)

    self.observations = self.observations.select(i_sel)
    self.observations_original = self.observations_original.select(i_sel)
    if self.partiality is not None:
      self.partiality = self.partiality.select(i_sel)
      self.rs_set = self.rs_set.select(i_sel)
      self.rh_set = self.rh_set.select(i_sel)
      self.bragg_angle_set = self.bragg_angle_set.select(i_sel)
      self.alpha_angle_set = self.alpha_angle_set.select(i_sel)
      self.spot_pred_x_mm_set = self.spot_pred_x_mm_set.select(i_sel)
      self.spot_pred_y_mm_set = self.spot_pred_y_mm_set.select(i_sel)

  def outlier_filter(self, sigma_I_cutoff=2.5):
    i_sel = flex.bool(np.absolute((self.observations.data()-np.mean(self.observations.data()))/np.std(self.observations.data())) < sigma_I_cutoff)

    self.observations = self.observations.select(i_sel)
    self.observations_original = self.observations_original.select(i_sel)
    if self.partiality is not None:
      self.partiality = self.partiality.select(i_sel)
      self.rs_set = self.rs_set.select(i_sel)
      self.rh_set = self.rh_set.select(i_sel)
      self.bragg_angle_set = self.bragg_angle_set.select(i_sel)
      self.alpha_angle_set = self.alpha_angle_set.select(i_sel)
      self.spot_pred_x_mm_set = self.spot_pred_x_mm_set.select(i_sel)
      self.spot_pred_y_mm_set = self.spot_pred_y_mm_set.select(i_sel)

  def common_filter(self, observations_hdlr_other):

    obs_hdlr = observations_handler(self.iparams)
    obs_hdlr.copy_from(self)

    obs_hdlr_other = observations_handler(self.iparams)
    obs_hdlr_other.copy_from(observations_hdlr_other)

    obs_match, obs_other_match = obs_hdlr.observations.common_sets(obs_hdlr_other.observations,assert_is_similar_symmetry=False)

    ma_obs_original = obs_hdlr.observations.customized_copy(data=obs_hdlr.observations_original.indices())
    ma_partiality = obs_hdlr.observations.customized_copy(data=obs_hdlr.partiality)
    ma_rs_set = obs_hdlr.observations.customized_copy(data=obs_hdlr.rs_set)
    ma_rh_set = obs_hdlr.observations.customized_copy(data=obs_hdlr.rh_set)
    ma_bragg_angle_set = obs_hdlr.observations.customized_copy(data=obs_hdlr.bragg_angle_set)
    ma_alpha_angle_set = obs_hdlr.observations.customized_copy(data=obs_hdlr.alpha_angle_set)
    ma_spot_pred_x_mm_set = obs_hdlr.observations.customized_copy(data=obs_hdlr.spot_pred_x_mm_set)
    ma_spot_pred_y_mm_set = obs_hdlr.observations.customized_copy(data=obs_hdlr.spot_pred_y_mm_set)

    ma_obs_original_other = obs_hdlr.observations.customized_copy(data=obs_hdlr_other.observations_original.indices())
    ma_partiality_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.partiality)
    ma_rs_set_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.rs_set)
    ma_rh_set_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.rh_set)
    ma_bragg_angle_set_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.bragg_angle_set)
    ma_alpha_angle_set_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.alpha_angle_set)
    ma_spot_pred_x_mm_set_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.spot_pred_x_mm_set)
    ma_spot_pred_y_mm_set_other = obs_hdlr_other.observations.customized_copy(data=obs_hdlr_other.spot_pred_y_mm_set)

    obs_match, obs_other_match = obs_hdlr.observations.common_sets(obs_hdlr_other.observations,assert_is_similar_symmetry=False)
    ma_partiality, ma_partiality_other = ma_partiality.common_sets(ma_partiality_other,assert_is_similar_symmetry=False)
    ma_rs_set, ma_rs_set_other = ma_rs_set.common_sets(ma_rs_set_other,assert_is_similar_symmetry=False)
    ma_rh_set, ma_rh_set_other = ma_rh_set.common_sets(ma_rh_set_other,assert_is_similar_symmetry=False)
    ma_bragg_angle_set, ma_bragg_angle_set_other = ma_bragg_angle_set.common_sets(ma_bragg_angle_set,assert_is_similar_symmetry=False)
    ma_alpha_angle_set, ma_alpha_angle_set_other = ma_alpha_angle_set.common_sets(ma_alpha_angle_set_other,assert_is_similar_symmetry=False)
    ma_spot_pred_x_mm_set, ma_spot_pred_x_mm_set_other = ma_spot_pred_x_mm_set.common_sets(ma_spot_pred_x_mm_set_other,assert_is_similar_symmetry=False)
    ma_spot_pred_y_mm_set, ma_spot_pred_y_mm_set_other = ma_spot_pred_y_mm_set.common_sets(ma_spot_pred_y_mm_set_other,assert_is_similar_symmetry=False)


    matches = miller.match_multi_indices(
                  miller_indices_unique=obs_hdlr.observations.indices(),
                  miller_indices=obs_match.indices())
    miller_indices_original = flex.miller_index([obs_hdlr.observations_original.indices()[pair[0]] for pair in matches.pairs()])

    matches = miller.match_multi_indices(
                  miller_indices_unique=obs_hdlr_other.observations.indices(),
                  miller_indices=obs_match.indices())
    miller_indices_original_other = flex.miller_index([obs_hdlr_other.observations_original.indices()[pair[0]] for pair in matches.pairs()])

    obs_original_match = obs_match.customized_copy(indices=miller_indices_original)
    obs_original_other_match = obs_other_match.customized_copy(indices=miller_indices_original_other)

    if obs_hdlr.partiality is not None:
      obs_hdlr.set_params(observations=obs_match,
          observations_original=obs_original_match,
          partiality=ma_partiality.data(),
          rs_set=ma_rs_set.data(),
          rh_set=ma_rh_set.data(),
          bragg_angle_set=ma_bragg_angle_set.data(),
          alpha_angle_set=ma_alpha_angle_set.data(),
          spot_pred_x_mm_set=ma_spot_pred_x_mm_set.data(),
          spot_pred_y_mm_set=ma_spot_pred_y_mm_set.data())
      obs_hdlr_other.set_params(observations=obs_other_match,
          observations_original=obs_original_other_match,
          partiality=ma_partiality_other.data(),
          rs_set=ma_rs_set_other.data(),
          rh_set=ma_rh_set_other.data(),
          bragg_angle_set=ma_bragg_angle_set_other.data(),
          alpha_angle_set=ma_alpha_angle_set_other.data(),
          spot_pred_x_mm_set=ma_spot_pred_x_mm_set_other.data(),
          spot_pred_y_mm_set=ma_spot_pred_y_mm_set_other.data())
    else:
      obs_hdlr.set_params(observations=obs_match,
          observations_original=obs_original_match)
      obs_hdlr_other.set_params(observations=obs_other_match,
          observations_original=obs_original_other_match)


    return obs_hdlr, obs_hdlr_other
