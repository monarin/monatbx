from iotbx import reflection_file_reader
from cctbx.array_family import flex
import sys
import cPickle as pickle
from mod_partiality import partiality_handler
import math

def read_input(args):
  hklref = ''
  data = ''
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklref':
      hklref = pair[1]
    elif pair[0]=='data':
      data = pair[1]
  if hklref == '':
    print "Please provide input hkl file."
    exit()

  return hklref, data


if (__name__ == "__main__"):
  #0 .read input parameters and frames (pickle files)
  hklref, data = read_input(args = sys.argv[1:])
  #read mtz ref. file
  reflection_file = reflection_file_reader.any_reflection_file(hklref)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array = miller_arrays[0]
  print "hklref"
  miller_array.show_summary()
  #open pickle file
  observations_pickle = pickle.load(open(data,"rb"))
  observations = observations_pickle["observations"][0]
  obs_asu = observations.map_to_asu().average_bijvoet_mates().as_amplitude_array()
  obs_asu_scaled = miller_array.scale(obs_asu, resolution_dependent=False)
  print "data"
  obs_asu_scaled.show_summary()
  #calculate r_factor by bin for scaled data
  print "scaled data R-factors"
  ref_match, obs_match = miller_array.common_sets(obs_asu_scaled,assert_is_similar_symmetry=False)
  ref_match.setup_binner(n_bins=20)
  r1_factor_bin = ref_match.r1_factor(obs_match, use_binning=True, assume_index_matching=True)
  r1_factor_bin.show()
  r1_factor = ref_match.r1_factor(obs_match, use_binning=False, assume_index_matching=True)
  print 'Overall R-factor:', r1_factor
  #calculate r_factor for partiality-corrected scaled data
  pixel_size_mm = 0.079346
  crystal_init_orientation = observations_pickle["current_orientation"][0]
  wavelength = observations_pickle["wavelength"]
  detector_distance_mm = observations_pickle['distance']
  mm_predictions = pixel_size_mm*(observations_pickle['mapped_predictions'][0])
  xbeam = observations_pickle["xbeam"]
  ybeam = observations_pickle["ybeam"]
  alpha_angle_obs = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
  spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
  spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])
  ph = partiality_handler()
  r0 = ph.calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                               observations.indices(), wavelength)
  ry, rz, re, voigt_nu, rotx, roty = (0, 0, 0.003, 0.5, 0, 0)
  partiality_init, delta_xy_init, rs_init, rh_init = ph.calc_partiality_anisotropy_set(\
                                                          crystal_init_orientation.unit_cell(),
                                                          rotx, roty, observations.indices(),
                                                          ry, rz, r0, re, voigt_nu,
                                                          two_theta, alpha_angle_obs, wavelength,
                                                          crystal_init_orientation,
                                                          spot_pred_x_mm, spot_pred_y_mm,
                                                          detector_distance_mm, "Voigt",
                                                          False)
  I_full = observations.data()/partiality_init
  sigI_full = observations.sigmas()/partiality_init
  obs_full = observations.customized_copy(data=I_full, sigmas=sigI_full)
  obs_full_asu = obs_full.map_to_asu().average_bijvoet_mates().as_amplitude_array()
  obs_asu_scaled = miller_array.scale(obs_asu, resolution_dependent=False)
