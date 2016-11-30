'''
Fit second amplitude set to the first set.
Report CC, K, and B
'''

from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller
import sys
import matplotlib.pyplot as plt
import math

def read_input(args):
  hkla = None
  hklb = ''
  d_min = 0
  d_max = 99
  n_bins = 20
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hkla':
      hkla = pair[1]
    elif pair[0]=='hklb':
      hklb = pair[1]
    elif pair[0]=='d_min':
      d_min = float(pair[1])
    elif pair[0]=='d_max':
      d_max = float(pair[1])
    elif pair[0]=='n_bins':
      n_bins = int(pair[1])

  if hklb == '':
    print "Please provide input hkl files."
    exit()

  return hkla, hklb, d_min, d_max, n_bins


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hkla, hklb, d_min, d_max, n_bins = read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hkla)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array_a = miller_arrays[0]
  miller_array_a = miller_array_a.resolution_filter(d_min=d_min, d_max=d_max)
  print 'First reflection file:', hkla
  miller_array_a.show_summary()
  miller_array_a.setup_binner(n_bins=n_bins)
  f_ref_selected = miller_array_a.select(miller_array_a.data() > 0)
  f_ref_selected.use_binning_of(miller_array_a)
  expected_f_sq = flex.double(f_ref_selected.mean_sq(
      use_binning=True,
      use_multiplicities=True).data[1:-1])

  reflection_file = reflection_file_reader.any_reflection_file(hklb)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array_b = miller_arrays[0]
  miller_array_b = miller_array_b.resolution_filter(d_min=d_min, d_max=d_max)
  print 'Second reflection file:', hklb
  miller_array_b.show_summary()

  """
  f_obs_selected = miller_array_b.select(miller_array_b.data() > 0)
  f_obs_selected.use_binning_of(miller_array_a)
  mean_fobs_sq = flex.double(f_obs_selected.mean_sq(
      use_binning=True,
      use_multiplicities=True).data[1:-1])

  # compute <s^2> = <(sin(theta)/lambda)^2> in resolution shells
  stol_sq = f_obs_selected.sin_theta_over_lambda_sq()
  stol_sq.use_binner_of(f_obs_selected)
  mean_stol_sq = flex.double(stol_sq.mean(
      use_binning=True,
      use_multiplicities=True).data[1:-1])

  # fit to straight line
  x = mean_stol_sq
  y = flex.log(mean_fobs_sq / expected_f_sq)
  fit = flex.linear_regression(x, y)
  assert fit.is_well_defined()
  fit_y_intercept = fit.y_intercept()
  fit_slope = fit.slope()
  wilson_intensity_scale_factor = math.exp(fit_y_intercept) # intensity scale factor
  wilson_k = math.sqrt(wilson_intensity_scale_factor) # conversion to amplitude scale factor
  wilson_b = -fit_slope / 2
  fit_correlation = flex.linear_correlation(x, y).coefficient()
  print 'G=%6.2f B=%6.2f'%(wilson_intensity_scale_factor, wilson_b)

  #write out the scaled intensity
  F_scaled = miller_array_b.data()/(wilson_k * flex.exp(-2 * wilson_b * miller_array_b.sin_theta_over_lambda_sq().data()))
  sigF_scaled = miller_array_b.sigmas()/(wilson_k * flex.exp(-2 * wilson_b * miller_array_b.sin_theta_over_lambda_sq().data()))
  miller_array_b_scaled = miller_array_b.customized_copy(\
                data=F_scaled, sigmas=sigF_scaled)
  """

  miller_array_b_scaled, scale_factor = miller_array_a.scale(miller_array_b, resolution_dependent=True)

  mtz_dataset = miller_array_b_scaled.as_mtz_dataset(column_root_label="IOBS")
  hklb_farr = hklb.split('/')
  hklb_fname = hklb_farr[len(hklb_farr)-1]
  hklb_fname_arr = hklb_fname.split('.')
  hklb_fname_only = hklb_fname_arr[0]
  mtz_dataset.mtz_object().write(file_name=hklb_fname_only+'_scaled.mtz')

  #write as cns file
  f_cns = open(hklb_fname_only+'_scaled.hkl', 'w')
  miller_array_b_scaled.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()
  """
  common_a, common_b = miller_array_a.common_sets(miller_array_a, miller_array_b)
  plt.plot(mean_stol_sq, flex.log(expected_f_sq), linestyle='-', linewidth=2.0, c='k')
  plt.plot(mean_stol_sq, flex.log(mean_fobs_sq), linestyle='-', linewidth=2.0, c='r')
  plt.plot(mean_stol_sq, flex.log(mean_fobs_sq/(wilson_intensity_scale_factor * flex.exp(-2 * wilson_b * mean_stol_sq))), linestyle='-', linewidth=2.0, c='b')
  plt.xlabel('stol_sq')
  plt.ylabel('<I>')

  plt.show()
  """
  """
  offset = 0
  flex_fmodel = miller_array_a.amplitudes().data()
  flex_phifmodel = miller_array_a.phases().data()

  miller_set=miller.set(
            crystal_symmetry=miller_array_a.crystal_symmetry(),
            indices=miller_array_a.indices(),
            anomalous_flag=False)
  miller_array_out = miller_set.array(
            data=flex_fmodel).set_observation_type_xray_amplitude()

  mtz_dataset = miller_array_out.as_mtz_dataset(column_root_label="FP")

  for data,lbl,typ in [(flex_phifmodel, "PHIFMODEL", "P")]:
    mtz_dataset.add_miller_array(miller_array_out.array(data=data),
            column_root_label=lbl,
            column_types=typ)

  mtz_dataset.mtz_object().write(file_name='3u3e_fmodel_offset_01.mtz')



  #get matching amplitudes
  matches_template = miller.match_multi_indices(
                  miller_indices_unique=miller_array_a.indices(),
                  miller_indices=miller_array_b.indices())

  I_a = flex.double([miller_array_a.data()[pair[0]] for pair in matches_template.pairs()])
  I_b = flex.double([miller_array_b.data()[pair[1]] for pair in matches_template.pairs()])
  print 'Found ', len(I_a), ' reflections matched.'

  #fit b to a
  fit = flex.linear_regression(I_a, I_b)
  assert fit.is_well_defined()
  fit_y_intercept = fit.y_intercept()
  fit_slope = fit.slope()

  I_b_scaled = I_b - fit_y_intercept
  """
