'''
Fit second amplitude set to the first set.
Report CC, K, and B
'''

from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller
import sys
import matplotlib.pyplot as plt

def read_input(args):
  hklin = ''
  hklout = ''
  b_factor = 0
  offset = 0
  n_bins = 20
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklin':
      hklin = pair[1]
    elif pair[0]=='hklout':
      hklout = pair[1]
    elif pair[0]=='b_factor':
      b_factor = float(pair[1])
    elif pair[0]=='n_bins':
      n_bins = int(pair[1])
    elif pair[0]=='offset':
      offset = float(pair[1])

  if hklin == '' or hklout == '':
    print "Please provide input and output hkl file."
    exit()

  return hklin, hklout, b_factor, n_bins, offset


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hklin, hklout, b_factor, n_bins, offset = read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array_a = miller_arrays[0]
  print 'First reflection file:', hklin
  miller_array_a.show_summary()
  print 'B-factor', b_factor
  print 'Offset', offset
  I_obs = miller_array_a.data()

  #apply offset
  I_obs_scaled = I_obs - offset

  #apply b-factor
  sin_theta_over_lambda_sq = miller_array_a.sin_theta_over_lambda_sq().data()
  I_obs_scaled = I_obs *  flex.exp(-2 * b_factor *sin_theta_over_lambda_sq)
  miller_array_out = miller_array_a.customized_copy(data=I_obs_scaled).set_observation_type_xray_intensity()
  mtz_dataset = miller_array_out.as_mtz_dataset(column_root_label="IOBS")
  mtz_dataset.mtz_object().write(file_name=hklout)

  #write as cns file
  hklout_arr = hklout.split('.')
  f_cns = open(hklout_arr[0]+'.hkl', 'w')
  miller_array_out.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()

  #plot the before and after
  binner_template_asu = miller_array_a.setup_binner(n_bins=n_bins)
  binner_template_asu_indices = binner_template_asu.bin_indices()
  for i in range(1, n_bins+1):
    i_binner = (binner_template_asu_indices == i)
    I_obs_bin = I_obs.select(i_binner)
    I_obs_scaled_bin = miller_array_out.data().select(i_binner)

    print binner_template_asu.bin_d_range(i)[1], flex.mean(I_obs_bin), flex.mean(I_obs_scaled_bin)
