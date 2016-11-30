from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller
import sys

def read_input(args):
  hklin = ''
  hklout = 'reflection_out'
  d_limit = 0
  n_bins = 20
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklin':
      hklin = pair[1]
    elif pair[0]=='hklout':
      hklout = pair[1]
    elif pair[0]=='d_limit':
      d_limit = float(pair[1])
    elif pair[0]=='n_bins':
      n_bins = int(pair[1])

  if hklin == '':
    print "Please provide input hkl file."
    exit()

  return hklin, hklout, d_limit, n_bins


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hklin, hklout, d_limit, n_bins = read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()

  miller_array = miller_arrays[0]
  miller_array.show_summary()
  d_min = flex.min(miller_array.d_spacings().data())
  d_max = flex.max(miller_array.d_spacings().data())
  d_spacings = miller_array.d_spacings().data()

  #randomize amplitudes for reflections < d_min
  #random range is between 0 and max(intensity)
  #filter resolution
  i_sel_res = miller_array.resolution_filter_selection(d_min=0, d_max=d_limit)

  import random
  import numpy as np
  import math
  I_shuffle = miller_array.data()[:]
  sigI_shuffle = miller_array.sigmas()[:]
  I_sel = I_shuffle.select(i_sel_res)
  lower_I = np.median(I_sel) - (2.5*math.sqrt(np.median(I_sel)))
  upper_I = np.median(I_sel) + (2.5*math.sqrt(np.median(I_sel)))

  i_seq = flex.int(range(len(miller_array.data())))
  i_seq_sel = i_seq.select(i_sel_res)

  print 'Select ', len(i_seq_sel), ' reflections'
  for ii in range(len(i_seq_sel)):
    i = i_seq_sel[ii]
    tmp_I = random.uniform(lower_I, upper_I)
    print ii, i, miller_array.indices()[i], d_spacings[i],  I_shuffle[i], tmp_I
    I_shuffle[i] = tmp_I
    sigI_shuffle[i] = math.sqrt(tmp_I)*4

  print np.median(I_sel), np.std(I_sel), lower_I, upper_I
  miller_array_out = miller_array.customized_copy(data=I_shuffle,
        sigmas=sigI_shuffle)

  #write as mtz file
  mtz_dataset_out = miller_array_out.as_mtz_dataset(column_root_label="IOBS")
  mtz_dataset_out.mtz_object().write(file_name=hklout+'.mtz')

  #write as cns file
  f_cns = open(hklout+'.hkl', 'w')
  miller_array.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()

  #report intensity
  miller_array_template_asu = miller_array_out.complete_set().resolution_filter(
      d_min=d_min, d_max=d_max)
  binner_template_asu = miller_array_template_asu.setup_binner(n_bins=n_bins)
  binner_template_asu_indices = binner_template_asu.bin_indices()

  txt_out = 'Summary for '+hklout+'.mtz\n'
  txt_out += 'Bin Resolution Range       Completeness         <I>     <I/sigI>\n'
  txt_out += '------------------------------------------------------------------\n'
  for i in range(1,n_bins+1):
    i_binner = (binner_template_asu_indices == i)
    miller_indices_bin = miller_array_template_asu.indices().select(i_binner)

    matches_template = miller.match_multi_indices(
                  miller_indices_unique=miller_indices_bin,
                  miller_indices=miller_array_out.indices())

    I_bin = flex.double([miller_array_out.data()[pair[1]] for pair in matches_template.pairs()])
    sigI_bin = flex.double([miller_array_out.sigmas()[pair[1]] for pair in matches_template.pairs()])

    completeness = len(I_bin)/len(miller_indices_bin)

    txt_out += '%02d %7.2f - %7.2f %5.1f %6.0f / %6.0f %10.2f %10.2f\n' \
          %(i, binner_template_asu.bin_d_range(i)[0], binner_template_asu.bin_d_range(i)[1], completeness*100, \
          len(I_bin), len(miller_indices_bin), flex.mean(I_bin), flex.mean(I_bin/sigI_bin))

  print txt_out
