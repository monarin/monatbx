from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller
import sys
import numpy as np
import random
import math

def read_input(args):
  hklin = ''
  hklout = 'reflection_out'
  b_factor = 0
  wavelength = 0
  n_bins = 20
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklin':
      hklin = pair[1]
    elif pair[0]=='hklout':
      hklout = pair[1]
    elif pair[0]=='b_factor':
      b_factor = float(pair[1])
    elif pair[0]=='wavelength':
      wavelength = float(pair[1])
    elif pair[0]=='n_bins':
      n_bins = int(pair[1])

  if hklin == '':
    print "Please provide input hkl file."
    exit()

  return hklin, hklout, b_factor, wavelength, n_bins


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hklin, hklout, b_factor, wavelength, n_bins = read_input(args = sys.argv[1:])
  print 'Input parameters:'
  print 'hklin=', hklin
  print 'hklout=', hklout
  print 'b_factor=', b_factor
  print 'n_bins=', n_bins

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()

  miller_array = miller_arrays[1]
  miller_array.show_summary()
  d_min = flex.min(miller_array.d_spacings().data())
  d_max = flex.max(miller_array.d_spacings().data())
  d_spacings = miller_array.d_spacings().data()

  #get an exponential increase of missing reflections as a function of resolution
  binner = miller_array.setup_binner(n_bins=n_bins)
  binner_indices = binner.bin_indices()
  stol_sq = miller_array.sin_theta_over_lambda_sq()
  stol_sq.use_binner_of(miller_array)
  mean_stol_sq = flex.double(stol_sq.mean(
        use_binning=True,
        use_multiplicities=True).data[1:-1])
  percent_keep = np.exp(-2*b_factor*mean_stol_sq)

  miller_indices_sel = flex.miller_index()
  I_sel = flex.double()
  sigI_sel = flex.double()
  moment_I = flex.double()
  moment_I_sel = flex.double()
  for i in range(binner.n_bins_used()):
    i_binner = (binner_indices == i)
    miller_array_bin = miller_array.select(i_binner)
    if len(miller_array_bin.data()) > 0:
      n_refl_bin = len(miller_array_bin.data())
      n_kept_refl = int(round(n_refl_bin * percent_keep[i]))


      import operator
      I_bin_sort_index= [p for (p,q) in sorted(enumerate(miller_array_bin.data()), key=operator.itemgetter(1))]
      #I_bin_sort_index.reverse()    #uncomment for missing weak

      i_sel = I_bin_sort_index[:n_kept_refl]
      #for k in range(10):
      #  print miller_array_bin.data()[I_bin_sort_index[k]]


      #i_sel = random.sample(range(n_refl_bin), n_kept_refl)

      i_sel_bool = flex.bool([False] * n_refl_bin)
      for i_i_sel in i_sel:
        i_sel_bool[i_i_sel] = True

      miller_indices_sel.extend(miller_array_bin.indices().select(i_sel_bool))
      I_sel.extend(miller_array_bin.data().select(i_sel_bool))
      sigI_sel.extend(miller_array_bin.sigmas().select(i_sel_bool))

      moment_I_bin = np.mean(miller_array_bin.data()**2)/(np.mean(miller_array_bin.data())**2)
      moment_I_sel_bin = np.mean(miller_array_bin.data().select(i_sel_bool)**2)/(np.mean(miller_array_bin.data().select(i_sel_bool))**2)

      moment_I.append(moment_I_bin)
      moment_I_sel.append(moment_I_sel_bin)

      print '%6.2f - %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'%(binner.bin_d_range(i)[0], binner.bin_d_range(i)[1], 100, moment_I_bin, percent_keep[i]*100, moment_I_sel_bin, np.mean(miller_array_bin.data().select(i_sel_bool)))

  print '     TOTAL      %6.2f %6.2f %6.2f %6.2f'%(100, np.mean(moment_I), (len(I_sel)/len(miller_array.data()))*100, np.mean(moment_I_sel))




  """
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
  """
