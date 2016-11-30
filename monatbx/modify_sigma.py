from iotbx import reflection_file_reader
import sys
import matplotlib.pyplot as plt
import numpy as np
from cctbx.array_family import flex

hklin = sys.argv[1]
flag_plot = False
reflection_file = reflection_file_reader.any_reflection_file(hklin)
miller_arrays = reflection_file.as_miller_arrays()
for miller_array in miller_arrays:
  if miller_array.sigmas() is not None and miller_array.is_xray_intensity_array():
    miller_array.show_summary()
    #modify sigma by bin
    binner = miller_array.setup_binner(auto_binning=True)
    binner_indices = binner.bin_indices()
    w_sigp_coeff = [1.5, 1.0]
    #prepare output
    out_indices = flex.miller_index()
    out_data = flex.double()
    out_sigmas = flex.double()
    for i in range(1, binner.n_bins_used()+1):
      i_binner = (binner_indices == i)
      ma_bin = miller_array.select(i_binner)
      m = (w_sigp_coeff[1]-w_sigp_coeff[0])/(flex.min(ma_bin.data())-flex.max(ma_bin.data()))
      b = w_sigp_coeff[1] - (m * flex.min(ma_bin.data()))
      w_sigp = (m * ma_bin.data()) + b
      if flag_plot:
        sel_perm = flex.sort_permutation(ma_bin.data(), reverse=True)
        ma_bin_sorted = ma_bin.select(sel_perm)
        w_sigp_sorted = w_sigp.select(sel_perm)
        plt.subplot(3,1,1)
        plt.plot(ma_bin_sorted.data(), color='navy')
        plt.subplot(3,1,2)
        plt.plot(w_sigp_sorted, color='firebrick')
        plt.subplot(3,1,3)
        plt.plot(ma_bin_sorted.sigmas(), color='darksage')
        plt.plot(ma_bin_sorted.sigmas() * w_sigp_sorted, color='yellowgreen')
        plt.show()
      out_indices.extend(ma_bin.indices())
      out_data.extend(ma_bin.data())
      out_sigmas.extend(ma_bin.sigmas() * w_sigp)
    #write out mtz
    miller_array_out = miller_array.customized_copy(indices=out_indices, \
      data=out_data, sigmas=out_sigmas)
    mtz_dataset = miller_array_out.as_mtz_dataset(column_root_label="IOBS")
    hklin_farr = hklin.split('/')
    hklin_fname = hklin_farr[len(hklin_farr)-1]
    hklin_fname_arr = hklin_fname.split('.')
    hklin_fname_only = hklin_fname_arr[0]
    mtz_dataset.mtz_object().write(file_name=hklin_fname_only+'_sigmap_02.mtz')
