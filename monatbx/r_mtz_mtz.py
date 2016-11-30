'''
Scale second mtz to the first (linear scale only) and calculate r-factors.
Note that all intensity array will be converted to amplitude.
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
  #read input parameters and frames (pickle files)
  hkla, hklb, d_min, d_max, n_bins = read_input(args = sys.argv[1:])
  #first reflection file (reference)
  reflection_file = reflection_file_reader.any_reflection_file(hkla)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array_a = miller_arrays[0]
  miller_array_a = miller_array_a.resolution_filter(d_min=d_min, d_max=d_max)
  if miller_array_a.is_xray_intensity_array():
    miller_array_a = miller_array_a.enforce_positive_amplitudes()
  print 'First reflection file:', hkla
  miller_array_a.show_summary()
  #second reflection file
  reflection_file = reflection_file_reader.any_reflection_file(hklb)
  miller_arrays=reflection_file.as_miller_arrays()
  miller_array_b = miller_arrays[0]
  miller_array_b = miller_array_b.resolution_filter(d_min=d_min, d_max=d_max)
  if miller_array_b.is_xray_intensity_array():
    miller_array_b = miller_array_b.enforce_positive_amplitudes()
  print
  print 'Second reflection file:', hklb
  miller_array_b.show_summary()
  ma_common_a, ma_common_b = miller_array_a.common_sets(miller_array_b)
  ma_common_a.setup_binner(n_bins=n_bins)
  ma_common_b.use_binning_of(ma_common_a)
  """
  print
  print "I/sigI first reflection file"
  miller_array_a.setup_binner(n_bins=n_bins)
  miller_array_a.i_over_sig_i(use_binning=True).show()
  print "Overall=",miller_array_a.i_over_sig_i(use_binning=False)
  print
  print "I/sigI second reflection file"
  miller_array_b.setup_binner(n_bins=n_bins)
  miller_array_b.i_over_sig_i(use_binning=True).show()
  print "Overall=",miller_array_b.i_over_sig_i(use_binning=False)
  """
  #scale b to a
  miller_array_b_scaled = miller_array_a.scale(miller_array_b, resolution_dependent=False)
  #get common sets
  ma_common_a, ma_common_b = miller_array_a.common_sets(miller_array_b_scaled)
  #set up binning
  ma_common_a.setup_binner(n_bins=n_bins)
  #calculate R-factor
  print
  print "R-factors"
  r1_factor_bin = ma_common_a.r1_factor(ma_common_b, use_binning=True)
  r1_factor_bin.show()
  r1_factor = ma_common_a.r1_factor(ma_common_b, use_binning=False)
  print 'Overall R-factor:', r1_factor
  #calculate cc
  print
  print "Correlations"
  ma_cc = ma_common_a.correlation(ma_common_b, use_binning=True)
  ma_cc.show()
  cc_scalar = ma_common_a.correlation(ma_common_b, use_binning=False)
  print 'Overall CC:', cc_scalar.coefficient()
