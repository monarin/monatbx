from iotbx import reflection_file_reader
from cctbx.array_family import flex
import sys

def read_input(args):
  hklin = ''
  hklout = 'reflection_out.mtz'
  fraction_percent = 0.05
  n_cones = 1
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklin':
      hklin = pair[1]
    elif pair[0]=='hklout':
      hklout = pair[1]
    elif pair[0]=='fraction_percent':
      fraction_percent = float(pair[1])
    elif pair[0]=='n_cones':
      n_cones = int(pair[1])

  if hklin == '':
    print "Please provide input hkl file."
    exit()

  return hklin, hklout, fraction_percent, n_cones


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hklin, hklout, fraction_percent, n_cones = read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()

  miller_array = miller_arrays[0]
  miller_array.show_summary()

  ma_updated = miller_array.remove_cone(fraction_percent, axis_point_2=(1,1,1), negate=False)

  #write as mtz file
  mtz_dataset_out = ma_updated.as_mtz_dataset(column_root_label="IOBS")
  mtz_dataset_out.mtz_object().write(file_name=hklout)

  #write as cns file
  f_cns = open(hklout+'.hkl', 'w')
  ma_updated.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()
