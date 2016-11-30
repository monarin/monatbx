from iotbx import reflection_file_reader
from cctbx.array_family import flex
import sys

def read_input(args):
  hklin = ''
  hklout = 'reflection_out.mtz'
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklin':
      hklin = pair[1]
    elif pair[0]=='hklout':
      hklout = pair[1]

  if hklin == '':
    print "Please provide input hkl file."
    exit()

  return hklin, hklout


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hklin, hklout = read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()

  miller_array = miller_arrays[0]
  miller_array.show_summary()

  i_sel = (miller_array.data() >= 0)
  ma_positive = miller_array.select(i_sel)
  print "N_refl all=%8.0f positive=%8.0f negative=%8.0f"%(len(miller_array.indices()), \
      len(ma_positive.indices()), len(miller_array.indices())-len(ma_positive.indices()))

  #write as mtz file
  mtz_dataset_out = ma_positive.as_mtz_dataset(column_root_label="IOBS")
  mtz_dataset_out.mtz_object().write(file_name=hklout+'.mtz')

  #write as cns file
  f_cns = open(hklout+'.hkl', 'w')
  miller_array.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()
