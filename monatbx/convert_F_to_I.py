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
  hklin, hklout= read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()

  miller_array = miller_arrays[0]

  F_as_I = flex.sqrt(miller_array.data())
  sigF_as_sigI = flex.sqrt(miller_array.sigmas())
  miller_array_I = miller_array.customized_copy(data=F_as_I, sigmas=sigF_as_sigI)
  miller_array_I = miller_array_I.set_observation_type_xray_amplitude()

  for miller_index, d, F, sigF, I, sigI in zip(miller_array.indices(), miller_array.d_spacings().data(), miller_array.data(), miller_array.sigmas(), miller_array_I.data(), miller_array_I.sigmas()):
    print miller_index, d, F, sigF, I, sigI

  #write as mtz file
  mtz_dataset_out = miller_array_I.as_mtz_dataset(column_root_label="FOBS")
  mtz_dataset_out.mtz_object().write(file_name=hklout)
