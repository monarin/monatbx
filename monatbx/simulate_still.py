from iotbx import reflection_file_reader
from cctbx.array_family import flex
import sys

def read_input(args):
  hklref = ''
  data = ''
  flag_crystal_orientation = False
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklref':
      hklref = pair[1]
    elif pair[0]=='data':
      data = pair[1]
    elif pair[0]=='flag_crystal_orientation':
      flag_crystal_orientation = bool(pair[1])

  if hklref == '':
    print "Please provide input hkl file."
    exit()

  if len(d_upper_set) != len(d_lower_set):
    print "Please specify equal number of d_min and d_max parameters."
    exit()

  return hklin, hklout, d_upper_set, d_lower_set


if (__name__ == "__main__"):

  #0 .read input parameters and frames (pickle files)
  hklin, hklout, d_upper_set, d_lower_set = read_input(args = sys.argv[1:])

  reflection_file = reflection_file_reader.any_reflection_file(hklin)
  miller_arrays=reflection_file.as_miller_arrays()

  miller_array = miller_arrays[0]
  miller_array.show_summary()
