from iotbx import reflection_file_reader
import sys

def read_input(args):
  hklin = ''
  hklout = 'reflection_out.hkl'
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

  f_cns = open(hklout, 'w')
  miller_array.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()
