from iotbx import reflection_file_reader
from cctbx.array_family import flex
import sys

def read_input(args):
  hklin = ''
  hklout = 'reflection_out.mtz'
  d_upper_set = flex.double()
  d_lower_set = flex.double()
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='hklin':
      hklin = pair[1]
    elif pair[0]=='hklout':
      hklout = pair[1]
    elif pair[0]=='d_upper':
      d_upper_set.append(float(pair[1]))
    elif pair[0]=='d_lower':
      d_lower_set.append(float(pair[1]))

  if hklin == '':
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

  #start filtering out reflection within d_min and d_max
  # max res.
  # +
  # +
  # d_upper
  # -
  # -
  # d_lower
  # +
  # +
  # +
  # +
  # min res.

  for i in range(len(d_upper_set)):
    miller_indices = flex.miller_index()
    I_set = flex.double()
    sigI_set = flex.double()
    for miller_index, d, I, sigI in zip(miller_array.indices(), miller_array.d_spacings().data(),
                                      miller_array.data(), miller_array.sigmas()):
      if d > d_upper_set[i] or d < d_lower_set[i]:
        miller_indices.append(miller_index)
        I_set.append(I)
        sigI_set.append(sigI)
        print 'Accept', miller_index, d, I, sigI
      else:
        print 'Discard', miller_index, d, I, sigI

    n_before = len(miller_array.indices())
    miller_array = miller_array.customized_copy(indices = miller_indices,
          data = I_set,
          sigmas = sigI_set)

    print 'Length before: %6.0f, after: %6.0f'%(n_before, len(miller_array.indices()))



  #write as mtz file
  mtz_dataset_out = miller_array.as_mtz_dataset(column_root_label="IOBS")
  mtz_dataset_out.mtz_object().write(file_name=hklout+'.mtz')

  #write as cns file
  f_cns = open(hklout+'.hkl', 'w')
  miller_array.export_as_cns_hkl(file_object=f_cns)
  f_cns.close()
