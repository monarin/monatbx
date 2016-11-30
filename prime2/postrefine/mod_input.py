from __future__ import division
import iotbx.phil
from libtbx.utils import Usage, Sorry
import sys, os

master_phil = iotbx.phil.parse("""
data = None
  .type = path
  .multiple = True
  .help = Directory containing integrated data in pickle format.  Repeat to \
    specify additional directories.
run_no = None
  .type = str
  .help = Run no. is used as folder name that stores output files.
  .optional = False
title = None
  .type = str
  .help = Title of the run.
  .multiple = False
d_min = 0.1
  .type = float
  .help = Minimum resolution.
d_max = 99
  .type = float
  .help = Maximum resolution.
target_unit_cell = None
  .type = unit_cell
  .help = Target unit-cell parameters are used to discard outlier cells.
  .optional = False
target_space_group = None
  .type = str
  .help = Target space group.
  .optional = False
target_anomalous_flag = False
  .type = bool
  .help = Target anomalous flag (False = Not anomalous data)
  .optional = False
pixel_size_mm = None
  .type = float
  .help = Pixel size in mm. (MAR = 0.079346)
  .optional = False
gamma_e = 0.003
  .type = float
  .help = Initial spread of the energy spectrum (1/Angstrom).
flag_beam_divergence = False
  .type = bool
  .help = Default is not to refine beam divergence. Set to True to allow gammaX and gammaY refinement.
n_processors = 32
  .type = int
  .help = No. of processing units
  .optional = True
n_bins = 20
  .type = int
  .help = No. of bins
indexing_ambiguity
  .help = "Parameters used in resolving indexing ambiguity"
{
  flag_on = False
    .type = bool
    .help = Set to True to allow the program to read in polarity info. \
      from the pickle file specified in index_basis_in.
  index_basis_in = None
    .type = path
    .help = Pickle file storing polarity info or an mtz file. (output from Brehm & \
      Diederichs clustering algorithm).
  assigned_basis = h,k,l
    .type = str
    .help = In case index_basis_in given is an mtz file, you can specify a basis that each integration can be converted to.
  d_min = 10.0
    .type = float
    .help = In case index_basis_in given is an mtz file, you can pecify minimum resolution used to calculate correlation with the given mtz file.
  d_max = 20.0
    .type = float
    .help = In case index_basis_in given is an mtz file, you can pecify maximum resolution used to calculate correlation with the given mtz file.
}
cluster
  .help = Parameters used in image clustering.
{
  n_min_refl = 50
    .type = int
    .help = Minimum no. of common reflections for the determination of a connection.
  n_max_images_per_cluster = 10
    .type = int
    .help = Maximum no. of images per a cluster.
  n_max_clusters = 20
    .type = int
    .help = Maximum no. of clusters for image clustering.
}
hklisoin = None
  .type = path
  .help = Mtz file for the calculation of CCiso
""")

txt_help = """**************************************************************************************************

Prime: post-refinement and merging.

For more detail and citation, see Enabling X-ray free electron laser crystallography
for challenging biological systems from a limited number of crystals
"DOI: http://dx.doi.org/10.7554/eLife.05421".

Usage: prime.postrefine parameter.phil

With this command, you can specify all parameters required by prime in your parameter.phil file.
To obtain the template of these parameters, you can perform a dry run (simply run prime.postrefine).
You can then change the values of the parameters.

For feedback, please contact monarin@stanford.edu.

**************************************************************************************************

List of available parameters:
"""

def process_input(argv=None):
  if argv == None:
    argv = sys.argv[1:]

  user_phil = []
  for arg in sys.argv[1:]:
    if os.path.isfile(arg):

      user_phil.append(iotbx.phil.parse(open(arg).read()))
    elif (os.path.isdir(arg)) :
      user_phil.append(iotbx.phil.parse("""data=\"%s\"""" % arg))
    else :
      print arg
      if arg == '--help' or arg == '-h':
        print txt_help
        master_phil.show(attributes_level=1)
        exit()
      try :
        user_phil.append(iotbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  if (len(params.data) == 0):
    master_phil.show()
    raise Usage("Use the above list of parameters to generate your input file (.phil). For more information, run prime.postrefine -h.")

  #generate run_no folder
  if os.path.exists(params.run_no):
    raise Sorry("The run number %s already exists."%params.run_no)

  os.makedirs(params.run_no)
  os.makedirs(params.run_no+'/clusters')

  #capture input read out by phil
  from cStringIO import StringIO
  class Capturing(list):
    def __enter__(self):
      self._stdout = sys.stdout
      sys.stdout = self._stringio = StringIO()
      return self
    def __exit__(self, *args):
      self.extend(self._stringio.getvalue().splitlines())
      sys.stdout = self._stdout

  with Capturing() as output:
    working_phil.show()

  txt_out = 'prime.postrefine input:\n'
  for one_output in output:
    txt_out += one_output + '\n'

  return params, txt_out
