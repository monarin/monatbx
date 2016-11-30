from __future__ import division
from cctbx.array_family import flex
import numpy as np
import math
from cctbx import miller
from cctbx import crystal
import cPickle as pickle
from iotbx import reflection_file_reader

class polarity_handler(object):
  '''
  A wrapper class for organizer class
  '''

  def __init__(self, iparams):
    """
    Intialitze parameters
    """
    self.iparams = iparams

  def get_observations_non_polar(self, observations_original, observations_filename):

    """
    Determine polarity based on input data.
    The function still needs isomorphous reference so, if flag_polar is True,
    miller_array_iso must be supplied in input file.
    """
    observations_asu = observations_original.map_to_asu()

    if self.iparams.indexing_ambiguity.flag_on == False:
      return observations_asu

    pickle_filename_arr = observations_filename.split('/')
    if len(pickle_filename_arr) == 1:
      pickle_filename_only = pickle_filename_arr[0]
    else:
      pickle_filename_only = pickle_filename_arr[len(pickle_filename_arr)-1]

    cc_asu = 0
    cc_rev = 0
    if self.iparams.indexing_ambiguity.index_basis_in.endswith('mtz'):
      #use reference mtz file to determine polarity
      from mod_util import utility_handler
      util_hdlr = utility_handler()
      flag_mtz_found, miller_array_polar = util_hdlr.get_miller_array_from_mtz(self.iparams.indexing_ambiguity.index_basis_in)
      if self.iparams.target_anomalous_flag:
        miller_array_polar = miller_array_polar.generate_bijvoet_mates()
      miller_array_polar = miller_array_polar.resolution_filter(d_min=self.iparams.indexing_ambiguity.d_min, d_max=self.iparams.indexing_ambiguity.d_max)

      observations_asu = observations_original.map_to_asu()
      observations_rev = self.get_observations_non_polar(observations_original, self.iparams.indexing_ambiguity.assigned_basis)

      matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_polar.indices(),
                  miller_indices=observations_asu.indices())
      I_ref_match = flex.double([miller_array_polar.data()[pair[0]] for pair in matches.pairs()])
      I_obs_match = flex.double([observations_asu.data()[pair[1]] for pair in matches.pairs()])
      cc_asu = np.corrcoef(I_ref_match, I_obs_match)[0,1]
      n_refl_asu = len(matches.pairs())

      matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_polar.indices(),
                  miller_indices=observations_rev.indices())
      I_ref_match = flex.double([miller_array_polar.data()[pair[0]] for pair in matches.pairs()])
      I_obs_match = flex.double([observations_rev.data()[pair[1]] for pair in matches.pairs()])
      cc_rev = np.corrcoef(I_ref_match, I_obs_match)[0,1]
      n_refl_rev = len(matches.pairs())

      polar_hkl = 'h,k,l'
      if cc_rev > cc_asu:
        polar_hkl = self.iparams.indexing_ambiguity.assigned_basis

    else:
      #use basis in the given input file
      polar_hkl = 'h,k,l'
      basis_pickle = pickle.load(open(self.iparams.indexing_ambiguity.index_basis_in,"rb"))
      if observations_filename in basis_pickle:
        polar_hkl = basis_pickle[observations_filename]

    #return observations with correct polarity
    if polar_hkl == 'h,k,l':
      return observations_asu
    else:
      from cctbx import sgtbx
      cb_op = sgtbx.change_of_basis_op(polar_hkl)
      observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
      return observations_rev
