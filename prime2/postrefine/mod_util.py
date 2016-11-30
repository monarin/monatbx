from __future__ import division
from iotbx import reflection_file_reader

class utility_handler(object):
  """
  Author      : Uervirojnangkoorn, M.
  Created     : 8/15/2015
  A collection of utility functions
  """

  def __init__(self):
    """
    Constructor
    """

  def get_miller_array_from_mtz(self, mtz_filename):
    flag_hklisoin_found = False
    miller_array_iso = None
    if mtz_filename is not None:
      flag_hklisoin_found = True
      reflection_file_iso = reflection_file_reader.any_reflection_file(mtz_filename)
      miller_arrays_iso=reflection_file_iso.as_miller_arrays()
      is_found_iso_as_intensity_array = False
      is_found_iso_as_amplitude_array = False
      for miller_array in miller_arrays_iso:
        if miller_array.is_xray_intensity_array():
          miller_array_iso = miller_array.deep_copy()
          is_found_iso_as_intensity_array = True
          break
        elif miller_array.is_xray_amplitude_array():
          is_found_iso_as_amplitude_array = True
          miller_array_converted_to_intensity = miller_array.as_intensity_array()
      if is_found_iso_as_intensity_array == False:
        if is_found_iso_as_amplitude_array:
          miller_array_iso = miller_array_converted_to_intensity.deep_copy()
        else:
          flag_hklisoin_found = False

    return flag_hklisoin_found, miller_array_iso
