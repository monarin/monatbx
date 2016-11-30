from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime2.postrefine
'''
Author      : Uervirojnangkoorn, M.
Created     : 7/28/2015
Description : Commands linked to prime2 libraries.
'''
import os, sys
import numpy as np
import math
from datetime import datetime, time
from libtbx.easy_mp import pool_map
from cctbx.array_family import flex
from prime2.postrefine.mod_cluster import cluster_handler
from prime2.postrefine.mod_observations import observations_handler
from prime2.postrefine.mod_postrefine import postrefine_handler

def get_ma_match_dict(image_seq, image_dict, iparams):
  postref_hdlr = postrefine_handler(iparams)
  observations_hdlr, txt_out = postref_hdlr.scale_0(image_dict[image_seq], image_seq)
  observations_full_dict = observations_hdlr.get_full_observations_as_dict()
  return observations_full_dict

def read_pickles(data):
  frame_files = []
  for p in data:
    if os.path.isdir(p) == False:
      #check if list-of-pickle text file is given
      pickle_list_file = open(p,'r')
      pickle_list = pickle_list_file.read().split("\n")
      for pickle_filename in pickle_list:
        if os.path.isfile(pickle_filename):
          frame_files.append(pickle_filename)
    else:
      for pickle_filename in os.listdir(p):
        if pickle_filename.endswith('.pickle'):
          frame_files.append(p+'/'+pickle_filename)

  #check if pickle_dir is given in input file instead of from cmd arguments.
  if len(frame_files)==0:
    print 'No pickle files found.'
    exit()

  return frame_files


if (__name__ == "__main__"):
  #capture starting time
  time_global_start=datetime.now()
  import logging
  logging.captureWarnings(True)
  formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
  console_handler = logging.StreamHandler()
  console_handler.setLevel(logging.ERROR)
  console_handler.setFormatter(formatter)
  logging.getLogger().addHandler(console_handler)
  logging.getLogger('py.warnings').addHandler(console_handler)
  logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', level=logging.DEBUG)

  #0. read input parameters and frames (pickle files)
  from prime2.postrefine.mod_input import process_input
  iparams, txt_input = process_input(sys.argv[:1])
  print txt_input

  image_file_set = read_pickles(iparams.data)
  image_seq_set = range(len(image_file_set))
  image_dict = {}
  for image_file, image_seq in zip(image_file_set, image_seq_set):
    image_dict[image_seq] = image_file

  #1. cluster input images
  txt_step1 = 'STEP 1: clustering input images\n'
  print txt_step1

  # get a dictionary of image clusters using miller index as key
  cluster_hdlr = cluster_handler(iparams)
  ma_image_seq_dict = cluster_hdlr.get_dict()
  ma_data_dict = cluster_hdlr.get_dict()

  def get_ma_match_dict_wrapper(arg):
    return get_ma_match_dict(arg, image_dict, iparams)

  step1_results = pool_map(
            iterable=image_seq_set,
            func=get_ma_match_dict_wrapper,
            processes=iparams.n_processors)


  for miller_index_key in ma_image_seq_dict.keys():
    for result in step1_results:
      if result is not None:
        if miller_index_key in result:
          ma_image_seq_dict[miller_index_key].append(result[miller_index_key][0])
          ma_data_dict[miller_index_key].append(result[miller_index_key][1])
    print miller_index_key, np.mean(ma_data_dict[miller_index_key])/np.std(ma_data_dict[miller_index_key]), np.median(ma_data_dict[miller_index_key])/np.std(ma_data_dict[miller_index_key])

  txt_main = txt_input + txt_step1
  f = open(iparams.run_no+'/log.txt', 'w')
  f.write(txt_main)
  f.close()
