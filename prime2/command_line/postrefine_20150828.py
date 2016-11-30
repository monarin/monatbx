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
from prime2.postrefine.mod_organize import organize_handler
from prime2.postrefine.mod_postrefine import postrefine_handler
from prime2.postrefine.mod_observations import observations_handler
from prime2.postrefine.mod_merge import merge_handler
import matplotlib.pyplot as plt
from cctbx import miller

def get_observations_handler_mproc(frame_no, frame_files, iparams):
  obs_hdlr = observations_handler(iparams)
  flag_ok, txt_obs = obs_hdlr.init_params_from_file(frame_files[frame_no])
  if flag_ok == False:
    return None
  else:
    return obs_hdlr

def get_common_images_mproc(obs_hdlr_no, observations_handler_set, iparams):
  orgh = organize_handler(iparams)
  result = orgh.get_common_images(observations_handler_set[obs_hdlr_no], obs_hdlr_no, observations_handler_set)
  return result

def scale_cluster_mproc(cluster_no, cluster_paths, iparams):
  frame_files = read_pickles([cluster_paths[cluster_no]])

  #scale and merge
  observations_hdlr_set = []
  txt_scale = ' SCALE CLUSTER %3.0f\n'%(cluster_no+1)
  postrefine_hdlr = postrefine_handler(iparams)
  for i_frame in range(len(frame_files)):
    observations_hdlr_0, txt_out = postrefine_hdlr.scale_0(frame_files[i_frame])
    txt_scale += '\n' + txt_out
    if observations_hdlr_0 is not None:
      observations_hdlr_0.refresh()
      observations_hdlr_set.append(observations_hdlr_0)

  merge_hdlr = merge_handler(iparams)
  observations_hdlr_ref, observations_hdlr_selected_set, txt_out = merge_hdlr.merge_multi_observations_handler(observations_hdlr_set, flag_show_summary=False)
  txt_scale += '\n' + txt_out
  if observations_hdlr_ref is None:
    print txt_scale
    return None, None, txt_scale + '\n'

  print txt_scale
  return observations_hdlr_ref, observations_hdlr_selected_set, txt_scale + '\n'

def postrefine_mproc(obs_hdlr_no, observations_hdlr_ref, observations_hdlr_set, iparams):
  observations_hdlr = observations_hdlr_set[obs_hdlr_no]
  postrefine_hdlr = postrefine_handler(iparams)
  observations_hdlr_out, txt_out = postrefine_hdlr.postrefine_with_reference(observations_hdlr_ref, observations_hdlr)
  observations_hdlr_out.refresh()
  print txt_out
  return observations_hdlr_out, txt_out + '\n'

def scale(observations_hdlr_ref, observations_hdlr, observations_hdlr_set, iparams):
  #scale two observations and return the merged observations and the
  #updated observations set.
  txt_scale = ''
  postrefine_hdlr = postrefine_handler(iparams)
  observations_hdlr, txt_out = postrefine_hdlr.scale_with_reference(observations_hdlr_ref, observations_hdlr=observations_hdlr, use_binning=True)
  txt_scale += txt_out

  if observations_hdlr is not None:
    merge_hdlr = merge_handler(iparams)
    observations_hdlr_merge, txt_out = merge_hdlr.merge_observations_handler(observations_hdlr_ref, observations_hdlr)
    txt_scale += '\n' + txt_out

    for observations_hdlr_in_set in observations_hdlr_set:
      observations_hdlr_in_set.set_params(postref_params=observations_hdlr.postref_params)
      observations_hdlr_in_set.refresh()
  else:
    observations_hdlr_merge = observations_hdlr_ref
    observations_hdlr_set = None

  print txt_scale
  return observations_hdlr_merge, observations_hdlr_set, txt_scale + '\n'

def scale_observations_mproc(observations_no, observations_hdlr_ref, iparams, observations_hdlr_set=None, observations_file_set=None):
  #scale an observation set to the reference set
  txt_scale = ''
  postrefine_hdlr = postrefine_handler(iparams)
  if observations_hdlr_set is not None:
    observations_hdlr, txt_out = postrefine_hdlr.scale_with_reference(observations_hdlr_ref, observations_hdlr=observations_hdlr_set[observations_no], use_binning=True)
  elif observations_file_set is not None:
    observations_hdlr, txt_out = postrefine_hdlr.scale_with_reference(observations_hdlr_ref, filename=observations_file_set[observations_no], use_binning=True)

  txt_scale += txt_out
  print txt_scale
  return observations_hdlr, txt_scale + '\n'

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

  frame_files = read_pickles(iparams.data)
  frames = range(len(frame_files))
  txt_step0 = 'STEP 0: read and analyze input observations\n'
  print txt_step0

  def get_observations_handler_mproc_wrapper(arg):
    return get_observations_handler_mproc(arg, frame_files, iparams)

  results = pool_map(
            iterable=frames,
            func=get_observations_handler_mproc_wrapper,
            processes=iparams.n_processors)

  observations_handler_set = []
  for result in results:
    if result is not None:
      observations_handler_set.append(result)
  txt_out = ' found %6.0f files with %6.0f files pass the analysis test.\n'%(len(frames), len(observations_handler_set))
  print txt_out
  txt_step0 += txt_out

  #1. cluster input images
  txt_step1 = 'STEP 1: clustering input images\n'
  print txt_step1

  def get_common_images_mproc_wrapper(arg):
    return get_common_images_mproc(arg, observations_handler_set, iparams)

  results = pool_map(
            iterable=range(len(observations_handler_set)),
            func=get_common_images_mproc_wrapper,
            processes=iparams.n_processors)

  conn_images_set = []
  n_conn_list = []
  for result in results:
    if result is not None:
      conn_images_list, n_conn = result
      conn_images_set.append(conn_images_list)
      n_conn_list.append(n_conn)
    else:
      conn_images_set.append([])
      n_conn_list.append(0)

  orgh = organize_handler(iparams)
  txt_cluster = orgh.generate_clusters(observations_handler_set, conn_images_set, n_conn_list)
  print txt_cluster
  txt_step1 += txt_cluster

  #2. scale images in each cluster
  txt_step2 = 'STEP 2: scale images in each cluster\n'
  print txt_step2
  path_to_cluster = iparams.run_no+'/clusters'
  postref_hdlr = postrefine_handler(iparams)
  cluster_paths = []
  for cluster_lst_file in os.listdir(path_to_cluster):
    observations_hdlr_ref = None
    if cluster_lst_file.endswith('.lst'):
      cluster_paths.append(path_to_cluster+'/'+cluster_lst_file)

  def scale_cluster_mproc_wrapper(arg):
    return scale_cluster_mproc(arg, cluster_paths, iparams)

  results = pool_map(
            iterable=range(len(cluster_paths)),
            func=scale_cluster_mproc_wrapper,
            processes=iparams.n_processors)

  observations_hdlr_scaled_set = []
  i_found_good_cluster = 0
  for i_found_good_cluster in range(len(results)):
    observations_hdlr_ref, observations_hdlr_set, txt_out = results[0]
    txt_step2 += txt_out
    if observations_hdlr_ref is not None:
      for obs_hdlr in observations_hdlr_set:
        observations_hdlr_scaled_set.append(obs_hdlr)
      break

  for i_cluster in range(i_found_good_cluster+1, len(results)):
    if results[i_cluster] is not None:
      observations_hdlr, observations_hdlr_set, txt_out = results[i_cluster]
      txt_step2 += txt_out
      if observations_hdlr is not None:
        observations_hdlr_ref, observations_hdlr_set, txt_out = scale(observations_hdlr_ref, observations_hdlr, observations_hdlr_set, iparams)
        txt_step2 += txt_out
        if observations_hdlr_set is not None:
          for obs_hdlr in observations_hdlr_set:
            observations_hdlr_scaled_set.append(obs_hdlr)

  #do first merge
  merge_hdlr = merge_handler(iparams)
  observations_hdlr_scaled_set_merge, observations_hdlr_selected_set, txt_out = merge_hdlr.merge_multi_observations_handler(observations_hdlr_scaled_set, flag_show_summary=True)
  txt_step2 += txt_out
  """
  #use the first merge to remove outlier frames - do second merge
  def scale_observations_mproc_wrapper(arg):
    return scale_observations_mproc(arg, observations_hdlr_scaled_set_merge, iparams, observations_hdlr_set=observations_hdlr_selected_set)

  results = pool_map(
            iterable=range(len(observations_hdlr_selected_set)),
            func=scale_observations_mproc_wrapper,
            processes=iparams.n_processors)
  observations_hdlr_selected_set = []
  for result in results:
    if result is not None:
      obs_hdlr, txt_out = result
      if obs_hdlr is not None:
        observations_hdlr_selected_set.append(obs_hdlr)
  observations_hdlr_scaled_set_merge, observations_hdlr_selected_set, txt_out = merge_hdlr.merge_multi_observations_handler(observations_hdlr_selected_set, flag_show_summary=True)
  txt_step2 += txt_out

  #use the second merge to scale all other observations
  frame_files = read_pickles([iparams.run_no + '/clusters/cluster.unc'])
  def scale_observations_mproc_wrapper(arg):
    return scale_observations_mproc(arg, observations_hdlr_scaled_set_merge, iparams, observations_file_set=frame_files)

  results = pool_map(
            iterable=range(len(frame_files)),
            func=scale_observations_mproc_wrapper,
            processes=iparams.n_processors)
  for result in results:
    if result is not None:
      obs_hdlr, txt_out = result
      if obs_hdlr is not None:
        observations_hdlr_selected_set.append(obs_hdlr)

  observations_hdlr_scaled_set_merge, observations_hdlr_selected_set, txt_out = merge_hdlr.merge_multi_observations_handler(observations_hdlr_selected_set, flag_show_summary=True)
  txt_step2 += txt_out
  """
  #3. post-refine all images
  txt_step3 = 'STEP 3: post-refine all images\n'
  """
  print txt_step3
  def postrefine_mproc_wrapper(arg):
    return postrefine_mproc(arg, observations_hdlr_scaled_set_merge, observations_hdlr_scaled_set, iparams)

  results = pool_map(
            iterable=range(len(observations_hdlr_scaled_set)),
            func=postrefine_mproc_wrapper,
            processes=iparams.n_processors)
  observations_hdlr_postrefined_set = []
  for result in results:
    if result is not None:
      observations_hdlr_out, txt_out = result
      txt_step3 += txt_out
      observations_hdlr_postrefined_set.append(observations_hdlr_out)
  merge_hdlr = merge_handler(iparams)
  observations_hdlr_postrefined_set_merge, txt_out = merge_hdlr.merge_multi_observations_handler(observations_hdlr_postrefined_set, flag_show_summary=True)
  txt_step3 += txt_out + '\n'
  print txt_out
  """
  time_global_end=datetime.now()
  time_global_spent=time_global_end-time_global_start
  txt_out_time_spent = 'Total calculation time: '+'{0:.2f}'.format(time_global_spent.seconds)+' seconds\n'
  print txt_out_time_spent

  txt_main = txt_input + txt_step1 + txt_step2 + txt_step3 + txt_out_time_spent
  f = open(iparams.run_no+'/log.txt', 'w')
  f.write(txt_main)
  f.close()
