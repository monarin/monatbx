from __future__ import division
from cctbx.array_family import flex
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from cctbx import miller
import cPickle as pickle
from sets import Set
from mod_observations import observations_handler

class organize_handler(object):
  """
  A wrapper class for organizer class
  """

  def __init__(self, iparams):
    """
    Intialitze parameters
    """
    self.iparams = iparams

  def get_common_images(self, obs_hdlr_main, obs_hdlr_no, obs_hdlr_set):
    """
    Determine no. of common images and return a list of those images (indices).
    """
    conn_images_list = []
    n_conn = 0
    sum_refl = 0
    for obs_hdlr_other, obs_hdlr_no_other in zip(obs_hdlr_set, range(len(obs_hdlr_set))):
      if obs_hdlr_no != obs_hdlr_no_other:
        try:
          ma_x, ma_y = obs_hdlr_main.observations.common_sets(obs_hdlr_other.observations)
          if len(ma_x.indices()) > self.iparams.cluster.n_min_refl:
            n_conn += 1
            conn_images_list.append(obs_hdlr_no_other)
            sum_refl += len(ma_x.indices())
        except Exception:
          pass

    avg_common_refl = 0
    if n_conn > 0:
      avg_common_refl = sum_refl/n_conn

    print ' {0:40} ==> N_CONN:{1:5d} N_REFL:{2:5d} AVG_COMMON_REFL:{3:6.2f}'.format(obs_hdlr_main.pickle_filename_only, \
          n_conn, len(obs_hdlr_main.observations.indices()), avg_common_refl)
    return conn_images_list, n_conn

  def generate_clusters(self, obs_hdlr_set, conn_images_set, n_conn_list):
    """
    Get n_cluster from images based on their common reflections.
    Each file is assigned with a list of n_common_images and cc_images with
    length = n_images
    """
    #select frame with maximum no. of connected frames as the first cluster
    #continue down the ranking list by discarding duplicates.
    i_sort = np.argsort(np.array(n_conn_list))
    i_sort = i_sort[::-1]
    frames_cluster_set = []
    for i in i_sort:
      if len(frames_cluster_set) < self.iparams.cluster.n_max_clusters:
        #check if frame i is already a part of the frame clusters.
        flag_frame_ok = True
        for frames in frames_cluster_set:
          if len(np.extract(frames == i, frames)) > 0:
            flag_frame_ok = False
            break

        if flag_frame_ok:
          #discard frames that are already part of previous clusters.
          conn_frames_unused = []
          for conn_frame in np.array(conn_images_set[i]):
            flag_subframe_ok = True
            for frames in frames_cluster_set:
              if len(np.extract(frames == conn_frame, frames)) > 0:
                flag_subframe_ok = False
                break
            if flag_subframe_ok:
              conn_frames_unused.append(conn_frame)

          conn_frames = np.array(conn_frames_unused)
          conn_frames_of_this = np.array([n_conn_list[i_conn_frame] for i_conn_frame in conn_frames])
          i_sort_sub = np.argsort(conn_frames_of_this)
          i_sort_sub = i_sort_sub[::-1]
          conn_frames_sort = np.array([conn_frames[i_i_sort_sub] for i_i_sort_sub in i_sort_sub])

          conn_frames_sel_main = np.array([i])
          if len(conn_frames_sort) < self.iparams.cluster.n_max_images_per_cluster:
            conn_frames_sel_other = conn_frames_sort[:]
          else:
            conn_frames_sel_other = conn_frames_sort[:self.iparams.cluster.n_max_images_per_cluster]
          conn_frames_sel = np.concatenate((conn_frames_sel_main, conn_frames_sel_other))
          if len(conn_frames_sel) > 1:
            frames_cluster_set.append(conn_frames_sel)

    #determine frames that are not connected.
    frames_selected_set = Set()
    frames_all_set = Set(range(len(obs_hdlr_set)))
    for frames_cluster in frames_cluster_set: frames_selected_set = frames_selected_set.union(Set(frames_cluster))
    frames_left_set = frames_all_set.difference(frames_selected_set)

    #output list of image paths
    sum_frames_selected = 0
    txt_out = '\nSummary of image clustering.\n'
    cn_i = 0
    for frames_cluster in frames_cluster_set:
      txt_out += ' cluster (%3.0f) --> N_FRAMES: %6.0f\n'%(cn_i+1, len(frames_cluster))
      sum_frames_selected += len(frames_cluster)
      txt_frame = ''
      for i_frame in frames_cluster:
        txt_frame += obs_hdlr_set[i_frame].pickle_filename + '\n'
      f = open(self.iparams.run_no+'/clusters/cluster_'+str(cn_i+1)+'.lst', 'w')
      f.write(txt_frame)
      f.close()
      cn_i += 1

    #output unconnected frames
    txt_frame = ''
    for i_frame in frames_left_set:
      txt_frame += obs_hdlr_set[i_frame].pickle_filename + '\n'
    f = open(self.iparams.run_no+'/clusters/cluster.unc', 'w')
    f.write(txt_frame)
    f.close()

    txt_out += '%6.0f images in %3.0f clusters (%6.0f images not classified)\n'%(sum_frames_selected, len(frames_cluster_set), len(frames_left_set))
    return txt_out
