"""
Author      : Uervirojnangkoorn, M.
Created     : 8/2/2015
Description : Least-square refinement Ir-Io.

"""
from __future__ import division
from scipy import optimize
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation, basis_type
from cctbx import statistics
from mod_partiality import partiality_handler
from mod_observations import observations_handler

class leastsqr_handler(object):
  '''
  A wrapper class for least-squares refinement
  '''

  def __init__(self, iparams):
    """
    Intialitze parameters
    """
    self.iparams = iparams
    self.sigma_I_cutoff = 2.5
    self.scale_0_d_max = 5.0

  def prep_input(self, params, cs):
    """
    From crystal system cs, determine refined parameters
    """
    a, b, c, alpha, beta, gamma = params
    if cs == 'Triclinic':
      x0 = params
    elif cs == 'Monoclinic':
      x0 = np.array([a,b,c,beta])
    elif cs == 'Orthorhombic':
      x0 = np.array([a,b,c])
    elif cs == 'Tetragonal':
      x0 = np.array([a,c])
    elif cs == 'Trigonal' or cs == 'Hexagonal':
      x0 = np.array([a,c])
    elif cs == 'Cubic':
      x0 = np.array([a])

    return x0

  def prep_output(self, params, cs, output_type='unit_cell'):
    if cs == 'Triclinic':
      xopt = params
    elif cs == 'Monoclinic':
      if output_type == 'unit_cell':
        xopt = np.array([params[0],params[1],params[2],90,params[3],90])
      else:
        xopt = np.array([params[0],params[1],params[2],0,params[3],0])
    elif cs == 'Orthorhombic':
      if output_type == 'unit_cell':
        xopt = np.array([params[0],params[1],params[2],90,90,90])
      else:
        xopt = np.array([params[0],params[1],params[2],0,0,0])
    elif cs == 'Tetragonal':
      if output_type == 'unit_cell':
        xopt = np.array([params[0],params[0],params[1],90,90,90])
      else:
        xopt = np.array([params[0],params[0],params[1],0,0,0])
    elif cs == 'Trigonal' or cs == 'Hexagonal':
      if output_type == 'unit_cell':
        xopt = np.array([params[0],params[0],params[1],90,90,120])
      else:
        xopt = np.array([params[0],params[0],params[1],0,0,0])
    elif cs == 'Cubic':
      if output_type == 'unit_cell':
        xopt = np.array([params[0],params[0],params[0],90,90,90])
      else:
        xopt = np.array([params[0],params[0],params[0],0,0,0])

    return xopt

  def func(self, params, *args):
    I_r = args[0]
    observations_hdlr = args[1]
    const_params = args[2]
    refine_mode = args[3]
    cs = observations_hdlr.observations_original.crystal_symmetry().space_group().crystal_system()

    if refine_mode == 'scale_factor':
      G = params[0]
      b11, b22, b33, b12, b13, b23 = self.prep_output(params[1:], cs, output_type='b_tensor')
      rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = const_params
    elif refine_mode == 'all_params':
      rotx, roty, ry, rz, r0, re = params[:6]
      G, b11, b22, b33, b12, b13, b23 = const_params
      uc_params = self.prep_input(params[6:], cs)
      a, b, c, alpha, beta, gamma = self.prep_output(uc_params, cs)

    try:
      uc = unit_cell((a,b,c,alpha,beta,gamma))
    except Exception:
      return None

    #refresh observations handler
    observations_hdlr.set_params(postref_params=[G, b11, b22, b33, b12, b13, b23, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma])
    obs_full = observations_hdlr.calc_full_observations()

    return (I_r - obs_full.data())

  def func_scale(self, params, *args):
    mean_I_r = args[0]
    mean_I_o = args[1]
    mean_stol_sq = args[2]
    G, B = params
    mean_I_o_scaled = mean_I_o/(G * flex.exp(flex.double(-2 * B * mean_stol_sq)))
    return (mean_I_r - mean_I_o_scaled)

  def good_unit_cell(self, uc_params):
    flag_good_uc = False
    uc_tol = 5
    if (abs(uc_params[0]-self.iparams.target_unit_cell.parameters()[0]) \
        <= (uc_tol*self.iparams.target_unit_cell.parameters()[0]/100) \
                  and abs(uc_params[1]-self.iparams.target_unit_cell.parameters()[1]) \
                  <= (uc_tol*self.iparams.target_unit_cell.parameters()[1]/100) \
                  and abs(uc_params[2]-self.iparams.target_unit_cell.parameters()[2]) \
                  <= (uc_tol*self.iparams.target_unit_cell.parameters()[2]/100) \
                  and abs(uc_params[3]-self.iparams.target_unit_cell.parameters()[3]) \
                  <= (uc_tol*self.iparams.target_unit_cell.parameters()[3]/100) \
                  and abs(uc_params[4]-self.iparams.target_unit_cell.parameters()[4]) \
                  <= (uc_tol*self.iparams.target_unit_cell.parameters()[4]/100) \
                  and abs(uc_params[5]-self.iparams.target_unit_cell.parameters()[5]) \
                  <= (uc_tol*self.iparams.target_unit_cell.parameters()[5]/100)):
      flag_good_uc = True
    return flag_good_uc


  def optimize_scalefactors(self, observations_hdlr_ref, observations_hdlr, use_binning=False):
    """
    Keep other partiality parameters constant and refine only scale factors.
    """
    observations_hdlr_sel = observations_handler(self.iparams)
    observations_hdlr_sel.copy_from(observations_hdlr)
    observations_hdlr_sel.resolution_filter(d_min=self.iparams.d_min, d_max=self.scale_0_d_max)
    observations_hdlr_sel.outlier_filter(sigma_I_cutoff=self.sigma_I_cutoff)

    refine_mode = 'scale_factor'
    if observations_hdlr.postref_params is None:
      a,b,c,alpha,beta,gamma = observations_hdlr.observations.unit_cell().parameters()
      postref_params = [1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,self.iparams.gamma_e,a,b,c,alpha,beta,gamma]
    else:
      postref_params = observations_hdlr.postref_params[:]

    if use_binning:
      #use binning
      obs_ref = observations_hdlr_ref.observations.deep_copy()
      obs_ref = obs_ref.resolution_filter(d_min=self.iparams.d_min, d_max=self.scale_0_d_max)
      binner = obs_ref.setup_binner(n_bins=30)
      mean_I_ref = obs_ref.mean(
        use_binning=True,
        use_multiplicities=True).data[1:-1]

      obs = observations_hdlr_sel.calc_full_observations()
      obs.use_binning_of(obs_ref)
      mean_I_obs = obs.mean(
        use_binning=True,
        use_multiplicities=True).data[1:-1]

      stol_sq = obs_ref.sin_theta_over_lambda_sq()
      stol_sq.use_binning_of(obs_ref)
      mean_stol_sq = stol_sq.mean(
        use_binning=True,
        use_multiplicities=True).data[1:-1]

      mean_I_ref_sel = flex.double()
      mean_I_obs_sel = flex.double()
      mean_stol_sq_sel = flex.double()
      for i in range(len(mean_I_obs)):
        if mean_I_ref[i] is not None and mean_I_obs[i] is not None and mean_stol_sq[i] is not None:
          mean_I_ref_sel.append(mean_I_ref[i])
          mean_I_obs_sel.append(mean_I_obs[i])
          mean_stol_sq_sel.append(mean_stol_sq[i])

      xinp = np.array([1,0])
      xopt, cov_x, infodict, mesg, ier = optimize.leastsq(self.func_scale, xinp,
                                                            args=(mean_I_ref_sel, mean_I_obs_sel, mean_stol_sq_sel, None),
                                                            full_output=True)
      G, biso = xopt
      b11, b22, b33, b12, b13, b23 = (biso,biso,biso,0,0,0)
    else:
      #use common reflections
      obs_hdlr_ref, obs_hdlr = observations_hdlr_ref.common_filter(observations_hdlr_sel)
      I_r = obs_hdlr_ref.observations.data()

      cs = observations_hdlr.observations_original.crystal_symmetry().space_group().crystal_system()
      G, b11, b22, b33, b12, b13, b23, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = observations_hdlr.postref_params
      const_params = [rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma]
      b_tensor_indep = self.prep_input([b11, b22, b33, b12, b13, b23], cs)
      xinp = np.append(G, b_tensor_indep)
      xopt, cov_x, infodict, mesg, ier = optimize.leastsq(self.func, xinp,
                                                            args=(I_r, obs_hdlr, const_params, refine_mode),
                                                            full_output=True)
      G = xopt[0]
      b11, b22, b33, b12, b13, b23 = self.prep_output(xopt[1:], cs, output_type='b_tensor')

    postref_params[:7] = [G,b11, b22, b33, b12, b13, b23]
    observations_hdlr.set_params(postref_params=postref_params)

    cc_init, cc_final, r_init, r_final, r_xy_init, r_xy_final, cc_iso_init, cc_iso_final  = (0,0,0,0,0,0,0,0)
    if True:
      obs_ref, obs = observations_hdlr_ref.observations.common_sets(observations_hdlr.observations,
        assert_is_similar_symmetry=False)
      obs_ref, obs_full = observations_hdlr_ref.observations.common_sets(observations_hdlr.calc_full_observations(),
        assert_is_similar_symmetry=False)
      cc_init = obs_ref.correlation(other=obs, use_binning=False, assert_is_similar_symmetry=False).coefficient()
      cc_final = obs_ref.correlation(other=obs_full, use_binning=False, assert_is_similar_symmetry=False).coefficient()
      r_init = flex.sum((obs_ref.data()-obs.data())**2)/flex.sum(obs_ref.data()**2)
      r_final = flex.sum((obs_ref.data()-obs_full.data())**2)/flex.sum(obs_ref.data()**2)
    #except Exception:
    #  pass
    observations_hdlr.set_params(stats=[cc_init, cc_final, r_init, r_final, r_xy_init, r_xy_final, cc_iso_init, cc_iso_final])


  def optimize(self, observations_hdlr_ref, observations_hdlr):
    """
    Refine all parameters
    """
    refine_mode = 'all_params'
    #get matching observations_handler
    obs_hdlr_ref, obs_hdlr = observations_hdlr_ref.common_filter(observations_hdlr)
    I_r = obs_hdlr_ref.observations.data()

    cs = observations_hdlr.observations_original.crystal_symmetry().space_group().crystal_system()
    G, b11, b22, b33, b12, b13, b23, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = observations_hdlr.postref_params
    const_params = [G, b11, b22, b33, b12, b13, b23]
    partiality_params = [rotx, roty, ry, rz, r0, re]
    cs = observations_hdlr.observations_original.crystal_symmetry().space_group().crystal_system()
    uc_indep = self.prep_input(uc_params, cs)
    xinp = np.append(partiality_params, uc_indep)
    xopt, cov_x, infodict, mesg, ier = optimize.leastsq(self.func, xinp,
                                                          args=(I_r, obs_hdlr, const_params, refine_mode),
                                                          full_output=True, maxfev=100)
    CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final  = (0,0,0,0,0,0,0,0)
    uc_all = self.prep_output(xopt[6:], cs)
    if self.good_unit_cell(uc_all):
      postref_params = observations_hdlr.postref_params[:]
      xopt_all = np.append(xopt[:6], uc_all)
      postref_params[7:] = xopt_all
      observations_hdlr.set_params(postref_params=postref_params)
