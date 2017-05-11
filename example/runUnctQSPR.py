"""
A module containing functions to calculate uncertainty in QSPR models
"""

import numpy as np
from importlib import reload

from qspr import stracorn, stracorn_cc, stracorn_lp
from uncertainty import hybmdl
reload(stracorn)
reload(stracorn_cc)
reload(stracorn_lp)

dat_Plp = np.loadtxt("qspr/Kow_Plp.txt")
dat_Kcc = np.loadtxt("qspr/Kow_Kcc.txt")
dat_Ksc = np.loadtxt("qspr/Kow_Ksc.txt")

paras0 = np.log( np.array([4.2, 0.31, 0.69]) )
bnds = ((-10, 10), (-10, 10), (-10, 10))

sig2_y = np.array([0.05])
sig2_z = np.array([0.05, 0.05])

Xy = dat_Ksc[:,0].reshape((-1, 1))
Y = dat_Ksc[:,1].reshape((-1,1))
Xz = (dat_Plp[:,0].reshape((-1,1)), dat_Kcc[:,0].reshape((-1,1)))
Z = (dat_Plp[:,1].reshape((-1,1)), dat_Kcc[:,1].reshape((-1,1)))

paras = np.empty_like (paras0)
np.copyto(paras, paras0)

rlt_plugin = hybmdl.PluginMain(qspr_K_sc_plugin, qspr_K_cc_lp, Xy, Y, Xz, Z, paras0, sig2_y, sig2_z, 10, bnds)
rlt_pred_plugin = hybmdl.pred(qspr_K_sc, qspr_K_cc_lp, Xy, rlt_plugin[0], rlt_plugin[1], rlt_plugin[2], rlt_plugin[3])

# first variable in Z
#rlt_pred_plugin_z0 = hybmdl.pred(qspr_K_sc, qspr_K_cc_lp, Xz[0], rlt_plugin[0], rlt_plugin[1], rlt_plugin[2], rlt_plugin[3])

#plt_x = Xz[0][:,0]
#plt_dat = Z[0][:,0]
#plt_pred = np.squeeze(rlt_pred_plugin_z0[0][:,0])
#plt_pred_h = plt_pred + 1.58*np.squeeze(np.sqrt(rlt_pred_plugin_z0[1][:,0,0]))
#plt_pred_l = plt_pred - 1.58*np.squeeze(np.sqrt(rlt_pred_plugin_z0[1][:,0,0]))
#plt.plot(plt_x, plt_dat, 'x', plt_x, plt_pred, 'o', plt_x, plt_pred_h, '.', plt_x, plt_pred_l, '.')


#plt.plot(Xy[:,0], np.squeeze(rlt_pred_plugin[2]), 'x', Xy[:,0], Y, 'o', Xy[:,0], np.squeeze(rlt_pred_plugin[2])+1.58*np.squeeze(np.sqrt(rlt_pred_plugin[3])), '.', Xy[:,0], np.squeeze(rlt_pred_plugin[2])-1.58*np.squeeze(np.sqrt(rlt_pred_plugin[3])), '.')


def qspr_K_cc_lp(theta, Kow):
    ''' Function to predict the volumetric partition coefficient of corneocyte:water and lipid:water
    Here used as a combined function of the LOW-LEVEL of the multi-level model
    Args:
      theta -- model parameters
      Kow, dim: [n_dat, 1]
    Rtns:
      Z -- dim: [n_dat, 2]; each row contains [Kcc, Klp]
    '''

    theta = np.squeeze(np.asarray(theta)) # to avoid some mysterious conversion of np array to np matrix by the optimiser        
    paras_cc = theta[:2]
    paras_lp = theta[2]

    n_dat = Kow.shape[0]    #print X.shape
    Z = np.zeros((n_dat, 2))
    Z[:,0] = stracorn_cc.compK(paras_cc, Kow)
    Z[:,1] = stracorn_lp.compK(paras_lp, Kow)         
    return Z

    
def qspr_K_sc(theta, Kow, Z):
    ''' Function to predict the volumetric partition coefficient between stratum corneum and water    
    Here used as a function of the TOP-LEVEL of the multi-level model
    Args:
      theta -- model parameters
      Kow, dim: [n_dat, 1]
      Z -- Kcc & Klp, dim: [n_dat, 2]
    Rtns:
      Y -- Ksc, predicted partition coefficient between stratum corneum and water
    '''
    Y = stracorn.compK_from_cc_lp(theta, Kow, Z)
    #Y1 = stracorn.compK(theta, Kow)
    return Y#(Y, Y1)

def qspr_K_sc_plugin(theta, Kow, func_low):
    ''' Function to predict the volumetric partition coefficient between stratum corneum and water    
    Here used as a test function of the top-level model
    Args:
      theta -- model parameters
      Kow, dim: [n_dat, 1]
    Rtns:
      Y -- Ksc_pred, predicted coefficient between stratum corneum (and water)
    '''
    K_cc_lp = func_low(theta, Kow)
    Y = qspr_K_sc(theta, Kow, K_cc_lp)
    return Y
