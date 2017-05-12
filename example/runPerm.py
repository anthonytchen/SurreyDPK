# -*- coding: utf-8 -*-
"""
A module containing files for calculating steady-state permeability
    and related optimisation for parameter estimation
"""
import os
import numpy as np
from scipy.optimize import minimize, basinhopping, brute
import matplotlib.pyplot as plt
from importlib import reload


# Import various skin 
from core import chemical
reload(chemical)
from core import config
reload(config)
from core import viaepd
reload(viaepd)
from core import dermis
reload(dermis)
from core import skin_setup
reload(skin_setup)

def compPerm(fn_conf) :
    """Compute steady state permeability
    """
    # Read the .cfg, i.e. configuration, file to set up simulation
    _conf = config.Config(fn_conf)
    _conf.D_sc_paras[:2] = np.array([9.47, 9.32e-8]) 

    # Setup the chemical
    _chem = chemical.Chemical(_conf)

    # Setup skin and create compartments
    _skin = skin_setup.Skin_Setup(_chem, _conf)
    _skin.createComps(_chem, _conf)

    # Simulation time (in seconds) and steps
    t_start, t_end, Nsteps = [0, 3600*48, 101]
    t_range = np.linspace(t_start, t_end, Nsteps)

    for i in range(Nsteps):
        flux_vh_sc = -_skin.compFlux([0,0], 3)[0]
        flux_sc_down = -_skin.compFlux([1,0], 3)[0]
        if np.fabs( (flux_vh_sc-flux_sc_down) / flux_vh_sc ) < 1e-3 :
            break
        elif i == Nsteps-1 :
            raise ValueError('Simulation time too short to reach steady-state; re-run the simulation with longer time.')
        
        print('Time = ', t_range[i], 'Flux vh_sc= ', '{:.3e}'.format(flux_vh_sc), \
              'Flux sc_down=', '{:.3e}'.format(flux_sc_down) )
    
        # Simulate
        _skin.solveMoL(t_range[i], t_range[i+1])
        
    return flux_vh_sc / _conf.init_conc_vh



def compPermObj(paras) :
    """Compute the objective function for optimisation
    The objective function is to minimise the least square error
    between model prediction and data
    """
    #def objFunCalib(paras, t_range, blood_conc_data, dose_factor=1.0):

def calibPerm():
    """Calibration of permeability model by adjusting parameters
    """
    
    # optimising the depth of vehicle
    #bnds = ((-0.3, 0), (50, 200))
    #x0 = [0, 100];
    #res_brute = brute(objFunCalib, bnds, args=(t_range, blood_conc_data, dose_factor), Ns = 50, finish = None, full_output=True, disp=True)
    #res = minimize(objFunCalib, res_brute[0], args=(t_range, blood_conc_data, dose_factor), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 100})
    # Options for the optimiser: SLSQP, L-BFGS-B, TNC
