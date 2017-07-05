# -*- coding: utf-8 -*-
"""
A module containing files for calculating kinetics
    for in-vitro experiments
"""
import os
import numpy as np
#from scipy.optimize import minimize, basinhopping, brute
import matplotlib.pyplot as plt
from importlib import reload

# Using multiple processers
N_PROCESS = 8

# Import various skin 
from core import chemical
reload(chemical)
from core import config
reload(config)
from core import vehicle
reload(vehicle)
from core import viaepd
reload(viaepd)
from core import dermis
reload(dermis)
from core import skin_setup
reload(skin_setup)

def compDPK(fn_conf, chem=None, disp=1) :
    """Compute DPK
    Args:
        fn_conf -- the .cfg file, which gives the configuration of the simulation
        chem -- if given, it overrides the values given in fn_conf
    """
    # Read the .cfg, i.e. configuration, file to set up simulation
    _conf = config.Config(fn_conf)

    # Setup the chemical
    if chem is not None:
        _chem = chem
    else:
        _chem = chemical.Chemical(_conf)

    # Setup skin and create compartments
    _skin = skin_setup.Skin_Setup(_chem, _conf)
    _skin.createComps(_chem, _conf)
    #print(_skin.comps[0].meshes[0].Kw, _skin.comps[0].meshes[0].D)
    #print(_skin.comps[1].meshes[0].Kw, _skin.comps[1].meshes[0].D)    
    #print(_skin.comps[1].meshes[11].Kw, _skin.comps[1].meshes[11].D)    

    # Simulation time (in seconds) and steps
    #t_start, t_end, Nsteps = [0, 3600*25, 26]
    t_start, t_end, Nsteps = [0, 3600, 61]
    #t_start, t_end, Nsteps = [0, 3600, 3601]    
    t_range = np.linspace(t_start, t_end, Nsteps)  
    #t_range = np.r_[np.linspace(0, 1350, 2), np.linspace(1380, 1440, 61),\
    #                np.linspace(1450, 3600, 44)]
    #Nsteps = len(t_range)
    
    nComps = _skin.nxComp*_skin.nyComp
    total_mass = np.sum( _skin.compMass_comps() )
    for i in range(Nsteps-1):
        
        mass = _skin.compMass_comps()
        m_v = _skin.comps[0].getMass_OutEvap()
        #total_mass = np.sum(mass)
        mass_down = _skin.comps[0].compMassIrregMeshDown(_skin.comps[0].meshes[0], _skin.comps[0].meshes[0].conc)
        
        if disp >= 2:
            np.set_printoptions(precision=4)
            print('Time = ', t_range[i], ', total mass: ', total_mass, '% mass: ', m_v/total_mass, mass/total_mass )
            #print('\t K_vw = ', _skin.comps[0].meshes[0].Kw, ', % D_v = ', _skin.comps[0].meshes[0].D )
            #print('\t mass_down = ', mass_down )
            
        # Create directory to save results
        newpath = './simu/' + str(t_range[i])
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        
        # Save current concentrations        
        for j in range(nComps):
            fn = newpath + '/comp' + str(j) + '_' + _conf.comps_geom[j].name        
            _skin.comps[j].saveMeshConc(True, fn)
            
        # Simulate
        _skin.solveMoL(t_range[i], t_range[i+1])
        
        
    return mass


