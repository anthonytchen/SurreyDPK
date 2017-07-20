# -*- coding: utf-8 -*-
"""
A module containing files for calculating kinetics
    for in-vitro experiments
"""
import os, sys
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

def compDPK(fn_conf, chem=None, sc_Kw_paras=None, sc_D_paras=None, disp=1, wk_path='./simu/') :
    """Compute DPK
    Args:
        fn_conf -- the .cfg file, which gives the configuration of the simulation
        chem -- if given, it overrides the values given in fn_conf
        wk_path -- path to save simulation results
    """
    # Read the .cfg, i.e. configuration, file to set up simulation
    _conf = config.Config(fn_conf)
    if sc_Kw_paras is not None:
        _conf.Kw_sc_paras = sc_Kw_paras
    if sc_D_paras is not None:
        _conf.D_sc_paras = sc_D_paras  
        
    # Setup the chemical
    if chem is not None:
        _chem = chem
    else:
        _chem = chemical.Chemical(_conf)

    # Setup skin and create compartments
    _skin = skin_setup.Skin_Setup(_chem, _conf)
    _skin.createComps(_chem, _conf)     

    # Simulation time (in seconds) and steps
    t_start, t_end, Nsteps = [0, 3600*24, 145]
    #t_start, t_end, Nsteps = [0, 1800, 181]
    t_range = np.linspace(t_start, t_end, Nsteps)  
    #t_range = np.r_[np.linspace(0, 1000, 2), np.linspace(1200, 1800, 21),\
    #                np.linspace(1800, 3600, 21),np.linspace(7200, 3600*24, 23)]
    #Nsteps = len(t_range)
    
    nComps = _skin.nxComp*_skin.nyComp
    total_mass = np.sum( _skin.compMass_comps() )
    
    # Create directory to save results
    newpath = wk_path
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    fn = wk_path + 'MassFrac.csv'
    saveMass(total_mass, fn, b_1st_time=True)    
    
    for i in range(Nsteps):
        
        mass = _skin.compMass_comps()
        m_v = _skin.comps[0].getMass_OutEvap()
        m_all = np.insert(mass, 0, m_v) / total_mass
        
        if disp >= 2:
            np.set_printoptions(precision=2)
            print('Time = ', t_range[i], '% mass: ', m_all)            
            
        # Create directory to save results
        newpath = wk_path + str(t_range[i])
        if not os.path.exists(newpath):
            os.makedirs(newpath)
            
        # Save fraction of mass in all compartments
        fn = wk_path + 'MassFrac.csv'
        saveMass(np.insert(m_all, 0, t_range[i]), fn)
        
        # Save current concentrations        
        for j in range(nComps):
            fn = newpath + '/comp' + str(j) + '_' + _conf.comps_geom[j].name        
            _skin.comps[j].saveMeshConc(True, fn)
        
        if i == Nsteps-1:
            break
        
        # Simulate
        _skin.solveMoL(t_range[i], t_range[i+1])    
    
    #return mass

def saveMass(nparray, fn, b_1st_time=False) :
        """ Save mass and fractions to file
        Args: 
            nparray -- the data to be saved
            b_1st_time -- if True, write to a new file; otherwise append to the existing file
        """
        if b_1st_time :
            file = open(fn, 'w')
        else :
            file = open(fn, 'a')
        
        if type(nparray) is np.ndarray:
            nd = len(nparray)
            for i in range(nd):
                file.write("{:.6e}".format(nparray[i]))            
                if i<nd-1:
                    file.write(",")
                #print("nparray", file)
        else:
            file.write( "{:.6e}".format(nparray) )
        file.write('\n')
        file.close()

        
def compDPK_KwVar(fn_conf, sc_Kw_paras=None, wk_path='./simu/') :        
    """ Simulate with the uncertainty in Kw for lipid and corneocytes """
    _conf = config.Config(fn_conf)
    _chem = chemical.Chemical(_conf)        
        
    # Run uncertainty quantification 
    import example.runUnctQSPR as rUq
    reload(rUq)
    #print(np.array([_chem.mw]))
    Kw_var, Kw_base = rUq.compUnct_Ksc(np.array([[_chem.mw]]))
    print('Base case, Klp= ', 10**Kw_base[0][0][1], ' Kcc= ', 10**Kw_base[0][0][0])
    #print(Kw_var[0][0])
    m_cc, m_lp = Kw_var[0][0]
    print('Re-calibrated, Klp= ', 10**m_lp, ' Kcc= ', 10**m_cc)
    sd_cc, sd_lp = np.sqrt( np.diag(Kw_var[1][0]) ) 
    
    # Set up KwDParas to be passed to simulation
    from core.stracorn import KwDParas
    sc_Kw_paras = KwDParas()
    sc_Kw_paras.lp.option = 'VALE'
    sc_Kw_paras.cc.option = 'VALE'
    
    N = 50
    for i in range(N):  # N MC samples
        sc_Kw_paras.lp.value = np.array( [10**np.random.normal(m_lp, sd_lp)] )
        sc_Kw_paras.cc.value = np.array( [10**np.random.normal(m_cc, sd_cc)] )
        wk_path_i = wk_path + 'rep_' + str(i) + '/'
        print('\t Rep ', i, 'Klp= ', sc_Kw_paras.lp.value, ' Kcc= ', sc_Kw_paras.cc.value, '\n')
        compDPK(fn_conf, _chem, sc_Kw_paras, disp=3, wk_path=wk_path_i)