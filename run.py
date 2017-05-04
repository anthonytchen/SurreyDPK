# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:41:45 2017

@author: tc0008
"""
import numpy as np
import matplotlib.pyplot as plt
from importlib import reload

import chemical
import config
reload(config)
import viaepd
import dermis
import skin_setup
reload(skin_setup)


# read chemical config file
# read compartment structure file
# set up Kw/D etc 
# then set up skin to start simulation

#_conf = config.Config('Caffeine_test.cfg')
_conf = config.Config('Nicotine_test.cfg')
_chem = chemical.Chemical(_conf)

_skin = skin_setup.Skin_Setup(_chem, _conf)
_skin.createComps(_chem, _conf)

t_start = 0
t_end = 36000
Nsteps = 61

mass = np.zeros([Nsteps, 4])

t_range = np.linspace(t_start, t_end, Nsteps)
b_first = True

mass[0,:] = _skin.compMass_comps()
print('Time = ', t_range[0], 'Flux vh_sc= ', '{:.4e}'.format(_skin.compFlux([0,0], 3)[0]), \
       'Flux sc_down=', '{:.4e}'.format(_skin.compFlux([1,0], 3)[0]) )
for i in range(Nsteps-1):
    
    _skin.solveMoL(t_range[i], t_range[i+1])
    if b_first:
        _skin.comps[0].saveMeshConc(True, 'tmp.txt')
        b_first = False
    else:
        _skin.comps[0].saveMeshConc(False, 'tmp.txt')
    _skin.comps[1].saveMeshConc(False, 'tmp.txt')
    _skin.comps[2].saveMeshConc(False, 'tmp.txt')
    _skin.comps[3].saveMeshConc(False, 'tmp.txt')
    
    mass[i+1,:] = _skin.compMass_comps()
    print('Time = ', t_range[i+1], 'Flux vh_sc= ', '{:.4e}'.format(_skin.compFlux([0,0], 3)[0]), \
       'Flux sc_down=', '{:.4e}'.format(_skin.compFlux([1,0], 3)[0]) )

plt.plot(t_range/60, mass[:,0]/mass[0,0], \
         t_range/60, mass[:,1]/mass[0,0], \
         t_range/60, mass[:,2]/mass[0,0], \
         t_range/60, mass[:,2]/mass[0,0])
plt.show()

#_skin.comps[2].saveMeshConc(False, 'tmp.txt')
#_skin.comps[3].saveMeshConc(False, 'tmp.txt')
#_skin.comps[4].saveMeshConc(False, 'tmp.txt')

#_vpd = viaepd.ViaEpd(10, 10, 1, 2, 2, None, None)
#_vpd.createMesh( _chem, 0, 0)