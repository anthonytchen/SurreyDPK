# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:41:45 2017

@author: tc0008
"""

from importlib import reload
import chemical
import config
reload(config)
import viaepd
import dermis
import skin_setup
reload(skin_setup)

import paras
reload(paras)


comps_setup = 'H,V,S,E,D'

_paras = paras.InputParas(10)
_paras.comps_paras = [paras.StraCornParas() for i in range(10)] 

_paras.comps_paras[0]
# read chemical config file
# read compartment structure file
# set up Kw/D etc 
# then set up skin to start simulation

_conf = config.Config('Caffeine_test.cfg')
_chem = chemical.Chemical(_conf)

_skin = skin_setup.Skin_Setup(_chem, _conf)
_skin.createComps(_chem, _conf)
_skin.solveMoL(0, 1)
_skin.comps[0].saveMeshConc(True, 'tmp.txt')
_skin.comps[1].saveMeshConc(False, 'tmp.txt')
_skin.comps[2].saveMeshConc(False, 'tmp.txt')
_skin.comps[3].saveMeshConc(False, 'tmp.txt')
_skin.comps[4].saveMeshConc(False, 'tmp.txt')

#_vpd = viaepd.ViaEpd(10, 10, 1, 2, 2, None, None)
#_vpd.createMesh( _chem, 0, 0)