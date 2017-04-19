# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:41:45 2017

@author: tc0008
"""

import chemical
import config
import viaepd
import skin_setup
reload(skin_setup)

_conf = config.Config('Caffeine_test.cfg')
_chem = chemical.Chemical(_conf)

_skin = skin_setup.Skin_Setup(_chem, _conf)
_skin.createComps(_chem, _conf)
_skin.solveMoL(0, 1)
_skin.comps[0].saveMeshConc(True, 'tmp.txt')

#_vpd = viaepd.ViaEpd(10, 10, 1, 2, 2, None, None)
#_vpd.createMesh( _chem, 0, 0)