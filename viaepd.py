# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:28:52 2017

@author: tc0008
"""
import importlib
import numpy as np

import mesh
import comp
importlib.reload(comp)


class ViaEpd(comp.Comp):
    """Class definition for ViaEpd
    which is the viable epidermis, currently modelled as a homogenised media
    """
    
    def __init__(self, xlen, ylen, dz_dtheta, nx, ny, coord_sys, bdy_cond):
        comp.Comp.__init__(self)
        comp.Comp.setup(self, xlen, ylen, dz_dtheta, nx, ny, coord_sys, bdy_cond)
        #print('good')
        
        #self.b_homo_media = True # whether the compartment should be treated as homogeneous media
        
        
    def createMesh(self, chem, coord_x_start, coord_y_start) :
        """ Create mesh for this compartment
        Args:
                coord_x_start, coord_y_start: starting coordinates
        """
        init_conc = .0
        self.compParDiff(chem)
        comp.Comp.createMeshHomo(self, 'VE', chem, init_conc, coord_x_start, coord_y_start)
        self.meshes[0].setConc(1)

        
    def compParDiff(self, chem) :
        """ Compute the partition coefficient with respect to water
        and the diffusion coefficient
        """
        rou_lipid = 0.9
        rou_water = 1
        K_lip = rou_lipid / rou_water * chem.get_K_ow()**0.69 # Partition in SC lipid
  
        binding_factor = 0.65 + 0.325/chem.get_frac_unbound() + 0.025*chem.get_frac_non_ion()*K_lip

        comp.Comp.set_Kw(self, 0.6 * binding_factor)
        
        # c.f. L. Chen's Phar. Res. paper; -8.15 used because of unit (m2/s)
        #   Kasting's original paper used -4.15 because of unit (cm2/s)
        D_free = 10 ** ( -8.15-0.655*np.log10(chem.get_mw()) )
        comp.Comp.set_D(self, D_free / binding_factor)
        #print('self.D = ', self.D)

        
    def saveCoord(self, fn_x, fn_y) :
        comp.Comp.saveCoord(self, fn_x, fn_y, '.ve')
