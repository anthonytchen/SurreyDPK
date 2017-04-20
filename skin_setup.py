# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 09:57:44 2017

@author: tc0008
"""
import importlib
#import numpy as np
import stracorn
importlib.reload(stracorn)
import viaepd
importlib.reload(viaepd)
import dermis
importlib.reload(dermis)
import skin
importlib.reload(skin)

class Skin_Setup(skin.Skin):
    """ Class definition for Skin_Setup
    which intends to set up the compartments in simulation as instructed by user
    """
    def __init__(self, chem, conf):
        skin.Skin.__init__(self)
        self.comp_structure = conf.sComps
        
        if chem is None:
            skin.Skin.set_n_species(self, 0)
        elif type(chem) is not list :
            skin.Skin.set_n_species(self, 1)
        else :
            skin.Skin.set_n_species(self, len(chem)) # Number of chemcial species
        
        if conf.sComps.find('B') == 1 and conf.sComps.find('D') == 1 : 
            self.b_has_blood = True # Has blood compartment ('B') in dermis ('D')
        else :
            self.b_has_blood = False
    
    def createComps(self, chem, conf) :
        """ Create compartments
        Letter code:
            V: vehicle            S: stratum cornuem
            E: viable epidermis   D: dermis
            B: blood              H: Hair
        """
        
        # Read structure of the compartments
        tokens = list( filter(None, self.comp_structure.split(',')) )
        nrow = len(tokens)
        ncol = len(tokens[0])
        skin.Skin.createComps(self, nrow, ncol)        
        
        # Actually create the compartments
        
        dim_all = 0
        current_x = 0
        current_y = 0        
        
        for i in range(nrow) :
            
            for j in range(ncol) :
                
                # determine the boundary conditions for this compartment
                bdy_up = 'ZeroFlux' if i==0 else 'FromOther'
                #if i==0 : 
                #    bdy_up = 'ZeroFlux'
                #else:
                #    bdy_up = 'FromOther'
                bdy_down = 'ZeroConc' if i==nrow-1 else 'FromOther'
                bdy_left = 'Periodic'
                bdy_right = 'Periodic'
                assert(ncol==1) # otherwise the above periodic conditions are incorrect
                bdy_cond = [bdy_up, bdy_left, bdy_right, bdy_down]                
                 
                if tokens[i][j] == 'S':
                    comp = self.createSC(chem, current_x, current_y, conf.n_layer_x_sc, conf.n_layer_y_sc, 
                                         conf.offset_y_sc, bdy_cond)
                elif tokens[i][j] == 'E':
                    comp = self.createVE(chem, current_x, current_y, conf.x_len_ve, conf.y_len_ve, 
                                         conf.n_grids_x_ve, conf.n_grids_y_ve, bdy_cond)
                elif tokens[i][j] == 'D':
                    comp = self.createDE(chem, current_x, current_y, conf.x_len_de, conf.y_len_de, 
                                         conf.n_grids_x_de, conf.n_grids_y_de, bdy_cond)
                else :
                    pass
                skin.Skin.setComp(self, comp, i, j)
                dim_all += comp.get_dim()
                current_y += comp.get_y_length()
                
            current_x += comp.get_x_length()
            current_y = 0
            
        skin.Skin.set_dim_all(self, dim_all)
        
        # Link compartments through boundary conditions
        for i in range(nrow) :
            for j in range(ncol) :
                
                # down boundary
                if i == nrow-1 : # down-most compartmnet, its down boundary is zero-flux
                    n_dBdy = 0
                    mesh_dBdy = None
                else :
                    compDown = skin.Skin.getComp(self, i+1, j)
                    n_dBdy = compDown.get_ny()
                    mesh_dBdy = compDown.meshes[0:n_dBdy]
                    
                # right boundary
                if j==ncol-1 : # right-most compartment, its right boundary is zero-flux
                    n_rBdy = 0
                    mesh_rBdy = None
                else :
                    compRight = skin.Skin.getComp(self, i, j+1)
                    n_rBdy = compRight.get_nx()
                    mesh_rBdy = compRight.meshes[::compRight.get_ny()]
                    
                comp = skin.Skin.getComp(self, i, j)
                comp.createBdy(n_rBdy, n_dBdy)
                comp.setBdyMesh(mesh_rBdy, mesh_dBdy)
                #print('n_rBdy=', n_rBdy, ' n_dBdy=', n_dBdy)
                    
        
    ### (START OF) Create individual compartments ###
    
    def createSC(self, chem, 
                 coord_x_start, coord_y_start, n_layer_x, n_layer_y, offset_y,
                 bdyCond) :
        """ Create stratum corneum """                
        sc = stracorn.StraCorn(self.dz_dtheta, n_layer_x, n_layer_y, offset_y, self.coord_sys, bdyCond)
        sc.createMesh(chem, coord_x_start, coord_y_start)
        return sc
        
    def createVE(self, chem, 
                 coord_x_start, coord_y_start, xlen, ylen, n_grids_x, n_grids_y,
                 bdyCond) :
        """ Create viable epidermis """        
        viable_epidermis = viaepd.ViaEpd(xlen, ylen, self.dz_dtheta, 
                                         n_grids_x, n_grids_y, self.coord_sys, bdyCond)
        viable_epidermis.createMesh(chem, coord_x_start, coord_y_start)
        return viable_epidermis
    
    def createDE(self, chem, 
                 coord_x_start, coord_y_start, xlen, ylen, n_grids_x, n_grids_y,
                 bdyCond) :
        """ Create dermis """        
        derm = dermis.Dermis(xlen, ylen, self.dz_dtheta, 
                             n_grids_x, n_grids_y, self.coord_sys, bdyCond, self.b_has_blood)
        derm.createMesh(chem, coord_x_start, coord_y_start)
        return derm
        
    ### (END OF) Create individual compartments ###