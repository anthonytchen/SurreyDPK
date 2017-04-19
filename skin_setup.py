# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 09:57:44 2017

@author: tc0008
"""
import importlib
import numpy as np
import viaepd
importlib.reload(viaepd)
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
                if tokens[i][j] == 'E':
                    comp = self.createVE(chem, current_x, current_y, conf.x_len_ve, conf.y_len_ve, 
                                         conf.n_grids_x_ve, conf.n_grids_y_ve, ['ZeroFlux']*4)
                skin.Skin.setComp(self, comp, i, j)
                dim_all += comp.get_dim()
                current_y += comp.get_y_length()
                
            current_x += comp.get_x_length()
            current_y = 0
            
        skin.Skin.set_dim_all(self, dim_all)
        
        # Link compartments through boundary conditions
        for i in range(nrow) :
            for j in range(ncol) :
                if tokens[i][j] == 'E':
                    comp = skin.Skin.getComp(self, i, j)
                    comp.createBdy(0, 0)
                    comp.setBdyMesh(None, None)
                    
        
    ### (START OF) Create individual compartments ###
    
    def createVE(self, chem, 
                 coord_x_start, coord_y_start, xlen, ylen, n_grids_x, n_grids_y,
                 bdyCond) :
        """ Create viable epidermis """        
        viable_epidermis = viaepd.ViaEpd(xlen, ylen, self.dz_dtheta, 
                                         n_grids_x, n_grids_y, self.coord_sys, None)
        viable_epidermis.createMesh(chem, coord_x_start, coord_y_start)
        return viable_epidermis
    
    ### (END OF) Create individual compartments ###