# -*- coding: utf-8 -*-
'''
Created on Tue Apr 18 22:01:17 2017

@author: tc0008
'''

import warnings

class Config:
    ''' This class intends to read from a text configuration
    file the properties of the chemical(s), the setup of the skin compartments 
    and their dimensions, and the simulation parameters,
    and then feeds them to the appropriate classes (e.g. Chemical, Skin etc.)
    '''
    
    def __init__(self, fn_config):
        """         """
        ## Default values
        self.y_len_ve = 40e-6
        self.n_grids_y_ve = 2
        self.y_len_de = self.y_len_ve
        self.y_len_vehicle = self.y_len_ve
        self.n_grids_y_de = 2
		
        with open(fn_config, 'r') as f:
            lines = f.readlines()
        f.close()
        for lin in lines:
            tokens = list( filter(None, lin.split()) )
            if len(tokens) != 0 :
                self.readTokens(tokens)
        
    
    def readTokens(self, tokens):
        
        #print(tokens)

        if tokens[0][0] == '#' : # comments, ignore
            return            
        
        # setup compartments
        elif tokens[0] == 'COMPARTMENT_SETUP': 
            self.sComps = tokens[1]

        ### parameters relating to chemical(s)
        elif tokens[0] == 'CHEM_NO' : # number of compounds
            self.nChem = int(tokens[1])
        elif tokens[0] == 'CHEM_MW' : # molecular weight
            self.mw = float(tokens[1])
        elif tokens[0] == 'CHEM_KOW' : # partition coefficient between octanol and water
            self.K_ow = float(tokens[1])
        elif tokens[0] == 'CHEM_PKA' : # pKa -- acide dissociation constant
            self.pKa = float(tokens[1])
        elif tokens[0] == 'CHEM_NONION' : # fraction of solute non-ionised at pH 7.4
            self.frac_non_ion = float(tokens[1])
        elif tokens[0] == 'CHEM_UNBND' : # fraction of solute unbound in a 2.7% albumin solution at pH 7.4    
            self.frac_unbound = float(tokens[1])
        elif tokens[0] == 'CHEM_ACIDBASE' :
            self.acid_base = tokens[1]

        ### vehicle specific parameters
        elif tokens[0] == 'INFINITE_VH' : # 
            self.bInfSrc = bool(tokens[1])
        elif tokens[0] == 'AREA_VH' :
            self.area_vehicle = float(tokens[1])
        
        ### Inital conditions
        elif tokens[0] == 'INIT_CONC_VH' : # vehicle
            self.init_conc_vh = float(tokens[1])
        elif tokens[0] == 'INIT_CONC_SC' : # stratum corneum
            self.init_conc_sc = float(tokens[1])
        elif tokens[0] == 'INIT_CONC_VE' : # viable epidermis
            self.init_conc_ve = float(tokens[1])
        elif tokens[0] == 'INIT_CONC_DE' : # dermis
            self.init_conc_de = float(tokens[1])
        elif tokens[0] == 'INIT_CONC_HF' : # hair follicle
            self.init_conc_hf = float(tokens[1])
        elif tokens[0] == 'INIT_CONC_BD' : # blood
            self.init_conc_bd = float(tokens[1])
            
            
        # partition and diffusion coefficients
        elif tokens[0] == 'VEH_INIT_CONC' : # initial concentration in vehicle
            self.conc_vehicle = float(tokens[1])
        elif tokens[0] == 'VEH_DX' : 
            self.x_len_vehicle = float(tokens[1])
        elif tokens[0] == 'VEH_AREA' : 
            self.area_vehicle = float(tokens[1])
        elif tokens[0] == 'VEH_INFINITE' :
            self.bInfSrc = bool(tokens[1])
        elif tokens[0] == 'N_GRIDS_X_VH' :
            self.n_grids_x_vh = int(tokens[1])
        elif tokens[0] == 'N_GRIDS_Y_VH' :
            self.n_grids_y_vh = int(tokens[1])
        
        elif tokens[0] == 'SKIN_N_LAYER_X_SC' :
            self.n_layer_x_sc = int(tokens[1])
        elif tokens[0] == 'SKIN_N_LAYER_Y_SC' :
            self.n_layer_y_sc = int(tokens[1])
        elif tokens[0] == 'SKIN_OFFSET_Y_SC' :
            self.offset_y_sc = float(tokens[1])
            
        elif tokens[0] == 'SKIN_N_GRIDS_X_VE' :
            self.n_grids_x_ve = int(tokens[1])
        elif tokens[0] == 'SKIN_N_GRIDS_X_DE' :
            self.n_grids_x_de = int(tokens[1])
        elif tokens[0] == 'SKIN_LEN_X_VE' :
            self.x_len_ve = float(tokens[1])
        elif tokens[0] == 'SKIN_LEN_X_DE' :
            self.x_len_de = float(tokens[1])
        elif tokens[0] == 'SKIN_PARTITION_DE2BD' : # dermis to blood partition
            self.partition_dermis2blood = float(tokens[1])
        elif tokens[0] == 'SKIN_CLEAR_BD' : # blood clearance rate
            self.Kclear_blood = float(tokens[1])
            
        elif tokens[0] == 'SKIN_N_GRIDS_X_SB_SUR' : 
            self.n_grids_x_sb_sur = int(tokens[1])
        elif tokens[0] == 'SKIN_LEN_X_SB_SUR' : 
            self.x_len_sb_sur = float(tokens[1])
        elif tokens[0] == 'SKIN_N_GRIDS_Y_SB_SUR' :
            self.n_grids_y_sb_sur = int(tokens[1])
        elif tokens[0] == 'SKIN_N_GRIDS_X_SB_HAR' :
            self.n_grids_x_sb_har = int(tokens[1])
        elif tokens[0] == 'SKIN_N_GRIDS_Y_SB_HAR' :
            self.n_grids_y_sb_har = int(tokens[1])
        elif tokens[0] == 'SKIN_LEN_Y_SB_HAR' :
            self.y_len_sb_har = float(tokens[1])


        # name not found
        else :
            warnings.warn('Unrecognised line in config file')
            
    def assignComps(self):
        #tokens = list( filter(None, self.comp_structure.split(',')) )
        #nrow = len(tokens)
        #ncol = len(tokens[0])