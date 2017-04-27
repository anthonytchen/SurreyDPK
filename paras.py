# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 10:46:13 2017

@author: tc0008
"""
import importlib

import stracorn
importlib.reload(stracorn)

class CompParas:
    """ Parameters generic to compartments """
    def __init__(self, xlen=None, ylen=None, n_meshes_x=None, n_meshes_y=None, 
                 Kw=-1, D=-1, init_conc=0):
        self.xlen = xlen
        self.ylen = ylen
        self.n_meshes_x = n_meshes_x
        self.n_meshes_y = n_meshes_y
        self.Kw = Kw
        self.D = D
        self.init_conc = init_conc

class StraCornParas(CompParas):
    """ Parameters specific to stratum corneum """
    def __init__(self, n_layer_x, n_layer_y, n_meshes_x_lp, n_meshes_y_lp, offset_y=0,
                 Kw=None, D=None, init_conc=None):
        self.n_layer_x = n_layer_x
        self.n_layer_y = n_layer_y
        self.n_meshes_x_lp = n_meshes_x_lp
        self.n_meshes_y_lp = n_meshes_y_lp
        self.offset_y = offset_y
        
        # The only purpose of running this is to calculate the dimensions and number of
        #   meshes in SC from the supplied information
        sc = stracorn.StraCorn(None, n_layer_x, n_layer_y, offset_y, None, None)
        
        CompParas.__init__(self, sc.get_x_length(), sc.get_y_length(), sc.get_nx(), sc.get_ny(), 
                           Kw, D, init_conc)        
    
class InputParas:
    def __init__(self, ncomps):
        self.ncomps = ncomps
        self.comps_paras = None