# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:28:52 2017

@author: tc0008
"""
import importlib
#import numpy as np

from core import comp
importlib.reload(comp)


class Vehicle(comp.Comp):
    """Class definition for Vehicle
    which is the delivery vehicle, currently modelled as a homogenised media
    """
    
    def __init__(self, xlen, ylen, dz_dtheta, nx, ny, init_conc, Kw, D,
                 coord_sys, bdy_cond, b_inf_source=False):
        comp.Comp.__init__(self)
        comp.Comp.setup(self, xlen, ylen, dz_dtheta, nx, ny, coord_sys, bdy_cond)
        
        self.eta = 7.644E-4 # water viscosity at 305 K (32 deg C) (Pa s)
        self.b_inf_source = b_inf_source
        
        self.init_conc = init_conc
        comp.Comp.set_Kw(self, Kw)
        comp.Comp.set_D(self, D)

        
    def createMesh(self, chem, coord_x_start, coord_y_start) :
        """ Create mesh for this compartment
        Args:
                coord_x_start, coord_y_start: starting coordinates
        """
        self.compParDiff(chem)
        comp.Comp.createMeshHomo(self, 'VH', chem, self.init_conc, coord_x_start, coord_y_start)
        
        
    def compParDiff(self, chem) :
        """ Compute the partition coefficient with respect to water
        and the diffusion coefficient
        """
        if self.Kw < 0:
            Kw = 1 # caution: only placeholder and needs refining
            comp.Comp.set_Kw(self, Kw)
        
                    
        if self.D < 0: # calculation of diffusivity according to the Stoke-Eistein equation
            D = comp.Comp.compDiff_stokes(self, self.eta, chem.r_s)
            comp.Comp.set_D(self, D)
        
        
        #return (Kw, D)
                

    def compODEdydt(self, t, y, args=None):
        """ The wrapper function for computing the right hand side of ODEs
        """
        dydt = comp.Comp.compODEdydt_diffu (self, t, y, args)
        
        # If infinite source, concentration doesn't change, but above dydt calculation 
        #   is still needed since calling compODEdydt_diffu will calculate the 
        #   flux across boundaries properly
        #print(self.b_inf_source)
        #exit
        if self.b_inf_source :
            dydt.fill(0)            
            
        return dydt
        
    def saveCoord(self, fn_x, fn_y) :
        comp.Comp.saveCoord(self, fn_x, fn_y, '.vh')
