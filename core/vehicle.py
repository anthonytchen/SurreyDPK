# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:28:52 2017

@author: tc0008
"""
import importlib, sys
import numpy as np

from core import comp
importlib.reload(comp)


class Vehicle(comp.Comp):
    """Class definition for Vehicle
    which is the delivery vehicle, currently modelled as a homogenised media
    """
    
    def __init__(self, xlen, ylen, dz_dtheta, nx, ny, init_conc, Kw, D,
                 coord_sys, bdy_cond, b_inf_source=False,
                 rho_solute=1e3, rho_solvent=1e3, phase_solute='LIQUID',
                 k_evap_solvent=0, k_evap_solute=0, solubility=1e10):
        comp.Comp.__init__(self)
        comp.Comp.setup(self, xlen, ylen, dz_dtheta, nx, ny, coord_sys, bdy_cond)
        
        self.eta = 7.644E-4 # water viscosity at 305 K (32 deg C) (Pa s)
        self.b_inf_source = b_inf_source
        
        # evaporative mass transfer coefficient for solvent and solute
        self.b_vary_vehicle = True
        self.k_evap_solvent = k_evap_solvent
        self.k_evap_solute = k_evap_solute
        self.solubility = solubility
        self.rho_solute = rho_solute
        self.rho_solvent = rho_solvent
        self.phase_solute = phase_solute
        
        self.conc_solvent = rho_solvent
        self.depth_vehicle = xlen
        self.mass_out_evap = 0
        self.mass_out_phase = 0
        
                
        self.init_conc = init_conc
                
        comp.Comp.set_Kw(self, Kw)
        comp.Comp.set_D(self, D)
        
    def getMass_OutEvap(self):
        return self.mass_out_evap
    def getMass_OutPhase(self):
        return self.mass_out_phase
        
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
        if self.b_vary_vehicle is False:
            dydt = comp.Comp.compODEdydt_diffu (self, t, y, args)
            
            # If infinite source, concentration doesn't change, but above dydt calculation 
            #   is still needed since calling compODEdydt_diffu will calculate the 
            #   flux across boundaries properly
            #print(self.b_inf_source)
            #exit
            if self.b_inf_source :
                dydt.fill(0)            

        else :
            # Evaporation of solvent and solute
            #   Currently a crude approximation
            #   and only implemented for a homogeneous vehicle compartment                
            
            dim = self.nx*self.ny
        
            assert(self.nx==1 and self.ny==1)
            # y[0] - solute conc, y[1] - solvent conc, y[2] - h,
            # y[3] - solute mass out due to evarporation
            # y[4] - solute mass out due to over-solubility
            
            # todo: evaporation is a function of the total liquid concentration
            #  not just the concentration in the solution (but also in separate solute phase)
            tmp = self.k_evap_solute/self.rho_solute*y[0] \
                    + self.k_evap_solvent/self.rho_solvent*y[1]

            dhdt = - tmp
            dy0dt = ( y[0]*tmp - y[0]*self.k_evap_solute ) / y[2]
            dy1dt = ( y[1]*tmp - y[1]*self.k_evap_solvent ) / y[2]
            dy3dt = y[0]*self.k_evap_solute*self.compTotalArea(3)
            dy4dt = 0
            
            if y[0]>self.solubility:
                k_phase_out = 1e-10
                c_diff = y[0] - self.solubility
                V = self.compTotalArea(3) * y[2]
                dy4dt = k_phase_out * c_diff
                #print('dy4dt=', dy4dt, 'dy0dt=', dy0dt)
                dy0dt -= k_phase_out * c_diff / V
            
            if self.phase_solute is 'LIQUID':
                # update partiton coefficient
                
                # partition and volume of the active's liquid phase
                P1 = self.rho_solute/self.solubility
                V1 = y[4] / self.rho_solute

                # partition and volume of the solvent phase
                P2 = self.Kw
                V2 = self.compTotalArea(3) * y[2]

                K_vw = P1*V1/(V1+V2) + P2*V2/(V1+V2)
                self.setMeshes_Kw(K_vw)

            else: # undissolved solute is solid
                min_h = y[4] / self.rho_solute / self.compTotalArea(3)
                if y[2]<min_h : # too little liquid, stop diffusion
                    self.setMeshes_Kw(1e10)
            
            #print('y=', y, 'V=', self.compTotalArea(3) * y[2])
            #print('dydt=', [dy0dt, dy1dt, dhdt, dy3dt, dy4dt])
            #sys.exit()
            self.meshes[0].dx = y[2] # this is needed to use the correct vehicle depth
            dydt = comp.Comp.compODEdydt_diffu (self, t, y, args)                     
            
            dydt += np.array([dy0dt, dy1dt, dhdt, dy3dt, dy4dt])
            #print('dydt=', dydt)
            #sys.exit()

            #following are tmp code for solids; to be completed
            #total_mass = self.init_conc * self.compTotalVolume()
            #min_h = total_mass / density / self.compTotalArea(3)
            #if h > min_h :
            #    if y[0]<self.solubility :
            #        dydt_evap = ( self.k_evap_solvent - self.k_evap_solute ) * y[0] / h
            #        dydt[0] += dydt_evap
            #else: # turn off vehicle by setting a very small diffusivity
            #    self.setMeshes_D(1e-100)
            #print(y, dydt, dydt_evap)
        
        return dydt
        
    def saveCoord(self, fn_x, fn_y) :
        comp.Comp.saveCoord(self, fn_x, fn_y, '.vh')
        
        
    def getMeshConc(self) :
        """ This function name is a misnomer but meant to be consistent with 
        the same function in class comp
        The function returns the concentration from all meshes
        AND also the variables due to varying vehicle
        into a single numpy array
        """
        if self.b_vary_vehicle is False:
            return comp.Comp.getMeshConc(self)
        else:
            dim = self.dim
            y = np.zeros(self.get_dim())
            y[:dim] = comp.Comp.getMeshConc(self)
            y[dim:] = [self.conc_solvent, self.depth_vehicle,\
                       self.mass_out_evap, self.mass_out_phase]
            return y

    def setMeshConc_all(self, conc) :
        """ Similar to the above, this function name is a misnomer but meant to be 
        consistent with the same function in class comp
        The function sets the concentration for all meshes
        AND also the variables due to varying vehicle
        """
        if self.b_vary_vehicle is False:
            comp.Comp.setMeshConc_all(self,conc)
        else:
            dim = self.dim
            comp.Comp.setMeshConc_all(self,conc[:dim])
            self.conc_solvent, self.depth_vehicle, \
                self.mass_out_evap, self.mass_out_phase = conc[dim:]
                        
            assert(self.nx==1 and self.ny==1)
            self.x_length = self.depth_vehicle
            self.meshes[0].dx = self.x_length

    def get_dim(self):
        if self.b_vary_vehicle is False:
            return comp.Comp.get_dim(self)
        else:
            return comp.Comp.get_dim(self)+4

    def saveMeshConc(self, b_1st_time, fn) :
        """ Save mesh concentrations to file
        Args: b_1st_time -- if True, write to a new file; otherwise append to the existing file
        """
        comp.Comp.saveMeshConc(self, b_1st_time, fn)
        if self.b_vary_vehicle is True:            
            file = open(fn, 'a')
            file.write( "{:.6e}\n".format( self.conc_solvent ) )
            file.write( "{:.6e}\n".format( self.depth_vehicle ) )
            file.write( "{:.6e}\n".format( self.mass_out_evap ) )
            file.write( "{:.6e}\n".format( self.mass_out_phase ) )
            file.close()
            