# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:57:38 2017

@author: tc0008
"""
import numpy as np
import importlib
#import scipy as sp
from scipy.integrate import ode
import chemical
import comp
importlib.reload(comp)

class Skin:
    """ Class definition for Skin
    which can contain a number of compartments
    """
    
    def __init__(self, coord_sys = None, dz_dtheta = None, bdy_cond = None):
        """      """        
        self.comps = None # List of compartments
        self.dz_dtheta = 0.01
        self.coord_sys = 'Cartesian'
        self.nxComp = 0
        self.nyComp = 0
        self.dim_all = 0
        self.nSpecies = 0
        
    def createComps(self, nrow, ncol):
        self.nxComp = nrow
        self.nyComp = ncol
        self.comps = [comp.Comp() for i in range(nrow*ncol)]
        
    def set_n_species(self, n):
        self.nSpecies = n
    def set_dim_all(self, dim):
        self.dim_all = dim
    def get_dim_all(self):
        return self.dim_all
    def setComp(self, comp, idx_comp_x, idx_comp_y) :
        self.comps[ idx_comp_x*self.nyComp + idx_comp_y ] = comp
    def getComp(self, idx_comp_x, idx_comp_y) :
        return self.comps[ idx_comp_x*self.nyComp + idx_comp_y ]
        
    ### (START OF) Class methods dealing with ODE computation ###
    
    def compODEdydt_diffu (self, t, y, args=None):
        """Compute the right-hand side of the ODEs, i.e. dydt, due to diffusion
        """
        f = np.zeros(y.shape) #
        #print(y)
        
        for k in range(self.nSpecies) :
            idx = k*self.dim_all
            
            ## 1. reset mass transfer between compartments
            for i in range(self.nxComp) :
                for j in range(self.nyComp) :
                    idx_comp = i*self.nyComp+j
                    current_comp = self.comps[idx_comp]
                    current_comp.setBdyMassInOutZero()
                    
            ## 2. Calculate the right hand side of the differential equations
            for i in range(self.nxComp) :
                for j in range(self.nyComp) :
                    idx_comp = i*self.nyComp+j
                    current_comp = self.comps[idx_comp]

                    # 2.1. Update the concentration in boundary meshes according to y[]                    
                    compBdyRight, concBdyRight = self.getBdyRight(current_comp, y, idx, i, j)
                    compBdyDown, concBdyDown = self.getBdyDown(current_comp, y, idx, i, j)
                    current_comp.setBdyConc(concBdyRight, concBdyDown)
                    
                    # 2.2. Call the compartment specific ODE functions
                    #   for certain compartments special treatment is needed
                    #if type(current_comp).__name__ == 'Dermis' :
                        #if (self.b_has_blood) # y[ k*m_dim_all+m_dim_all-1 ] contains the blood concentration, i.e. the last term in the differential equations
                        #((Dermis *)pComp)->updateBlood( y[ k*m_dim_all+m_dim_all-1 ] );
                    tmp_f = current_comp.compODEdydt_diffu(t, y+idx)
                    f[idx:idx+current_comp.get_dim()] = tmp_f
                    #print(y)

                    #if (m_bInfSrc)                              // infinite source, concentration doesn't change, but above dydt calculation is still needed
                    #memset(f+idx, 0, sizeof(double)*((Vehicle *)pComp)->m_dim);  //  since calling compODE_dydt will calculate the flux across boundaries properly
                    current_comp.passBdyMassOut(compBdyRight, compBdyDown)
                    idx += current_comp.get_dim()
                # for j
            # for i
            
            # Simulation is for a small skin area, but needs to multiple
            #  the mass transport due to blood flow by the actual topical application area
            #if (m_b_has_blood){
          
            #double factor = m_Vehicle_area / m_Dermis[k].compTotalArea(0);

            #todo: when more than one dermis compartments are involved, need to collate mass in-out of all dermis compartments

            #m_Blood[k].updateMassInOutDermis(m_Dermis[k].m_mass_into_dermis, m_Dermis[k].m_mass_outof_dermis, factor);
            #m_Blood[k].compODE_dydt(t, y+idx, f+idx);
        # for k
        
        return f
    
    def solveMoL(self, t_start, t_end) :
        """ Solving PDE using method of lines (MoL)
        """
  
        ## get current concentration, and set as initial conditions
        y0 = np.zeros(self.dim_all*self.nSpecies)
        for k in range(self.nSpecies) :
            idx = k*self.dim_all
            
            for i in range(self.nxComp) :
                for j in range(self.nyComp) :
                    idx_comp = i*self.nyComp+j
                    current_comp = self.comps[idx_comp]
                    dim = current_comp.get_dim()
                    y0[ idx : idx+dim ] = current_comp.getMeshConc()
            #if (m_b_has_blood)
            #  m_Blood[k].getGridsConc(y+idx, m_Blood[k].m_dim);
        # for k, each species

        ## Integration
        
        r = ode(self.compODEdydt_diffu).set_integrator('vode', method='bdf', with_jacobian=False)
        r.set_initial_value(y0, t_start).set_f_params(None)
        r.integrate( r.t + t_end-t_start )
        
        ## Extract calculated values 
        for k in range(self.nSpecies) :
            idx = k*self.dim_all
            
            for i in range(self.nxComp) :
                for j in range(self.nyComp) :
                    idx_comp = i*self.nyComp+j
                    current_comp = self.comps[idx_comp]
                    dim = current_comp.get_dim()
                    current_comp.setMeshConc_all( r.y[ idx : idx+dim ] )
            #    if (m_b_has_blood)
            # m_Blood[k].setGridsConc(pNVs+idx, m_Blood[k].m_dim);
        # for k, each species
    
    ### (END OF) Class methods dealing with ODE computation ###
    
    
    ### (START OF) Class methods dealing with boundaries ###
    
    def getBdyRight(self, compThis, y, idx_up2now, cIdx_i, cIdx_j):
        """ Get the boundary to the right of this compartment
        using values in y[]
        """
        if (cIdx_j == self.nyComp-1) : # rightmost compartment, no right boundary
            assert( compThis.get_bdyCond(2) != 'FromOther' )
            compBdyRight = None
            size = 0            
            concBdyRight = None
        else :            
            compBdyRight = self.comps[ cIdx_i*self.nyComp + cIdx_j + 1 ]
            size = compBdyRight.get_nx()

            idx = idx_up2now + compThis.get_dim()
            concBdyRight = np.zeros(size)
            for i in range(size):
                concBdyRight[i] = y[ idx + i*compBdyRight.get_ny() ]            
            
        return (compBdyRight, concBdyRight)
    
        
    def getBdyDown(self, compThis, y, idx_up2now, cIdx_i, cIdx_j):
        """ Get the boundary to the down of this compartment
        using values in y[]
        """
        if cIdx_i == self.nxComp-1 : # downmost compartment, no down boundary
            assert( compThis.get_bdyCond(3) != 'FromOther' )
            compBdyDown = None
            size = 0
            concBdyDown = None
        else :
            compBdyDown = self.comps[ (cIdx_i+1)*self.nyComp + cIdx_j ]
            size = compBdyDown.get_ny()
            
            # work out the index for the down boundary
            idx = idx_up2now
            
            i = cIdx_i
            j = cIdx_j
            while j < self.nyComp :
                idx += self.comps[ i*self.nyComp + j ].get_dim()
                j += 1
                
            i = cIdx_i+1
            j = 0
            while j < cIdx_j :
                idx += self.comps[ i*self.nyComp + j ].get_dim()
                j += 1
            
            # Fill in concBdyDown from y[]
            concBdyDown = np.zeros(size)
            for j in range(size):
                concBdyDown[j] = y[ idx + j ]            
            
        return (compBdyDown, concBdyDown)   
           
      
    ### (END OF) Class methods dealing with boundaries ###