# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:52:24 2017

@author: tc0008
"""

import importlib
import numpy as np

import mesh
importlib.reload(mesh)
import comp
importlib.reload(comp)
import point

class StraCorn(comp.Comp):
    """ Class definition for StraCorn, 
    which is the stratum corneum, currently modelled as a heterogeneous media
    """
    
    def __init__(self, dz_dtheta, n_layer_x, n_layer_y, offset_y, coord_sys, bdy_cond) :
        
        comp.Comp.__init__(self)
          
        # dimension related; c.f. Readme.docx for more details
        self.w = 8.0 # offset ratio, 8.0
        self.nx_grids_lipid = 1 # Number of x-grids for lipid layer, 2
        self.nx_grids_cc = 1    # Number of x-grids for corneocyte layer, 4
        self.ny_grids_lipid = 1 # Number of y-grids for lipid layer, 2
        self.ny_grids_cc_dn = 1 # Number of y-grids for dn-part of the offset corneocyte layer, 2
        
        self.T = 309 # temperature (Kelvin)
        self.eta = 7.1E-4 # water viscosity at above temperature (Pa s),  
        
        self.geom_g = 0.075e-6  
        self.geom_d = 40e-6   
        self.geom_s = 0.075e-6  
        self.geom_t = 0.8e-6 

        self.geom_dm = self.w*(self.geom_d-self.geom_s)/(1+self.w)
        self.geom_dn = self.geom_d-self.geom_s-self.geom_dm
	
        # Vertical direction, lipid layer is at both top and bottom of the stratum corneum
        nx = (self.nx_grids_lipid+self.nx_grids_cc)*n_layer_x + self.nx_grids_lipid
        # Lateral direction, [dh] [s] [dm] [s], here d=dh+dm+s, w=dm/dh
        ny = int( ( self.ny_grids_lipid*2 + self.ny_grids_cc_dn + self.ny_grids_cc_dn*self.w ) * n_layer_y )
        
        xlen = n_layer_x*(self.geom_g+self.geom_t)+self.geom_g
        ylen = n_layer_y*(self.geom_d+self.geom_s)
        
        comp.Comp.setup(self, xlen, ylen, dz_dtheta, nx, ny, coord_sys, bdy_cond)

        # For the volume fraction calculation, we assume a Cartesian coordinate
        #   and the z-direction width m_dz_dtheta is directly used.
        #   This won't affect if other coordinates are used
        self.V_mortar = ( self.geom_g*(self.geom_d+self.geom_s)+self.geom_t*self.geom_s ) * dz_dtheta
        self.V_brick = self.geom_d*self.geom_t * dz_dtheta
        self.V_all = self.V_mortar + self.V_brick

        self.offset_y =  offset_y
        
        
    
    def createMesh(self, chem, coord_x_start, coord_y_start, water_frac_surface=.55) :
        """ Create mesh for this compartment
        Args:
                coord_x_start, coord_y_start: starting coordinates
                water_frac_surface : water content (w/w); saturation at .55; dry skin at .25
        """
        #init_conc = .0
        #self.compParDiff(chem)
        #comp.Comp.createMeshHomo(self, 'VE', chem, init_conc, coord_x_start, coord_y_start)
               
        # some initial settings
        bOffset = False
        cc_subtype = 0 # 0 = d_n; 1 = s; 2 = d_m, 3 = s;
  
        dx_lipid = self.geom_g / self.nx_grids_lipid
        dx_cc = self.geom_t / self.nx_grids_cc
        dy_lipid = self.geom_s / self.ny_grids_lipid
        dy_cc = self.geom_dn / self.ny_grids_cc_dn
        
        self.meshes = [mesh.Mesh() for i in range(self.nx*self.ny)] # organised in row dominant               
        
        # work out the starting point from given offset
        
        length = self.geom_d + self.geom_s
        len_vec = [self.geom_dn, self.geom_s, self.geom_dm, self.geom_s]

        while self.offset_y > length :
            self.offset_y -= length
        length = .0
        for i in range(4) :
            length += len_vec[i]
            if self.offset_y < length :
                break
        cc_subtype_offset = i
        length = self.offset_y - length + len_vec[i]

        if cc_subtype_offset == 0 or cc_subtype_offset == 2 :
            # corneocyte width
            idx_y_offset = int( length / dy_cc )
            dy_offset = dy_cc
        elif cc_subtype_offset == 1 or cc_subtype_offset == 3 :
            # lipid width
            idx_y_offset = int( length / dy_lipid )
            dy_offset = dy_lipid
        else :
            raise ValueError('Invalid subtype name')


        # Now prepare values for the loop to create meshes
        
        water_frac_sat = 0.55 # saturated water content (w/w)
        water_increment_per_x = (water_frac_sat - water_frac_surface) / self.x_length
    
        idx_x = 0
        idx_y = idx_y_offset
        cc_subtype = cc_subtype_offset

        coord_x = coord_x_start
        coord_y = coord_y_start
        # starting from lipid on the top layer
        current_point = point.Point(coord_x, coord_y, dx_lipid, dy_offset, 'LP', 'LP')
        init_conc = 0
        
        for i in range(self.nx) : # verticle direction up to down
            water_frac = water_frac_surface + (current_point.x_coord - coord_x_start) * water_increment_per_x

            for j in range(self.ny) : # lateral direction left to right
                idx = i*self.ny + j
                
                # assign type
                if current_point.x_type == 'LP' or current_point.y_type == 'LP' :
                    # Entire lipid layer (1st ==) or lateral lipid between two coreneocytes (2nd ==)
                    name = 'LP'
                else : # corneocyte
                    name = 'CC'
                #   and then create meshes
                Kw, D = self.compParDiff(name, chem, water_frac, water_frac_sat, 
                        self.V_mortar, self.V_brick, self.V_all, self.T, self.eta)                    
                self.meshes[idx].setup(name, chem, init_conc, Kw, D,
                        current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, 
                        self.dz_dtheta)
                
                ### update current_point, fairly complicated
                
                if j==self.ny-1 : # last element in the lateral direction, move down
                
                    idx_x += 1
                    coord_y = coord_y_start
                    
                    if current_point.x_type == 'LP' : 
                        coord_x += dx_lipid
                        
                        if idx_x == self.nx_grids_lipid : # reaching end of lateral lipid layer
                            if cc_subtype_offset == 0 or cc_subtype_offset == 2 :
                                current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'CC')
                            elif cc_subtype_offset == 1 :
                                if not bOffset:
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'CC')
                                else :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'LP')
                            elif cc_subtype_offset == 3 :
                                if not bOffset:
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'LP')
                                else :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'CC')
                            else :
                                raise ValueError('Invalid subtype name')	  
                            idx_x = 0 
                        else : # NOT reaching end of lateral lipid layer
                            current_point.setPoint(coord_x, coord_y, dx_lipid, dy_offset, 'LP', 'LP')
                            
                    elif current_point.x_type == 'CC' :
                        coord_x += dx_cc
                        
                        if idx_x == self.nx_grids_cc : # reaching end of lateral corneocyte layer
                            current_point.setPoint(coord_x, coord_y, dx_lipid, dy_offset, 'LP', 'LP')
                            idx_x = 0
                            bOffset = not bOffset
                        else : # NOT reaching end of lateral corneocyte layer
                            if cc_subtype_offset == 0 or cc_subtype_offset == 2 :
                                current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'CC')
                            elif cc_subtype_offset == 1 :
                                if not bOffset :   
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'CC')
                                else :            
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'LP')
                            elif cc_subtype_offset == 3 :
                                if not bOffset :   
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'LP')
                                else :           
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_offset, 'CC', 'CC')
                            else :
                                raise ValueError('Invalid subtype name')	
                                
                    else :
                        raise ValueError('Invalid current_point.x_type')
                            
                    idx_y = idx_y_offset
                    cc_subtype = cc_subtype_offset

                else : # not the last element in the lateral direction, thus move to the right
                
                    idx_y += 1
				
                    if current_point.x_type == 'LP' : # current row is lipid
                    
                        if cc_subtype == 0 : # now within dh
                            coord_y += dy_cc
                            if idx_y == self.ny_grids_cc_dn :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_lipid, 'LP', 'LP')
                                idx_y = 0
                                cc_subtype += 1
                            else :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_cc, 'LP', 'LP')
                        elif cc_subtype == 1 : # now within s
                            coord_y += dy_lipid
                            if idx_y == self.ny_grids_lipid :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_cc, 'LP', 'LP')
                                idx_y = 0
                                cc_subtype += 1
                            else :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_lipid, 'LP', 'LP')
                        elif cc_subtype == 2 : # now wtihin dm
                            coord_y += dy_cc
                            if idx_y == self.ny_grids_cc_dn*self.w :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_lipid, 'LP', 'LP')
                                idx_y = 0
                                cc_subtype += 1
                            else :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_cc, 'LP', 'LP')
                        elif cc_subtype == 3 : # now within the 2nd s
                            coord_y += dy_lipid
                            if idx_y == self.ny_grids_lipid :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_cc, 'LP', 'LP')
                                idx_y = 0
                                cc_subtype = 0
                            else :
                                current_point.setPoint(coord_x, coord_y, dx_lipid, dy_lipid, 'LP', 'LP')
                        else :
                            raise ValueError('Invalid cc_subtype')
                        
                    elif current_point.x_type == 'CC' : # current row is corneocyte
                    
                        if cc_subtype == 0 : # now within dh
                            coord_y += dy_cc
                            if idx_y == self.ny_grids_cc_dn :
                                if bOffset :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'LP')
                                else :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'CC')
                                idx_y = 0
                                cc_subtype += 1
                            else :
                                current_point.setPoint(coord_x, coord_y, dx_cc, dy_cc, 'CC', 'CC')
                        elif cc_subtype == 1 : # now within s
                            coord_y += dy_lipid
                            if idx_y == self.ny_grids_lipid :
                                current_point.setPoint(coord_x, coord_y, dx_cc, dy_cc, 'CC', 'CC')
                                idx_y = 0
                                cc_subtype += 1
                            else :
                                if bOffset :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'LP')
                                else :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'CC')
                        elif cc_subtype == 2 : # now wtihin dm
                            coord_y += dy_cc
                            if idx_y == self.ny_grids_cc_dn*self.w :
                                if bOffset :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'CC')
                                else :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'LP')
                                idx_y = 0
                                cc_subtype += 1
                            else :
                                current_point.setPoint(coord_x, coord_y, dx_cc, dy_cc, 'CC', 'CC')
                        elif cc_subtype == 3 : # now within the 2nd s
                            coord_y += dy_lipid
                            if idx_y == self.ny_grids_lipid :
                                current_point.setPoint(coord_x, coord_y, dx_cc, dy_cc, 'CC', 'CC')
                                idx_y = 0
                                cc_subtype = 0
                            else :
                                if bOffset :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'CC')
                                else :
                                    current_point.setPoint(coord_x, coord_y, dx_cc, dy_lipid, 'CC', 'LP')
                        else :
                            raise ValueError('Invalid cc_subtype')
                            
                    else :
                        raise ValueError('Invalid current_point.x_type')     
            # for j
        # for i


        
    def compParDiff(self, name, chem, mass_frac_water, mass_frac_water_sat, 
                    V_mortar, V_brick, V_all, T, eta) :
        """ Compute the partition coefficient with respect to water
        and the diffusion coefficient
        """
        rou_lipid = 0.9
        rou_keratin = 1.37
        rou_water = 1
        
        # 2. calulcate volume fraction of water in CC based on
        #    mass fraction of water in SC
        theta_b = self.compVolFracWater_cc( mass_frac_water, rou_lipid, rou_keratin, rou_water,
				 V_mortar, V_brick, V_all)
        # do the same for saturated water
        phi_b = self.compVolFracWater_cc( mass_frac_water_sat, rou_lipid, rou_keratin, rou_water,
				 V_mortar, V_brick, V_all)


        # 3. calculate diffusivity and partition coefficient

        K = 1.3806488 * 1E-23 # Boltzmann constant, Kg m^2 s^{-2}
        r_f = 3.5e-9 # keratin microfibril radius, 3.5 nm
        Dw = K*T/6/np.pi/eta/chem.r_s # diffusivity in water, Stoke-Eistein equation
        
        if name == 'LP' :
            # m_Kw = pow(m_K_ow, 0.7) # this is the old model
            Kw = rou_lipid / rou_water * (chem.K_ow ** 0.69)
            
            if chem.mw <= 380.0 :
                r_s_inA = chem.r_s*1e10 # unit in Angstrom
                D = 2 * 1E-9 * np.exp(-0.46*r_s_inA*r_s_inA)
            else : 
                D = 3 * 1E-13
        elif name == 'CC' : # corneocyte
            # the following not used:
            #   if (m_K_ow>10)
            #    K_kw = 5.6 * pow(m_K_ow, 0.27);
            #    else 
            #    K_kw = 0.5* ( 1 + pow(m_K_ow, 0.7) );
            K_kw = rou_keratin / rou_water * 4.2 * (chem.K_ow**0.31)
            Kw = (1-phi_b) * K_kw + theta_b

            # empirically fitted parameters
            alpha = 9.47
            beta = 9.32 * 1E-8
            lambdaa = 1.09
            gamma = -1.17
            
            r_s_inA = chem.r_s*1e10 # unit in A
            r_f_inA = r_f*1e10 # unit in A
            
            phi_f = 1 - theta_b
            k = beta*r_f_inA*r_f_inA* (phi_f**gamma)
            S = (r_s_inA+r_f_inA)/r_f_inA
            S = phi_f * S*S
            
            D = np.exp( -alpha*(S**lambdaa) ) / ( 1 + r_s_inA/np.sqrt(k) + r_s_inA*r_s_inA/3/k )
            D *= Dw
        
        #comp.Comp.set_Kw(self, Kw)
        #comp.Comp.set_D(self, Dw)
        
        return (Kw, D)
        
    def compVolFracWater_cc(self, mass_frac_water, rou_lipid, rou_keratin, rou_water, 
                            V_mortar_geometry, V_brick_geometry, V_all_geometry):
        """ Compute the volume fraction of water in corneocyte
            based on the water content (mass fraction of water) of the stratum corneum
        """
        f_l = 0.125 # dry mass fraction of SC lipid and keratin  
        f_k = 1 - f_l
        
        # mass fraction of lipid and keratin
        mass_lipid = (1 - mass_frac_water) * f_l
        mass_keratin = (1 - mass_frac_water) * f_k  

        V_all = mass_lipid/rou_lipid + mass_keratin/rou_keratin + mass_frac_water/rou_water
        V_lipid = mass_lipid/rou_lipid / V_all
        V_keratin = mass_keratin/rou_keratin / V_all  
        
        V_water_mortar = V_mortar_geometry/V_all_geometry - V_lipid
        V_water_brick = V_brick_geometry/V_all_geometry - V_keratin  
        
        vol_frac_water_cc = V_water_brick / V_brick_geometry * V_all_geometry
  
        return vol_frac_water_cc


    def compODEdydt(self, t, y, args=None):
        """ The wrapper function for computing the right hand side of ODEs
        """
        return comp.Comp.compODEdydt_diffu (self, t, y, args)
        
    def saveCoord(self, fn_x, fn_y) :
        comp.Comp.saveCoord(self, fn_x, fn_y, '.sc')
