# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:03:37 2018

@author: USER 
"""


# Codes are free to use. Do whatever you want

from __future__ import absolute_import 

"""Read raw data"""

####################### LIBRARY #############################

# exceptions library
from exceptions import (Excel_Load_Exception,
                        Data_Load_Exception,
                        Data_Format_Exception)

# Python stdlib imports
import os
import pickle

# data processing library
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# pyrod library

from read.read_parameters import initialization_parameters
from read.read_raw_data import initialization_rhkl

####################### CONSTANT ############################

# constant 

R0 = 2.82e-5 # Thompson scattering length, the radi of electron in Angs
NA = 6.022e23 # Avogadros number corresponding to g/cm3 e23

# atomic_mass.xlsx and atomic_scattering_factor.pickle

ATOMIC_MASS = os.path.abspath(os.path.dirname('atomic_mass.xlsx')) +\
                     '/base/atomic_mass.xlsx'
ASF = os.path.abspath(os.path.dirname('atomic_scattering_factor.pickle')) +\
                     '\\base\\atomic_scattering_factor.pickle'

try:
    ATOM_TABLE = pd.read_excel(ATOMIC_MASS,
                               sheet_name = 'Sheet1',
                               index_col = 0)
    ATOM_LIST = ATOM_TABLE.index.values.tolist()
    
    f = open(ASF,'rb')
    ASF_TABLE = pickle.load(f)
    f.close()
    
except Excel_Load_Exception:
    print('load atomic form factor data base load fail')

# load ctr optimised data

CTR_PATH = os.path.abspath(os.path.dirname('ctr_optimised_result.pickle')) +\
                     '/data/ctr_optimised_result.pickle'

ctr_optimised_result = {}

try:
    with open(CTR_PATH, 'rb') as f:
        unpickler = pickle.Unpickler(f)
        ctr_optimised_result = unpickler.load()
    
except Data_Load_Exception:
    print('Please optimise and save the CTR data first!')

# save xrr optimised data path

XRR_PATH = os.path.abspath(os.path.dirname('xrr_optimised_result.pickle')) +\
                     '/data/xrr_optimised_result.pickle'
                     
####################### FUNCTIONS ###########################

# elements = [A,B]
# result density - g/cm3
def density(elements,lattice):
    
#elements = ['Sr','Ti']
#lattice = [3.905,3.905,3.905,np.pi/2,np.pi/2,np.pi/2]
    
    """theroy density of ABO3 pervoskite oxide"""
    
    try:
        
        A_mass = ATOM_TABLE.at[elements[0],'mass']
        B_mass = ATOM_TABLE.at[elements[1],'mass']
        O_mass = ATOM_TABLE.at['O','mass']
        
        a = lattice[0]
        b = lattice[1]
        c = lattice[2]
        
        alpha = lattice[3]*np.pi/180
        beta  = lattice[4]*np.pi/180
        gamma = lattice[5]*np.pi/180
        
    except Data_Format_Exception:
        print('Parameters error!')
        
    lattice_volume = a*b*c*np.sqrt(1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)-\
                                   np.cos(alpha)**2 - \
                                   np.cos(beta)**2 - \
                                   np.cos(gamma)**2)
    A_unitcell = A_mass
    B_unitcell = B_mass
    O_unitcell = 3*O_mass
    
    na = NA/1e23
    
    density_A = A_unitcell/(lattice_volume*na*0.1)
    density_B = B_unitcell/(lattice_volume*na*0.1)
    density_O = O_unitcell/(lattice_volume*na*0.1)
    
    density_pervoskite = density_A + density_B + density_O
        
    return density_pervoskite, density_A, density_B, density_O
    
# wave length should be Angs
# energy is the x-ray energy kev
# elements = [A,B]
def imag_f(elements, energy):
    
    """imag part of x-ray scattering factor for a certin element"""
    
    A_energys = ASF_TABLE[elements[0]]['E(kev)']
    B_energys = ASF_TABLE[elements[1]]['E(kev)']
    O_energys = ASF_TABLE['O']['E(kev)']
    
    A_imagf = ASF_TABLE[elements[0]]['imag_f']
    B_imagf = ASF_TABLE[elements[1]]['imag_f']
    O_imagf = ASF_TABLE['O']['imag_f']
    
    def f2(energy, energys, imag_f):
        
        atom_f2 = np.interp(energy, energys, imag_f)
        
        return atom_f2
    
    A_f2 = f2(energy, A_energys, A_imagf)
    B_f2 = f2(energy, B_energys, B_imagf)
    O_f2 = f2(energy, O_energys, O_imagf)
    
    return A_f2, B_f2, O_f2

# energy is x-ray energy kev
# elements = [A,B]
# atttenuation coefficiment unit is Angs
def attenuation_coefficient(elements, densities, f2, energy):
    
    wave_length = 12.38/energy # unit Angs
    na = NA/1e23
    
    A_ac = (densities[1]*na/ATOM_TABLE.at[elements[0],'mass'])*2*R0*wave_length*f2[0]
    B_ac = (densities[2]*na/ATOM_TABLE.at[elements[1],'mass'])*2*R0*wave_length*f2[1]
    O_ac = (densities[3]*na/ATOM_TABLE.at['O','mass'])*2*R0*wave_length*f2[2]
    
    ac = A_ac + B_ac + O_ac
    
    return ac

######################## CLASSS #############################
        
class XRR(initialization_parameters, initialization_rhkl):
    
    def __init__(self, parameters, experiment_data, energy):
        
        initialization_parameters.__init__(self, parameters)
        initialization_parameters._var_list(self)
        initialization_parameters._var_table(self)
        initialization_rhkl.__init__(self, experiment_data)
            
        self.energy = energy  # unit - kev
        self.wave_length = 12.38/energy # unit-Angs
        self.coefficient = {}
        
        self.iq,self.i = self._read_rhkl(self)
        self.lattice_c = self.var_list['substrate']['lattice'].at['constant','c']
        self.q = self.iq*2*np.pi/self.lattice_c
        
        for key in self.var_list:
            
            element = [self.var_list[key]['ions'].columns[3], self.var_list[key]['ions'].columns[0]]
            lattice = self.var_list[key]['lattice'].as_matrix().tolist()[0]
            
            densities = density(element, lattice)
            f2 = imag_f(element, energy)
            ac = attenuation_coefficient(element, densities, f2, energy)
            
            self.coefficient[key] = [element,lattice,densities,ac]
            
        # stack list
        
        self.scattering_factor = [] # density*R0 + 1j*attenuation_coefficient
        self.d_space = [] # stack thickness
        self.roughness = [] # roughness at each layer
        
        for layeri in range(len(self.var_table['posz_list'])):
            
            self.d_space.append(self.var_table['posz_table'][0,0]*self.lattice_c*0.98)
            self.roughness.append(0)
            
            slab = self.var_table['slab_list'][layeri]
            dens = self.coefficient[slab][2][0]
            ac = self.coefficient[slab][3]
            
            self.scattering_factor.append(dens*R0+1j*ac)
            
        self.scattering_factor.append(1e-21) # the scattering factor of vacuum
        self.roughness.append(0)
        self.roughness[0] = 0
        self.thickness = np.sum(self.d_space)
        
        # rt - ratio
        self.rt = 1
        
    # re initilization the factors-d_space, scattering_factor, roughness, thichness
    def re(self):
        
        # stack list
        
        self.scattering_factor = [] # density*R0 + 1j*attenuation_coefficient
        self.d_space = [] # stack thickness
        self.roughness = [] # roughness at each layer
        
        for layeri in range(len(self.var_table['posz_list'])):
            
            self.d_space.append(self.var_table['posz_table'][0,0]*self.lattice_c*0.98)
            self.roughness.append(0)
            
            slab = self.var_table['slab_list'][layeri]
            dens = self.coefficient[slab][2][0]
            ac = self.coefficient[slab][3]
            
            self.scattering_factor.append(dens*R0+1j*ac)
            
        self.scattering_factor.append(1e-21) # the scattering factor of vacuum
        self.roughness.append(0)
        self.roughness[0] = 0
        self.thickness = np.sum(self.d_space)

    def disp(self):
        
        # list all the self parameters
        
        print("--------------------------------------------")
        print("densities:\n")
        for key in self.coefficient:
            print("    " + key + ": %s\n" %self.coefficient[key][2][0])
            
        print("scattering factor: %s - %s\n" % (self.scattering_factor[ 0], 
                                                self.scattering_factor[-1]))
        print("roughness: %s - %s\n" % (self.roughness[ 0],
                                        self.roughness[-1]))
        print("d_space: %s - %s Angs\n" % (self.d_space[ 0],
                                           self.d_space[-1]))
        print("thichness: %s Angs\n" % np.sum(self.d_space))
        print('--------------------------------------------')
        
    # homogeneous slab reflection
    def homo_slab(self):
        
        slabs = list(self.var_list.keys())
        slabs.remove('substrate')
        
        # check if the slab is homogeneous
        if len(slabs) != 1:
            raise Data_Format_Exception('Not homogeneous slab!')
        
        dens = self.coefficient[slabs[0]][2][0]
        thickness = np.sum(self.d_space)
            
        reflect_inten = -1j*(4*np.pi*dens*R0*thickness/self.q)*\
                            (np.sin(self.q*thickness/2)/(self.q*thickness/2))*\
                            np.exp(1j*self.q*thickness/2)
                                
        return self.iq, reflect_inten
    
    # Parratt reflectivities 
    # Maybe a little difficult while optimizing parameter, parratt model is accurate
    # Recomed method
    def parratt(self):
        
        """"Elements of Modern X-ray Physics" by Jens Als-Nielsen and Des McMorrow,
        Calculates: Parratt reflectivity of a multilayer"""
        
        k = 2*np.pi/self.wave_length # diffraction vector#
        layer_num = len(self.scattering_factor) # layer number, vacuum is added
        
        #----- Calculate refractive index n of each layer
        delta = self.wave_length**2*np.real(self.scattering_factor)/(2*np.pi)
        beta  = self.wave_length*np.imag(self.scattering_factor)/(4*np.pi)
        # relfractive index
#        nu = 1 - delta + 1j*beta
        
        #----- Wavevector transfer in each layer
        trans_vector = np.zeros([layer_num+1, len(self.q)], dtype = np.complex)
        trans_vector[0,:] = self.q
        
        for i in range(layer_num):
            trans_vector[i+1,:] = np.sqrt(self.q**2 - 8*k**2*delta[i] + 1j*8*k**2*beta[i])
            
        #----- Reflection coefficients (no multiple scattering)
        reflect_coe = np.zeros([layer_num, len(self.q)], dtype = np.complex)
        
        for i in range(layer_num):
            reflect_coe[i,:] = ((trans_vector[i,:] - trans_vector[i+1,:])/\
                                (trans_vector[i,:] + trans_vector[i+1,:]))*\
                                np.exp(-0.5*trans_vector[i,:]*trans_vector[i+1,:]*self.roughness[i])
                                
        #----- Reflectivity from first layer
        reflectivity = np.zeros([layer_num-1, len(self.q)], dtype = np.complex)
        
        phase1 = np.exp(1j*trans_vector[layer_num-1]*self.d_space[layer_num-2])
        
        if layer_num > 1:
            reflectivity[0,:] = (reflect_coe[layer_num-2,:] + \
                                 reflect_coe[layer_num-1,:]*phase1)/\
                                (1 + reflect_coe[layer_num-2,:]*\
                                 reflect_coe[layer_num-1,:]*phase1)
        if layer_num > 2:
            for i in range(1,layer_num-1):
                
                phasei = np.exp(1j*trans_vector[layer_num-i-1,:]*self.d_space[layer_num-i-2])
                
                reflectivity[i,:] = (reflect_coe[layer_num-i-2,:] +\
                                     reflectivity[i-1,:]*phasei)/\
                                    (1 + reflect_coe[layer_num-i-2,:]*\
                                     reflectivity[i-1,:]*phasei)
        
        #------ Intensity reflectivity
        
        # should be reminded here! The data is not squared! To keep uniform with ctr fitting
        
        if layer_num == 1:
            reflect_inten = reflect_coe[0,:]
        else:
            reflect_inten = reflectivity[-1,:]
            
        return self.iq, reflect_inten
    
#    def p_xrr(q, inten):
#        
#        plt.plot(q, np.log(inten))

class xr(XRR):
    
    # the relax of surface, secface  and interface
    # parratt method only
    
    def xr_relax(self, interface = 0.9, secface = 1, surface = 1.03, rt = 1):
        
        XRR.re(self)
        
        self.rt = rt
        self.d_space[ 0] = interface
        self.d_space[-2] = secface
        self.d_space[-1] = surface
        
        q, inten = XRR.parratt(self)
        
        plt.plot(q, np.log(abs(self.rt*inten + ctr_optimised_result['substrate_ctr'])))
        plt.plot(q, np.log(abs(ctr_optimised_result['shkl'])))
        
        return q, inten
    
    # modulate the roughness at interface and surface
    # parratt method only
    def xr_roughness(self, interface = 0.1, surface = 0.1, rt = 1):
        
        XRR.re(self)
        
        self.rt = rt
        self.roughness[ 0] = interface
        self.roughness[-1] = surface
        
        q, inten = XRR.parratt(self)
        
        plt.plot(q, np.log(abs(rt*inten + ctr_optimised_result['substrate_ctr'])))
        plt.plot(q, np.log(abs(ctr_optimised_result['shkl'])))
        
        return q, inten
    
    # thickness modulation
    # homo slab method or parratt method
    def xr_thickness(self, miu = 0, mode = 'parratt', rt = 1):
        
        XRR.re(self)
        
        self.rt = rt
        self.d_space = (np.array(self.d_space) - miu).tolist()
        
        if mode == 'homo_slab':
            q, inten = XRR.homo_slab(self)
        elif mode == 'parratt':
            q, inten = XRR.parratt(self)
            
        plt.plot(q, np.log(abs(rt*inten + ctr_optimised_result['substrate_ctr'])))
        plt.plot(q, np.log(abs(ctr_optimised_result['shkl'])))
            
        return q, inten
    
    # export optimised xrr data to /data
    def save(self, mode = 'parratt'):
        
        xrr_optimised_result = {}
        
        if mode == 'homo_slab':
            q, inten = XRR.homo_slab(self)
        elif mode == 'parratt':
            q, inten = XRR.parratt(self)
            
        xrr_optimised_result['q'] = q
        xrr_optimised_result['inten'] = inten
        xrr_optimised_result['rt'] = self.rt
        
        f = open(XRR_PATH,'wb')
        pickle.dump(xrr_optimised_result, f)
        f.close()
        
        