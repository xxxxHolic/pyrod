# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 20:14:16 2018

@author: USER
"""

# Codes are free to use. Do whatever you want

from __future__ import absolute_import

"""Read raw data"""

####################### LIBRARY #############################

# exceptions library
from exceptions import (Data_Format_Exception,
                        Data_Match_Exception)

# Python stdlib imports
import datetime

# data processing library
import numpy as np

# pyrod library

####################### CONSTANT ############################

# constant 

####################### FUNCTIONS ###########################



######################## CLASSS #############################


# film property--roughenss, surface roughness and interface roughenss

class roughness(object):
    
    def __init__(self, iq):
        
        self.q = iq
    
    # r_beta -- robbinson roughness factor
    def robinson_roughness(self, r_beta):
        
        """Robinson roughness model: beta is the film fraction in first layer.
           In surface layer n, the fraction is beta^n. 
           This fuction is for cubic crystal, which is common for pervoskite
           
           [1]. Robinson, Ian K. "Crystal truncation rods and surface roughness." 
                Physical Review B 33.6 (1986): 3830."""
                
        if r_beta < 0:
            raise Data_Format_Exception('robinsion_roughness factor should be positive')
            
        r_factor = (1-r_beta)/((1-r_beta)**2 + 4*r_beta*(np.sin(np.pi*self.q)**2))**0.5

        return r_factor
    
    def interface_roughness(self, i_beta):
        
        """interface roughness model, refinement for small angle reflection
           exp(-q2*r2)"""
           
        if i_beta < 0:
            raise Data_Format_Exception('interface_roughness factor should be positive')
            
        i_factor = np.exp(-1*(i_beta*self.q))
        
        return i_factor
    
# film vacancy--gradient vacancy, layer surface 

class vacancy(object):
    
    def __init__(self, iq, var_list, var_table):
        
        self.q = iq
        
        self.l = var_list
        self.t = var_table
        
    def gradient_vacancy(self, alpha, va, vi):
        
        """
           Normally, A and O is easy to exist vancany in ABO3, and exponential function is 
           good fitting in XRD fitting. However, some materials have special vancany mdoel
           For example: the oxygen vancany distribtion of SCoO3 and SrFeO3, oxygen vancany 
           is very complex. Even for STO, linear oxygen vancany clustering is very common 
           in thermal threatment.
           ...............................
           
           [1].Druce, J., et al. "Surface termination and subsurface restructuring of 
               perovskite-based solid oxide electrode materials." Energy & Environmental 
               Science 7.11 (2014): 3593-3599.
        """
        
        # check and put vacancy at interface and surface in order
        if va < 0 or vi < 0 or alpha < 0:
            raise Data_Format_Exception('interface_roughness factor should be positive')
            
        vs = np.max([va, vi])
        ve = np.min([va, vi])
        
        # layers
        layer = range(len(self.t['posz_list']))
        
        # factors
        
        a = (vs - ve)/(np.exp(-alpha*(layer-1)) - 1)
        b = ve - a
        
        gradient_v = a + b*np.exp(-alpha*np.array(layer))
        
        return gradient_v
    
    def surface_vacancy(self, ions, vs, layers = [-1]):
        
        """vacancy for specified layers and ions. Normally, A site ion and O2- is vacancy"""
        
        if not len(ions) == len(vs):
            raise Data_Match_Exception('ions number should match vacancy number!')
            
        ocu_table = self.t['ocu_table']
        ion_table = self.t['ion_table']
        
        ion_index = 0
        
        for ion in ions:
            for layer in layers:
                
                layer_ion  = ion_table[-1*(1 + layer)]
                ion_locate = layer_ion.index(ion)
                ocu_table[-1*(1 + layer), ion_locate] = 1-vs[ion_index]
                
            ion_index += 1
            
        return ocu_table
        
# strain relex for epitaxial films. 50uc, normally.
# interface relax and surface relax

class relax(object):
    
    def __init__(self, posz_list):
        
        self.pz = posz_list
        
    def stain_relax(self, cl, bc, sc):
        
        # sr -- strain relax
        # sc -- strained lattice c
        # cl -- critalc layaers number
        # bc -- bulk lattice c
        
        ln = len(self.pz)
        sr = self.pz*((sc-bc)*np.exp(-np.array(range(ln)/cl) + bc))
        
        return sr
    
    # layer relax include interface relax and surface relax
    # surface relax is common and well studied
    # interface relax and second surface layer relax is also included
    
    def layer_relax(self, rv, layers):
        
        if not len(rv) == len(layers):
            raise Data_Match_Exception('relax value number should match layer number')
            
        lr = self.pz
        
        # input relax value in posz_list
        for layer in layers:
            lr[layer] = rv[layer]
            
        return lr
    
    