# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 23:33:19 2018

@author: USER
"""

# Codes are free to use. Do whatever you want

from __future__ import absolute_import

"""Read raw data"""

####################### LIBRARY #############################

# exceptions library

# Python stdlib imports

# data processing library
import numpy as np

# pyrod library

####################### CONSTANT ############################

# constant 

####################### FUNCTIONS ###########################

def contral_vars(varl,
          key_list = ['roughness',
                      'vacancy',
                      'absorption',
                      'dis',
                      'dw',
                      'ocu',
                      'pos',
                      'scale',
                      'lattice',
                      'lattice_abc']):
    
    """
       contral independent var list 
    """
    
    # ubr_index--which kind if variable dose this var contral and their values
    ubr_index = {}
    
    for slab in varl:
        
        # ubr_slab_index-- the index of independent var list of one slab
        ubr_slab_index = {}
        
        # for substrate, there are two different modulation parametes
        #scale: bragg peaks intensity
        #beta: surface roughness
        
        if slab == 'substrate':
            ubr_slab_index['scale'] = 1e-4
            ubr_slab_index['beta'] = 0
        
        for key in key_list:
            
            slab_dict = varl[slab]
            
            # roughness-- contral variable: beta
            if key == 'roughness':
                
                ubr_slab_index[key] = 1
                
            # vacancy-- contral variable:alpha,layers,vacmi,vacma
            elif key == 'vacancy':
                
                ubr_slab_index[key] = np.array([1,1,1,1])

            # absorption-- contral variable: absorption
            elif key == 'absorption':
                
                ubr_slab_index[key] = 1

           # scale: contral the bragg peak intensity
            elif key == 'scale':
                
                ubr_slab_index[key] = 0
                
            elif key == 'lattice_abc':
                
                ubr_slab_index[key] = np.array([1,1,1])
                
            else:
                
                ubr_slab_index[key] = np.mat(np.ones(slab_dict[key].shape))
            
        ubr_index[slab] = ubr_slab_index
    
    # this parameter is used to scale the intensity between experiment data and 
    # fitting data
    ubr_index['intensity'] = 1
    
    return ubr_index
    
def initialize_contral_vars(ubr_index,slabs,keys):
    
    """    
       for multi-dimensional minimum, variables should be a list or array 
    """
    
    # to change the pandas matrix in ubr_index to list
    def tolist(mat):
        
        # get the matrix size
        c,r = mat.shape
        # target list
        tlist = []
        
        # Noted! this loop determine the path from matrix to list
        # recover list to matrix recall the revse path
        # line1: from left to right-line2:.......
        for ci in range(c):
            for ri in range(r):
                tlist.append(mat[ci,ri])
        
        return tlist
    
    # the parameters are rank as a list
    slab_ubr = []
    # the index of slab and keys, used to index the parameters and recover
    slab_ubr_index = []

    # loop of slabs
    for slab in slabs:
        
        slab_index = []
        ubr_slab = ubr_index[slab]
        
        # loop of keys in slab
        for key in keys:
            
            ubr_key = ubr_slab[key]
            
            # this parameters are single number
            if key in ['absorption','beta','roughness','scale']:
                
                slab_index.append([key,1])
                slab_ubr.append(ubr_key)
                
            # vacancy is a list
            elif key in ['vacancy','lattice_abc']:
                
                slab_index.append([key,len(ubr_key)])
                for i in ubr_key:
                    slab_ubr.append(float(i))
                    
            # these parameters are matrix
            else:
                slab_index.append([key,ubr_key.size])
                tlist = tolist(ubr_key)
                for i in tlist:
                    slab_ubr.append(i)
        
        slab_ubr_index.append([slab,slab_index])
        
    return slab_ubr, slab_ubr_index
            
def indicator_contral_vars(slab_ubr,slab_index,slab,key):
    
    """
       Note! slab is str,not a list    
       indicator value
    """
    # indicators for slab, keys and variables
    indicator_slab = 0
    indicator_keys = 0
    indicator_ubrs = [0,0]
    
    # used to break outer loop
    exit_flag = 0
    
    for slab_part in slab_index:
        
        if slab in slab_part:
            
            # if ture, slab is indicated, indicator for slab = 0
            indicator_slab = 0
            # keep on indicate the keys
            
            for key_part in slab_part[1]:
                
                # if ture, key is indicated, indicator for key = 0
                if key in key_part:
                    
                    # the postion of ubrs: first place = 0
                    #                      second place = first place + variable number -1(start from zeros)
                    indicator_keys  = 0
                    indicator_ubrs[1] = indicator_ubrs[0] + key_part[1] - 1
                    # since all the vars have been indicated, break out from the loop
                    # and set exit _flag  = 1 to break outer loop
                    exit_flag = 1
                    break
                
                else:
                    
                    # if key haven't been found, key indicator + 1, search for next key
                    indicator_keys = indicator_keys + 1                    
                    # the start postion of ubrs should change + variable number
                    indicator_ubrs[0] = indicator_ubrs[0] + key_part[1]
                
                #v if exit_flag  = 1, variable have been indcated, break out from the loop
                if exit_flag:
                    break
                    
        else:
            
            # slab haven't been found, slab indicator + 1, search for next slab
            indicator_slab = indicator_slab + 1
            # indicator keys count
            indicator_keys = indicator_keys + len(slab_part[1]) - 1
            # indicator ubrs count
            for key_part in slab_part[1]: 
                indicator_ubrs[0] = indicator_ubrs[0] + key_part[1]
                
        #v if exit_flag  = 1, variable have been indcated, break out from the loop
        if exit_flag:
            break
                
#    indicated_ubrs = slab_ubr[indicator_ubrs[0]:indicator_ubrs[1]+1]
    
    return [indicator_ubrs[0],indicator_ubrs[1]+1]

# if contral variable is changed, update the variable value in ubr_index
def refresh_index(ubr_index,slab_ubr,slab_index):
    
    # locate slab
    for slabi in slab_index:
        # slab name
        slab_k = slabi[0]
        # keys in slab
        for keyi in slabi[1]:
            # parameters name
            key = keyi[0]
            # ir-- indicator ubr
            ir = indicator_contral_vars(slab_ubr,slab_index,slab_k,key)
            # iedr-- inidcated ubr
            iedr = slab_ubr[ir[0]:ir[1]]
            # ogr -- original ubr
            ogr = ubr_index[slab_k][key]
            
            # single vaule parameters
            if key in ['absorption',
                       'roughness',
                       'scale']:
                ubr_index[slab_k][key] = iedr
            # list vaule parameters                
            elif key in ['lattice_abc',
                         'vacancy']:                
                ubr_index[slab_k][key] = iedr
            # matrix vaule paramters
            elif key in ['dis','dw','lattice','ocu','pos']:
                iedr_m = np.reshape(iedr,ogr.shape)
                ubr_index[slab_k][key] = iedr_m
                
    return ubr_index

def bragg_mask(iq,w,weight,mode = 'yang',limit = 0):
    
    """    
       set bargg_mask for weight fitting
    """
    # calculate bragg peaks indice 
    qs = int(iq[ 0]) + 1        
    qe = int(iq[-1])
    
    # bg--bragg peaks
    bg = range(qs, qe + 1)
    
    # pb--bg postions
    pb = []
    # bp-- bragg peak
    for bp in bg:
        # indicate bragg peak postion
        pb.append(np.argmin(abs(iq - bp)))
        
    if pb[0] < w or (len(iq) - pb[-1]) < w:
        print('Warning!mask width is too large!!')
        
    # set empty mask
    yang_mask = np.zeros(len(iq))
    yin_mask = np.ones(len(iq))
    bragg_mask = np.zeros(len(iq))
    
    # cut window in empty mask
    for p in pb:
        # cut yang window
        if mode == 'yang':
            yang_mask[p-w:p+w] = 1
            yang_mask[0:int(limit)] = 1
            bragg_mask = yang_mask
        # cut yin window
        elif mode == 'yin':
            yin_mask[p-w:p+w] = 0
            yin_mask[0:int(limit)] = 0
            bragg_mask = yin_mask
            
    return bragg_mask*weight

######################## CLASSS #############################

    