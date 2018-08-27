# -*- coding: utf-8 -*-
"""
Copyright (c) 2010-2014 pyrod

@author: Han Xu

@mail: xuhan@mail.ustc.edu.cn

"""
# Codes are free to use. Do whatever you want

from __future__ import absolute_import

"""Read orginaized parameters"""

####################### LIBRARY #############################

# exceptions library
from exceptions import (Illegal_Filename_Exception, 
                        Illegal_Sheetname_Exception,
                        Excel_Load_Exception,
                        Data_Base_Exception)

# Python stdlib imports
import datetime
import re
import os

# data processing library
import numpy as np
import pandas as pd
import openpyxl as ox
import xlrd

# pyrod library

####################### CONSTANT ############################

# constant 
#regex of parameters, sheet names

PARAMETERS    = re.compile('^parameters_?\w*\.xlsx$')
SHEET_SLAB    = re.compile('^slab\d$')

# location of parameters in substrate,slabs
# [header,usecols,index_col,skip_footer]

PARAMETERS_LOCATE = {'properties':  [0,1,0,19],
                     'ions':        [6,5,0,16],
                     'lattice':     [4,6,0,18],
                     'ocu':         [8,5,0,14],
                     'dw':          [10,5,0,10],
                     'dis':         [14,5,0,6],
                     'nel':         [18,5,0,4],
                     'pos':         [20,5,0,0]} 

# atomic_form_factor.xlsx

ATOMIC_FORM_FACTOR = os.path.abspath(os.path.dirname('atomic_form_facotr.xlsx')) +\
                     '/base/atomic_form_factor.xlsx'

try:
    IONS_TABLE = pd.read_excel(ATOMIC_FORM_FACTOR,
                               sheet_name = 'atomic_form_factor',
                               index_col = 0)
    IONS_LIST = IONS_TABLE.index.values.tolist()
    
except Excel_Load_Exception:
    print('load atomic form factor data base load fail')

####################### FUNCTIONS ###########################

# load the parameters from parameter excel

def _load_parameters(file,parameter,SHEET):
    
    # paramters from PARAMETERS_LOCATE
    
    PARS = PARAMETERS_LOCATE[parameter]
    
    HEADER      = PARS[0]
    USECOLS     = PARS[1]
    INDEX_COL   = PARS[2]
    SKIP_FOOTER = PARS[3]
    
    file_path = os.path.abspath(os.path.dirname(file)) +\
                '/data/' + file    
    try:
        variable = pd.read_excel(file_path,
                                 sheet_name = SHEET,
                                 header = HEADER,
                                 usecols = USECOLS,
                                 index_col = INDEX_COL,
                                 skip_footer = SKIP_FOOTER)
    except:
        
        error = 'Error while loading parameters using pandas'
        raise Excel_Load_Exception(error)
    
    return variable

# check if ions are included in atomic form factor data base

def check_ions(ions):
    
    # input one ion, str
    if isinstance(ions, str):
        
        # ion is included
        if ions in IONS_LIST:
            pass
        # not included, raise error
        else:
            error = 'ion ' + ions + ' is not included'
            raise Data_Base_Exception(error)
    
    # input seveal ions, list
    elif isinstance(ions, list):
        
        ions_not_included = []
        flag = 1
        
        for ion in ions:
            
            #ion is included
            if ion in IONS_LIST:
                pass
            #ion is included, record
            else:
                ions_not_included.append(ion)
                flag = 0
        # not included, raise error
        if not flag:
            error = 'ions: ' + str(ions_not_included) + ' not included'
            raise Data_Base_Exception(error)

# append atomic form factor of ions to atomic_form_factor.xlsx
# a--atomic form factor
            
def append_ions(ion_a):
    
    # openpyxl read atomic_form_facotr
    # loaded xlsx
    try:
        pxlsx = ox.load_workbook(ATOMIC_FORM_FACTOR)
    except Excel_Load_Exception:
        print('openpyxl load excel error')
    # loaded sheet
    try:
        psheet = pxlsx['atomic_form_factor']
    except Illegal_Sheetname_Exception:
        print('wrong sheet name')
        
    # new appended ions location
    location = 'A' + str(len(IONS_LIST)+1) + ':' + 'J' + str(len(IONS_LIST)+1)
    pcell = psheet[location]
    
    # load new atomic form factor
    
    index = 0
    
    for row in pcell:
        for element in row:
            
            element.value = ion_a[index]
            index += 1
    
    # reload atomic_form_factor.xlsx
    try:
        ions_table = pd.read_excel(ATOMIC_FORM_FACTOR,
                               sheet_name = 'atomic_form_factor',
                               index_col = 0).index.values.tolist()
    except Excel_Load_Exception:
        print('reload atomic form factor data base load fail')
    
    return ions_table

######################## CLASSS #############################

# experiment raw data and parameters initialization
class initialization_parameters(object):
    
    def __init__(self,
                 parameters = 'parameters.xlsx'):

        self.parameters = parameters
        self.path = os.path.abspath(os.path.dirname(parameters)) + parameters
                    
        # check files
        
        self.check_parameters = 1
        
        # check sheets
        
        self.check_sheet_fir = 1
        self.check_sheet_las = 1
        self.check_sheet_slab = 1
        self.check_sheet_rank = 1
        
        self.var_table = {}
        self.var_list  = {}
        
        # layers_sum-- slab1 layer number, slab1+slab2 layer number, 1+2+3 layer number .....
        self.layers_sum = []
        # layers_n-- 1 layer number, 2 layer number ,3.....
        self.layers_n   = []
        # atoms_sum-- 1 atoms number, 2 atoms numbers......
        self.atoms_sum  = []
        # atoms_max-- max atoms number of all slabs
        self.atoms_max  = 0
        # layers_max-- max layers number of all slabs
        self.layers_max = 0

    def _check_parameters(self):
        
        # check parameters, xlsx
        self.check_parameters = PARAMETERS.match(self.parameters)
        
        if not self.check_parameters:
            error = 'parameters name is illegal.appropriate file name: parameters_xxx--.xlsx'
            raise Illegal_Filename_Exception(error)

    # check sheet names
    def _check_sheetname(self):
        
        """legal names: substrate slab1 slab2 slab3 .... factors"""
        
        # read sheet names
        sheet_names = xlrd.open_workbook(self.path, 
                                         on_demand=True).sheet_names()
        
        self.check_sheet_fir = sheet_names[ 0] == 'substrate'
        self.check_sheet_las = sheet_names[-1] == 'factors'
        slab_names = sheet_names[1:-1]
        
        # check slab names
        self.check_sheet_slab = 1
        self.check_sheet_rank = 1
        sheet_rank = []
        
        for slab in slab_names:
            self.check_sheet_slab = self.check_sheet_slab and SHEET_SLAB.match(slab)
            try:
                sheet_rank.append(int(slab.split('b')[1]))
            except:
                self.check_sheet_slab = 0
    
        self.check_sheet_rank = sheet_rank == np.linspace(1,len(sheet_rank),
                                                          len(sheet_rank),
                                                          dtype = int).tolist()    
        if not self.check_sheet_fir:
            error = 'The first sheet should be "substrate"'
            raise Illegal_Sheetname_Exception(error)
            
        elif not self.check_sheet_slab:
            error = 'appropriate slab names'
            raise Illegal_Sheetname_Exception(error)
            
        elif not self.check_sheet_rank:
            error = 'Slab names are not in right order'
            raise Illegal_Sheetname_Exception(error)
            
        elif not self.check_sheet_las:
            error = 'The las sheet should be "factors"'
            raise Illegal_Sheetname_Exception(error)

    # read all the parameters
    def _var_list(self):

        """ Read model parameters. construct model parameters list and table"""
        
        pars       = {}
        
        xls = xlrd.open_workbook(self.path, on_demand=True)

        for sheet_n in xls.sheet_names()[0:-1]:
            
            for key in PARAMETERS_LOCATE.keys():
                
                if key == 'properties':
                    
                    # properties-- atoms number , layers number and absorption
                    properties = _load_parameters(self.parameters,'properties',sheet_n)
                    
                    # atoms_num-- atoms number in one crysatl lattice
                    atoms_num  = properties.at['atoms_num','c1']
                    # absorption of one layers
                    absorption = properties.at['absorption','c1']
                    # layers_num-- layers number in one slab
                    layers_num = properties.at['layers_num','c1']
                    
                    self.layers_max = self.layers_max + layers_num
                    self.layers_sum.append(self.layers_max)
                    self.layers_n.append(layers_num)
                    self.atoms_sum.append(atoms_num)
                
                    if atoms_num > self.atoms_max:
                        self.atoms_max = atoms_num
                    
                    pars['atoms_num']  = atoms_num
                    pars['absorption'] = absorption
                    pars['layers_num'] = layers_num
                    
                else:
                    
                    pars[key] = _load_parameters(self.parameters,key,sheet_n)

                    
                pars['roughness'] = 1
            
            self.var_list[sheet_n] = pars
            
            pars = {}
                
        return self.var_list
#    
    # construct parameters tables
    def _var_table(self):
        
        
        self.var_table = {}
        # slab_index-- which slab    
        slab_index = 0
        # tot_layer-- total layers till the layer under calculation
        tot_layer = 0
#    
        # talbe--matrix total layers x max atoms
        x, y = int(self.layers_max), int(self.atoms_max)       
        
        posx_table = np.mat(np.zeros([x, y]))
        posy_table = np.mat(np.zeros([x, y]))
        posz_table = np.mat(np.zeros([x, y]))
        ocup_table = np.mat(np.zeros([x, y]))
        elec_table = np.mat(np.zeros([x, y]))
        ions_table = []
        dw_table   = []
        slab_list  = []
        
        # initalize lattice c modulation list
        posz_list = np.ones(x)
        posx_list = np.ones(x)
        posy_list = np.ones(x)
        
        # contral_table is used to multiply the contral variables
        contral_table = np.mat(np.zeros([x, y]))
   
        # first, we should locate and remove the zero value from layers_sum
        # for zero values equals to the inexistence of this slab
        layers_sum_exist = []
        layers_index_exist = []
        exist_slab = np.nonzero(np.array(self.layers_sum))[0]
        
        for s in exist_slab:
            layers_sum_exist.append(self.layers_sum[int(s)])
            layers_index_exist.append(int(s))
    
        for layer_index in range(x):
            
            # a-- determine which layer this layer_index close with
            # [5,6,10]-8 = [-3,-2,2]
            # "numpy.argmin return the first minimux index"
            # if result>0 slab index is a
            # elif result<0 slab index is a+1
            a = np.argmin(abs(np.array(layers_sum_exist) - layer_index -1))
            
            if (np.array(layers_sum_exist)-layer_index)[a] > 0:
                slab_index = layers_index_exist[a]
            elif (np.array(self.layers_sum)-layer_index)[a] <= 0:
                slab_index = layers_index_exist[a + 1]
            else:
                print('There must be something wrong! layer index should not bigger than total layer number!')
            
            # recognise which slab
            # index which slab
            xls = xlrd.open_workbook(self.path, 
                                     on_demand=True)
            slab_list.append(xls.sheet_names()[slab_index])
            par = self.var_list[xls.sheet_names()[slab_index]]
            
            # normalize the lattice matrix and cosntant by substrate lattice 
            # reduce_lattice = slab_lattice/substrate_lattice
            substrate_lattice = self.var_list['substrate']['lattice'].as_matrix()[0,0:3]
            slab_lattice = par['lattice'].as_matrix()[0,0:3]
            reduce_lattice = slab_lattice/substrate_lattice
            
            # postion x,y and z are constructed by reduce_lattice
            posx_table[layer_index,] = par['pos'].loc['pos_x'].as_matrix()*reduce_lattice[0]
            posy_table[layer_index,] = par['pos'].loc['pos_y'].as_matrix()*reduce_lattice[1]
            
            ocup_table[layer_index,] = par['ocu'].as_matrix()
            elec_table[layer_index,] = par['nel'].as_matrix()
            ions_table.append(par['ions'].as_matrix().tolist()[0])
            dw_table.append(par['dw'])
            
            # tot_layer-- the layer postion in the current slab
            tot_layer = tot_layer + reduce_lattice[2]
            # the atom postion in the current layer
            posz_table[layer_index] = tot_layer + \
                                      par['pos'].loc['pos_z'].as_matrix()*reduce_lattice[2] - \
                                      self.layers_sum[0]*1

        self.var_table = {'posx_table':   posx_table,
                          'posy_table':   posy_table,
                          'posz_table':   posz_table,
                          'posz_list':    posz_list,
                          'posx_list':    posx_list,
                          'posy_list':    posy_list,
                          'ocu_table':    ocup_table,
                          'nel_table':    elec_table,
                          'ion_table':    ions_table,
                          'dw_table':     dw_table,
                          'slab_list':    slab_list,
                          'contral_table':contral_table}
        
        return self.var_table
    
    