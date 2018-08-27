# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 18:01:27 2018

@author: USER
"""

# Codes are free to use. Do whatever you want

from __future__ import absolute_import

"""Read raw data"""

####################### LIBRARY #############################

# exceptions library
from exceptions import (Illegal_Filename_Exception, 
                        Data_Format_Exception)

# Python stdlib imports
import datetime
import re
import os

# data processing library
import numpy as np
import pandas as pd

# pyrod library

####################### CONSTANT ############################

# constant 
#regex of raw_data, parameters, sheet names

RAW_DATA_DAT  = re.compile('^\d\dL\.dat$')
RAW_DATA_XLSX = re.compile('^\d\dL\.xlsx$')

####################### FUNCTIONS ###########################


    
######################## CLASSS #############################

class initialization_rhkl(object):
    
    def __init__(self,
                 raw_data = '00L.xlsx'):
        
        self.raw_data = raw_data
        self.path = os.path.abspath(os.path.dirname(raw_data)) + raw_data
                    
        # check raw_data
        self.check_raw_data_dat = 1
        self.check_raw_data_xlsx = 1
        
    # check file names
    # raw_data and parameters
    def _check_filename(self):
    
        """legal names: xxL.xlsx or xxL.dat; parameters_xxx--.xlsx"""
        
        # check raw data, .dat or xlsx
        self.check_raw_data_dat = RAW_DATA_DAT.match(self.raw_data)
        self.check_raw_data_xlsx = RAW_DATA_XLSX.match(self.raw_data)
        
        if not (self.check_raw_data_dat or self.check_raw_data_xlsx):
            error = 'raw data name is illegal.appropriate file name: xxL.dat or xxL.xlsx'
            raise Illegal_Filename_Exception(error)
            

    def _read_rhkl(self, density = 100):
        
        """Read origin data.default reconstruct data density is 100"""
        
        try:
            # if data is origin data form from APS .dat
            if 'dat' in self.raw_data:
                
                # open od-origin data
                raw_data = open(self.path)
                
                qz = []
                ie = []
                
                # read in line form
                for line in raw_data:
                    
                    lst = line.split()
                    
                    # float data
                    qz.append(float(lst[0]))
                    ie.append(float(lst[1]))
            
            # if data is xlsx, maybe modulated
            elif 'xlsx' in self.raw_data:
                
                # read excel data as matrix
                rd = pd.read_excel(self.path).as_matrix()
                
                qz = rd[:,0].tolist()
                ie = rd[:,1].tolist()
                
        except Data_Format_Exception:
            print('Data format of raw data is illegal')
            
        # interpolant data density. ensure the integrality of bragg peak
        # qs--q start
        qs = qz[ 0]
        # qe--q end
        qe = qz[-1]
        
        # interpolant data with default intensity 100     
        iq0 = 0.0
        iqs = round(qs*100)/100
        iqe = round(qe*100)/100
        
        iq= np.linspace(iq0, 
                        iqe,
                        (iqe-iq0)/0.01+1)
        
        intensity = np.interp(iq,qz,ie)
        # the signal between iq0 and iqs is not detected
        intensity[0: int(iqs/0.01)] = 0
        
        return iq, intensity