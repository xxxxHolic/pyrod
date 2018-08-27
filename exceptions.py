# -*- coding: utf-8 -*-
"""
Copyright (c) 2010-2014 pyrod

@author: Han Xu

@mail: xuhan@mail.ustc.edu.cn

"""
# Codes are free to use

from __future__ import absolute_import

"""Definitions for pyrod shared exception classes."""

class Illegal_Filename_Exception(Exception):
    """Error for Illegal raw data and parameters names
       raw data name: hkl.dat or hkl.xlsx, hk is variable, l is str
       parameters name: parameters_xxxx.xlsx"""
       
class Illegal_Sheetname_Exception(Exception):
    """Error for Illegal sheet names
       sheet name: substrate,slab1,slab2,.. slabxx,factors"""

class Excel_Load_Exception(Exception):
    """Error for loading data"""

class Data_Base_Exception(Exception):
    """inputs are not stored in data base"""

class Data_Format_Exception(Exception):
    """Error for any data format inconsistencies."""

class Data_Match_Exception(Exception):
    """Error for data match."""

class Data_Load_Exception(Exception):
    """Error for data loading."""
