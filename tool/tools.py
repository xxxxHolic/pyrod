# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 21:50:58 2018

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
from math import factorial

# data processing library
import numpy as np

# pyrod library

####################### CONSTANT ############################

# constant 

####################### FUNCTIONS ###########################

'.......................optimise.........................'

# f - fitting data
# y - experiment data
# mask - mask data

def R_square(f, y, mask):
    
    if not len(f) == len(y) == len(mask):
        raise Data_Match_Exception('Please input equal length')
    
    def nplist(data):
        
        # check and transform data
        try:
            
            # check np array
            if isinstance(data, np.ndarray):
                pass
            # check list
            elif isinstance(data, list):
                rl = np.array(data)
            # check np mat
            elif isinstance(data, np.matrix):
                rl = np.asarray(data).reshape(-1)
            # for other unpoackable datatype
            else:
                # init a list first
                l = []
                # unpack raw data with for
                for e in data:
                    l.append(e)
                # trans to np array
                rl = np.array(l)
                
        # unknown type
        except Data_Format_Exception:
            
            print('unknown data type')
            
        return rl

    # tranform to np array; apply mask 
    rf, ry = nplist(f)*nplist(mask), nplist(y)*nplist(mask)

    # calculate r square
    ss_tot = np.sum((ry - np.sum(ry)/len(ry))**2)
    ss_res = np.sum((ry - rf)**2)
    
    r2 = 1 - ss_res/ss_tot
    
    return r2


def opt_step_brute(func,x0_range,grid_size = 10,step = 2):
    
    """
       Brute method is much too slow and big.
       However, its usefull and simple. To improve it, we try to step it
       
       x0_range: range of variable, [x1-,x1+],[x2-,x2+]
       currently,only two axes are avaialble
    """
    # current step is 3
    step = 3
    
    # grid_size and step have to be integer
    try:
        grid_size = int(grid_size)
        step = int(step)
        
    except ValueError:
        raise ValueError("grid_size and step have to be of type int")
    
    # one dimensional step brute method
    if len(x0_range) == 1:
        
        # store func(grid_data) result
        grid_list0 = []
        x0 = np.linspace(x0_range[0][0],x0_range[0][1],grid_size)
        
        # func(grid_data)
        for px in range(grid_size):
            grid_list0.append(func(x0[px]))
        # store min in step1
        min_idx = np.argmin(grid_list0)
        
        # continue step2
        grid_list1 = []
        x1 = x0[min_idx]
        delta = (abs(x0_range[0][1] - x0_range[0][0]))/grid_size
        
        x2 = np.linspace(x1-delta,x1+delta,grid_size)
        for sx in range(grid_size):
            grid_list1.append(func(x2[sx]))
            
        min_step2 = x2[np.argmin(grid_list1)]
        
    elif len(x0_range) == 2:
    
        # step1: grid the x0_range
        min_step1 = []
        au = np.linspace(x0_range[0][0],x0_range[0][1],grid_size)
        av = np.linspace(x0_range[1][0],x0_range[1][1],grid_size)
        
        # find minimum in xu and xv grid
        def grid_min(xu,xv):
            
            x0_grid = np.meshgrid(xu, xv)
            
            #grid list
            grid_list = np.mat(np.zeros([grid_size**2,3]))
            idx = 0
            
            # pu-- for postion in u axes
            for pu in range(grid_size):
                # pv--for postion in v axes
                for pv in range(grid_size):
                    
                    grid_list[idx,0] = x0_grid[0][pu,pv]
                    grid_list[idx,1] = x0_grid[1][pu,pv]
                    grid_list[idx,2] = func([x0_grid[0][pu,pv],
                                             x0_grid[1][pu,pv]])
                    idx = idx + 1
            # find the minimum in step1
            min_idx = np.argmin(grid_list[:,2])
            
            return grid_list[min_idx,:]
        
        # append the firt minimum before rocking
        min_step1.append(grid_min(au,av))
        
        # start rocking, try to avoid local minmum
        bu = au - (au[1]-au[0])/2
        bv = av - (av[1]-av[0])/2
        
        min_step1.append(grid_min(bu,bv))
        
        # step 2
        # step 2 new x range
        u_min = np.min([min_step1[0][0,0],
                        min_step1[1][0,0]])
        u_max = np.max([min_step1[0][0,0],
                        min_step1[1][0,0]])
        deta_u = u_max - u_min
        v_min = np.min([min_step1[0][0,1],
                        min_step1[1][0,1]])
        v_max = np.max([min_step1[0][0,1],
                        min_step1[1][0,1]])
        deta_v = v_max - v_min
        # new u and v
        cu = np.linspace(u_min-deta_u, u_min+deta_u, grid_size)
        cv = np.linspace(v_min-deta_v, v_min+deta_v, grid_size)
        
        min_step2 = grid_min(cu,cv).tolist()
    
    return min_step2
    
            
'......................smooth.........................'

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    
    """    
       Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
       The Savitzky-Golay filter removes high frequency noise from data.
       It has the advantage of preserving the original shape and
       features of the signal better than other types of filtering
       approaches, such as moving averages techniques.

      ----------
      .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
         Data by Simplified Least Squares Procedures. Analytical
         Chemistry, 1964, 36 (8), pp 1627-1639.
      .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
         W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
         Cambridge University Press ISBN-13: 9780521880688
    """
    
    # integer value
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
     
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
        
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    
    # precompute coefficients
    
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    
    # pad the signal at the extremes with
    # values taken from the signal itself
    
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
     
    return np.convolve( m[::-1], y, mode='valid')

######################## CLASSS #############################

    