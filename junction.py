# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 21:15:03 2018

@author: USER
"""

# Codes are free to use. Do whatever you want

from __future__ import absolute_import

"""junction of xrr and ctr data"""

####################### LIBRARY #############################

# exceptions library
from exceptions import (Data_Load_Exception)

# Python stdlib imports
import pickle
import os

# data processing library
import numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt

# pyrod library
from tool import control

####################### CONSTANT ############################

# data path 

CTR_PATH = os.path.abspath(os.path.dirname('ctr_optimised_result.pickle')) +\
            '/data/ctr_optimised_result.pickle'
            
XRR_PATH = os.path.abspath(os.path.dirname('xrr_optimised_result.pickle')) +\
            '/data/xrr_optimised_result.pickle'

FIT_PATH = os.path.abspath(os.path.dirname('fit.pickle')) +\
            '/data/fit.pickle'
            
# load data
            
ctr = {}
xrr = {}

try:
    
    # load ctr optimised result
    with open(CTR_PATH, 'rb') as c:
        unpickler = pickle.Unpickler(c)
        ctr = unpickler.load()
        
    # load xrr optimised result
    with open(XRR_PATH, 'rb') as x:
        unpickler = pickle.Unpickler(x)
        xrr = unpickler.load()
        
except Data_Load_Exception:
    print('Loading ctr or xrr optimised result error!')

# visulise data
PLINE = np.ones([len(ctr['shkl']), 3])

PLINE[:,0] = abs(ctr['substrate_ctr'] + ctr['slab_ctr'])
PLINE[:,1] = abs(xrr['rt']*xrr['inten'] + ctr['substrate_ctr'])
PLINE[:,2] = ctr['shkl']
        
####################### FUNCTIONS ###########################

# pick data point from a axes graph
def pick_points(line, number):
    
    # construct figure and axis handle
    fig, ax = plt.subplots()
    # add title
    ax.set_title('pick points', picker = True)
    # plot the , use 'o' to make it easyer to pick
#    l, = ax.plot(line, 'o','^','*', picker = number)
    l1,l2,l3, = ax.plot(np.log(line), 'o', picker = number)
    
    # x, y to store the xdata and ydata of picked points
    x = []
    y = []
    
    # the click function to operate event
    def click(event):
        
# The PickEvent which is passed to your callback is always fired with two attributes:
# mouseevent and artist
# Additionally, certain artists like Line2D and PatchCollection may attach additional 
# meta data like the indices into the data that meet the picker criteria 
# (e.g., all the points in the line that are within the specified epsilon tolerance)
        
        tline = event.artist
        ind = event.ind
        
        # pick xdata and ydata
        xdata = np.take(tline.get_xdata(), ind)[0]
        ydata = np.take(tline.get_ydata(), ind)[0]
        
        # store xdata and ydata
        x.append(xdata)
        y.append(ydata)
        
        # print the picked points
        print('x = ' + str(xdata))
        print('y = ' + str(ydata))
    
    # trigger the events
    fig.canvas.mpl_connect('pick_event', click)
   
    # loop to pause the figure, and pick the points
    while True:
        # before enough points are picked, pause
        if len(x) < number:
            # pause the figure
            plt.pause(0.5)
        # enough points are picked, break from the loop
        else:
            break
    
    plt.close(fig)
    
    return x,y

######################## CLASSS #############################

class connection(object):
    
    def __init__(self, bond):
        
        self.bond = bond
        self.data = np.zeros(len(ctr['q']))
        self.r_square = 0
    
    # show all the optimised data -- xrr, ctr and raw data
    def show(self):
        
        fig, ax = plt.subplots()
        
        line1, = ax.plot(ctr['q'], np.log(PLINE[:,0]), 'g')
        line2, = ax.plot(ctr['q'], np.log(PLINE[:,1]), 'y')
        line3, = ax.plot(ctr['q'], np.log(PLINE[:,2]), 'r')
    
        plt.legend((line1, line2, line3), ('ctr','xrr','raw data'))
    
    # connect xrr and ctr data near first bragg peak
    def connect(self, b = 0.98):
        
        self.bond = b
        
        con_data = np.zeros(len(ctr['q']))
        location = np.argmin(abs(ctr['q'] - self.bond))
        
        con_data[0: location] = PLINE[0: location, 1]
        con_data[location:, ] = PLINE[location:,   0]
        
        con_f = si.interp1d(ctr['q'], con_data)
        self.data = con_f(ctr['q'])
        
        plt.scatter(ctr['q'], np.log(PLINE[:,2]), s = 5, c = 'r')
        plt.plot(ctr['q'], np.log(PLINE[:,2]), c = 'r')
        plt.plot(ctr['q'], np.log(self.data))
    
    # calculate the r_square between raw data and fit y
    # limit: the location of start raw data
    def error(self, limit):
        
        mask = control.bragg_mask(ctr['q'], 6, 1, 'yin', limit)
        number = {i:mask.tolist().count(i) for i in mask}[1]
        
        yhat = np.log(abs(self.data))*mask
        ydat = np.log(abs(ctr['shkl']))*mask
        ybar = np.sum(ydat)/number
        
        ssreg = np.sum((yhat - ybar)**2)
        sstot = np.sum((ydat - ybar)**2)
        
        if ssreg >= sstot:
            r2 = sstot/ssreg
        else:
            r2 = ssreg/sstot
        
        self.r_square = r2
        
        return r2
        
    # save the final fitting data
    def save(self):
        
        fit = {'q': ctr['q'],
               'fit_y': self.data,
               'shkl': ctr['shkl'],
               'r_square': self.r_square}
        
        f = open(FIT_PATH, 'wb')
        pickle.dump(fit, f)
        f.close()