# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 00:07:08 2018

@author: USER
"""

# Codes are free to use. Do whatever you want

from __future__ import absolute_import

"""CTR data fitting"""

####################### LIBRARY #############################

# exceptions library

# Python stdlib imports
import pickle
import os

# data processing library
import numpy as np
import pandas as pd
import scipy.optimize as so
import matplotlib.pyplot as plt
import openpyxl as ox

# pyrod library

from read.read_parameters import initialization_parameters
from read.read_raw_data import initialization_rhkl

import tool.control as tc
import tool.tools as tt

from tool.properties import (roughness, vacancy, relax)
import structure

####################### CONSTANT ############################

# constant 
data_path = os.path.abspath(os.path.dirname('ctr_optimised_result.pickle')) +\
            '/data/ctr_optimised_result.pickle'
####################### FUNCTIONS ###########################

######################## CLASSS #############################

class var_refine(object):
    
# class for modulate the parameters maunally or opimisly

    def __init__(self, parameters, experiment_data):
        
        # init var_list and table
        self.parameters = parameters
        self.para_path = os.path.abspath(os.path.dirname(parameters)) + parameters
        self.experiment_data = experiment_data
        self.data_path = os.path.abspath(os.path.dirname(experiment_data)) + experiment_data
        
        init_para = initialization_parameters(self.parameters)
        init_para._check_parameters()
        init_para._check_sheetname()
        
        init_data = initialization_rhkl(self.experiment_data)
        init_data._check_filename()
        
        
        # real time read xlsx, real time refresh
        self.var_list = init_para._var_list()
        self.var_table = init_para._var_table()
        # import factors
        self.properties = pd.read_excel(self.para_path,sheet_name = 'factors',index_col = 0)
            
        # init experiment data iq and intensity
        self.iq,self.i = init_data._read_rhkl(experiment_data)
       
        # whick rod
        self.h = self.properties.at['h','data']
        self.k = self.properties.at['k','data']
        
        # inti parameters
        self.absorption = self.properties.at['absorption','data']
        self.roughness = self.properties.at['roughness','data']
        self.scale = self.properties.at['scale','data']
        self.intensity = self.properties.at['intensity','data']
        # interface roughness
        self.iroughness = self.properties.at['interface roughness','data']
        
        # contral parameters
        self.ubr_index = tc.contral_vars(self.var_list)
        
    def track_difference(self):
        
        var_list_rei = initialization_parameters._var_list(self.parameters)
        
        track_flag = 0
        
        try:
            
            for slabi in var_list_rei:
                for keyi in var_list_rei[slabi]:
                    
                    # compare every key value in new excel and old excel
                    if keyi in ['absorption','atoms_num','layers_num','roughness']:
                        # refreshed data and old data
                        value_rei = var_list_rei[slabi][keyi]
                        value_old = self.var_list[slabi][keyi]
                        
                        # if same,continue. if change,flag = 1
                        if np.equal(value_rei,value_old):
                            continue
                        else:
                            track_flag = 1
                    # these parameters are matrix
                    elif keyi in ['dis','dw','ions','lattice','nel','ocu','pos']:
                        # refreshed matrix and old matrix
                        if False in np.equal(var_list_rei[slabi][keyi],
                                             self.var_list[slabi][keyi]).as_matrix():
                            track_flag = 1
                            break
                        else:
                            continue

                    # beak from the inner loop
                    if track_flag:
                        break
                # break from the outer loop
                if track_flag:
                    break
                
        except:
            track_flag = 1
            
        properties = pd.read_excel(self.para_path,sheet_name = 'factors',index_col = 0)
        
        if False in np.equal(properties,self.properties).as_matrix():
            track_flag = 1
            
        # refresh
        if track_flag == 1:
            
            init_para = initialization_parameters(self.parameters)
            init_data = initialization_rhkl(self.experiment_data)
            
            self.var_list = init_para._var_list()
            self.var_table = init_para._var_table()
            # import factors
            self.properties = pd.read_excel(self.para_path,sheet_name = 'factors',index_col = 0)
                
            # init experiment data iq and intensity
            self.iq,self.i = init_data._read_rhkl(self.experiment_data)
           
            # whick rod
            self.h = self.properties.at['h','data']
            self.k = self.properties.at['k','data']
            
            # inti parameters
            self.absorption = self.properties.at['absorption','data']
            self.roughness = self.properties.at['roughness','data']
            self.scale = self.properties.at['scale','data']
            self.intensity = self.properties.at['intensity','data']
            # interface roughness
            self.iroughness = self.properties.at['interface roughness','data']
            # contral parameters
            self.ubr_index = tc.contral_vars(self.var_list)
            
        return self,track_flag
    
    def re(self):
        
        init_para = initialization_parameters(self.parameters)
        init_data = initialization_rhkl(self.experiment_data)
        
        self.var_list = init_para._var_list()
        self.var_table = init_para._var_table()
        
        init_para = initialization_parameters(self.parameters)
        init_data = initialization_rhkl(self.experiment_data)
        
        # import factors
        self.properties = pd.read_excel(self.para_path,sheet_name = 'factors',index_col = 0)
            
        # init experiment data iq and intensity
        self.iq,self.i = init_data._read_rhkl(self.experiment_data)
       
        # whick rod
        self.h = self.properties.at['h','data']
        self.k = self.properties.at['k','data']
        
        # inti parameters
        self.absorption = self.properties.at['absorption','data']
        self.roughness = self.properties.at['roughness','data']
        self.scale = self.properties.at['scale','data']
        self.intensity = self.properties.at['intensity','data']
        # interface roughness
        self.iroughness = self.properties.at['interface roughness','data']
            
        
    def disp(self):
        
        # list all the self parameters
        
        print("--------------------------------------------")
        print("parameters: %s\n" % self.parameters)
        print("experiment_data: %s\n" % self.experiment_data)
        print("h k qz: %s %s %s\n" % (int(self.h),
                                      int(self.k),
                                      [np.min(self.iq),np.max(self.iq)]))
        print("intensity: %s\n" % self.intensity)
        print("absorption: %s\n" % self.absorption)
        print("roughness: %s\n" % self.roughness)
        print("scale: %s\n" % self.scale)
        print("var_list: %s\n" %self.var_list.keys())
        print('var_table: %s' %self.var_table.keys())
        print("roughness interface: %s" %self.iroughness)
        print('--------------------------------------------')
        
    def sub_ctr(self):
        
        # the used keys
        key_list = ['absorption', 'roughness', 'scale']
        # input substrate ubr and substrate index
        subr,sindex = tc.initialize_contral_vars(self.ubr_index,
                                                 ['substrate'],
                                                 key_list)
        
        substrate_ctr = structure.substrate_ctr(self.var_list,
                                                self.iq,
                                                subr,
                                                sindex,
                                                self.absorption,
                                                self.h,
                                                self.h,
                                                self.scale,
                                                key_list) 
        
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        
        return np.multiply(self.intensity*substrate_ctr, np.mat(r).T)
    
    def slab_ctr(self):
        
        # used keys
        key_list = ['lattice_abc','dw']
        # input slab ubr and slab index
        slab_ubr,slab_index = tc.initialize_contral_vars(self.ubr_index,
                                                         np.unique(self.var_table['slab_list']).tolist(),
                                                         key_list)
        fctr = structure.film_ctr(self.var_list,
                                  self.var_table,
                                  slab_ubr,
                                  slab_index,
                                  self.iq,
                                  self.h,
                                  self.k,
                                  key_list)
        
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        # interface roughness
        ir = cr.interface_roughness(self.iroughness)
        
        return np.multiply(np.multiply(self.intensity*fctr, np.mat(r).T),np.mat(ir).T)
    
    def lattice_refine(self, p_mask = 13, n_mask = 3, weight = 1):
        
        # Prepare the data mask
        # yin mask to mask the bragg peaks 
        bn_mask = np.mat(tc.bragg_mask(self.iq,
                                       n_mask,1,
                                       mode = 'yin'))
        # yang mask to refine the signal around bragg peaks only
        bp_mask = np.mat(tc.bragg_mask(self.iq,
                                       p_mask,1,
                                       mode = 'yang'))
        bgm  = np.asarray(np.multiply(bn_mask, bp_mask).T).reshape(-1)*weight
        
        # calculate ss tot.For only seveal signal points are calculate
        # mean value is np.sum(abs(i)*bgm)/np.sum(bgm)
        ss_tot = np.sum((abs(self.i)*bgm - \
                         np.sum(abs(self.i)*bgm)/np.sum(bgm))**2)
        
        # initialize the variable about the lattice
        slab_list = np.unique(self.var_table['slab_list']).tolist()        
        subr,sindex = tc.initialize_contral_vars(self.ubr_index,
                                                 slab_list,
                                                 ['lattice_abc'])
        # substrate ctr
        # the used keys
        key_list = ['absorption', 'roughness', 'scale']
        # input substrate ubr and substrate index
        sub_ubr,sub_index = tc.initialize_contral_vars(self.ubr_index,
                                                       ['substrate'],
                                                       key_list)
        
        substrate_ctr = structure.substrate_ctr(self.var_list,
                                                self.iq,
                                                sub_ubr,
                                                sub_index,
                                                self.absorption,
                                                self.h,
                                                self.h,
                                                self.scale,
                                                key_list) 
    
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        
        def c_refine(c_var):
            
            # input lattice c variable into subr
            if len(slab_list) == 1:
                subr[2] = c_var
            elif len(slab_list) >= 2:
                for slabi in range(len(slab_list)):
                    subr[slabi*3-1] = c_var[slabi]
                
            fctr = structure.film_ctr(self.var_list,
                                      self.var_table,
                                      subr,
                                      sindex,
                                      self.iq,
                                      self.h,
                                      self.k,
                                      ['lattice_abc'])
            
            ss = np.multiply(np.mat(fctr) + substrate_ctr,np.mat(r).T)*self.intensity
            # calculate r square
            sa = np.asarray(ss).reshape(-1)
            ss_res = np.sum((abs(self.i)*bgm - abs(sa)*bgm)**2)
            varience = ss_res/ss_tot
            
            print(int(varience*1e4)/1e4)
            
            return varience
        
        # select properity optimized method for lattice_c optimize
            
        # if there are only one or two variables, step_brute method is a direct method
        if len(slab_list) == 1:
            print('lattice c optimising....\nOPT_STEP_BRUTE method is used')
            rec = tt.opt_step_brute(c_refine,[[0.8,1.2]],grid_size = 20)
        elif len(slab_list) == 2:
            print('lattice c optimising....\nOPT_STEP_BRUTE method is used')
            rec = tt.opt_step_brute(c_refine,[[0.8,1.2],[0.8,1.2]],grid_size = 10)
            
        # for larger variable number, Nelder-Mead method is more comparable
        elif len(slab_list) >= 3:
            print('lattice c optimising....\nNELDER_MEAD method is used')
            c0 = np.ones(len(slab_list)).tolist()
            resc = so.minimize(c_refine,c0,method = 'Nelder-Mead',tol=10)
            rec = resc.x
        
        # update subr
        if len(slab_list) == 1:
            subr[2] = rec
        elif len(slab_list) >= 2:
            for slabi in range(len(slab_list)):
                subr[slabi*3-1] = rec[0][slabi]

        # update ubr_index
        self.ubr_index = tc.refresh_index(self.ubr_index,subr,sindex)
                
        # plot the refined result
        
        fctr = structure.film_ctr(self.var_list,
                                  self.var_table,
                                  subr,
                                  sindex, 
                                  self.iq, 
                                  self.h, 
                                  self.k,
                                  ['lattice_abc'])
        
        ss = np.multiply(np.mat(fctr) + substrate_ctr,np.mat(r).T)
        plt.plot(self.iq,np.log(abs(ss*self.intensity)))
        plt.plot(self.iq,np.log(np.sqrt(self.i)))
    
    def dw_refine(self,
                  p_mask = 13,
                  n_mask = 3,
                  weight = 1,
                  key_list = ['dw']):
        
        # Prepare the data mask
        # yin mask to mask the bragg peaks 
        bn_mask = np.mat(tc.bragg_mask(self.iq,
                                       n_mask,1,
                                       mode = 'yin'))
        # yang mask to refine the signal around bragg peaks only
        bp_mask = np.mat(tc.bragg_mask(self.iq,
                                       p_mask,1,
                                       mode = 'yang'))
        bgm  = np.asarray(np.multiply(bn_mask, bp_mask).T).reshape(-1)*weight
        
        # calculate ss tot.For only seveal signal points are calculate
        # mean value is np.sum(abs(i)*bgm)/np.sum(bgm)
        ss_tot = np.sum((abs(self.i)*bgm - 
                         np.sum(abs(self.i)*bgm)/np.sum(bgm))**2)
        
        # initialize the variable about the lattice
        slab_list = np.unique(self.var_table['slab_list']).tolist()        
        ubr_index = tc.contral_vars(self.var_list)
        subr,sindex = tc.initialize_contral_vars(ubr_index,
                                                 slab_list,
                                                 ['dw'])
        
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        
        # substrate ctr
        # the used keys
        key_list = ['absorption', 'roughness', 'scale']
        # input substrate ubr and substrate index
        sub_ubr,sub_index = tc.initialize_contral_vars(self.ubr_index,
                                                       ['substrate'],
                                                       key_list)
        
        substrate_ctr = structure.substrate_ctr(self.var_list,
                                                self.iq,
                                                sub_ubr,
                                                sub_index,
                                                self.absorption,
                                                self.h,
                                                self.h,
                                                self.scale,
                                                key_list) 
    
        def d_refine(d_var):
            
            fctr = structure.film_ctr(self.var_list,
                                      self.var_table,
                                      d_var,
                                      sindex,
                                      self.iq,
                                      self.h,
                                      self.k,
                                      ['dw'])
            
            ss = np.multiply(np.mat(fctr) + 
                             substrate_ctr,np.mat(r).T)*self.intensity
            # calculate r square
            sa = np.asarray(ss).reshape(-1)
            ss_res = np.sum((abs(self.i)*bgm - abs(sa)*bgm)**2)
            varience = ss_res/ss_tot
            
            print(int(varience*1e4)/1e4)
            
            return varience
        
        d0 = subr
        re = so.minimize(d_refine,d0,method = 'Nelder-Mead',tol=1)
        
        self.ubr_index = tc.refresh_index(self.ubr_index,re.x,sindex)
        
        # plot the refined result
        
        fctr = structure.film_ctr(self.var_list,
                                  self.var_table,
                                  re.x,
                                  sindex,
                                  self.iq,
                                  self.h,
                                  self.k,
                                  ['dw'])
                
        ss = np.multiply(np.mat(fctr) + substrate_ctr,np.mat(r).T)
        plt.plot(self.iq,np.log(abs(ss*self.intensity)))
        plt.plot(self.iq,np.log(np.sqrt(self.i)))
        
    
    # return the updated ubr_index, some parameters are list in ubr_index
    def update_var(self):
        
        # pl--parameters location in excel file
        pl = {'absorption':'B3',
              'layers_num':'B4',
              'lattice':'B6:G6',
              'lattice_abc':'B6:D6',
              'ocu':'B10:F10',
              'dw':'B12:F14',
              'dis':'B16:F18',
              'pos':'B22:F24'}
        
        #p x--parameters xlsx workbook
        px = ox.load_workbook(self.para_path)
        
        # check the slab in ubr_index
        for slabi in self.ubr_index:
            
            # intensity is not a slab should be picked out
            if slabi == 'intensity':
                continue
            
            # start update the slab data
            else:
                
                # parameters loop
                for keyi in self.ubr_index[slabi]:
                    
                    # only several parameters are updated from ubr_index
                    # scale roughness vacancy veta are not
                    if keyi in ['scale','roughness','vacancy','beta']:
                        pass
                    
                    else:
                        # key is the parameters from ubr_index
                        key = np.mat(self.ubr_index[slabi][keyi])
                        # psheet is the slabi sheet
                        psheet = px[slabi]
                        # cell is the cell of parameters in slabi sheet
                        cell = psheet[pl[keyi]]
                        
                        # parameters absorptio is a singe value
                        if keyi == 'absorption':
                            # key_var is the value of parameters from var_list
                            key_var = self.var_list[slabi][keyi]
                            # key is the contral variable from ubr_index
                            # key_var is the parameter value from excel file
                            cell.value = key[0,0]*key_var
                        
                        elif keyi == 'lattice_abc':
                            
                            
#                            print(key)
                            
                            # only apply lattice abc
                            key_var = self.var_list[slabi]['lattice']
                            
#                            print(key_var)
                            
                            ci = 0
                            # loop cell and update value
                            # in openpyxl, cell can only be valued in this way
                            for row in cell:
                                for element in row:
                                    
                                    if ci <= 2:
                                        element.value = key[0, ci]*key_var.as_matrix()[0,ci]
                                        
                                        print(element.value)
                                        ci += 1
                                    else:
                                        pass
                                    
                        else:
                            # for matrix parameters
                            key_var = self.var_list[slabi][keyi]
                            
                            # the row index
                            r = 0
                            # the column index
                            c = 0
                            
                            for row in cell:
                                for element in row:
                                    # apply the contral variable to origin data
                                    element.value = key[r, c]*key_var.as_matrix()[r, c]
                                    
                                    # re init the column index
                                    if c == key.shape[1]-1:
                                        c = c - key.shape[1]+1
                                    else:
                                        # loop column
                                        c += 1
                                # loop row
                                r += 1
        try:
            px.save(self.para_path)
        except:
            print('Error!!Please close the excel file %s and try again'%self.parameters)
            
    def pos_modulate(self, bc, sc, cl):
        
        subr0 = self.var_table['posz_list']
        slab_index = 'posz'
        
        subr0 = relax.strain_relax(subr0, cl, bc, sc)
        key_list = ['pos']
        
        print(subr0)
        
        fctr = structure.film_ctr(self.var_list,
                                  self.var_table,
                                  subr0,
                                  slab_index,
                                  self.iq,
                                  self.h,
                                  self.k,
                                  key_list)
        
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        # interface roughness
        ir = cr.interface_roughness(self.iroughness)
        
        return np.multiply(np.multiply(self.intensity*fctr, np.mat(r).T),np.mat(ir).T)
    
    def surface_modulate(self,thrface = 0.9,secface = 0.9,surface = 1.1):
        
        subr0 = self.var_table['posz_list']
        slab_index = 'posz'
        
        subr0[-3] = thrface
        subr0[-1] = surface
        subr0[-2] = secface
        key_list = ['pos']
        
        self.var_table['posz_list'][ 0] = thrface
        self.var_table['posz_list'][-1] = surface
        self.var_table['posz_list'][ 1] = secface
        
#        print(subr0)
        
        fctr = structure.film_ctr(self.var_list,
                                  self.var_table,
                                  subr0,
                                  slab_index,
                                  self.iq,
                                  self.h,
                                  self.k,
                                  key_list)
        
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        # interface roughness
        ir = cr.interface_roughness(self.iroughness)
        
        return np.multiply(np.multiply(self.intensity*fctr, np.mat(r).T),np.mat(ir).T)
    
    def pos_refine(self):
        
        # the used keys
        sub_key_list = ['absorption', 'roughness', 'scale']
        # input substrate ubr and substrate index
        sub_ubr,sub_index = tc.initialize_contral_vars(self.ubr_index,
                                                       ['substrate'],
                                                       sub_key_list)
        
        substrate_ctr = structure.substrate_ctr(self.var_list,
                                                self.iq,
                                                sub_ubr,
                                                sub_index,
                                                self.absorption,
                                                self.h,
                                                self.h,
                                                self.scale,
                                                sub_key_list) 
        
        # roughness
        cr = roughness(self.iq)
        r = cr.robinson_roughness(self.roughness)
        
        s = np.multiply(self.intensity*substrate_ctr, np.mat(r).T)
        
        subr0 = self.var_table['posz_list']
        slab_index = 'posz'
        key_list = ['pos']
        
        mask = tc.bragg_mask(self.iq,3,1,mode = 'yin',limit = 24)
        
#        fig, ax = plt.subplots()
        
        def varience(subr):
            
            fctr = structure.film_ctr(self.var_list,
                                      self.var_table,
                                      subr,
                                      slab_index,
                                      self.iq,
                                      self.h,
                                      self.k,
                                      key_list)
            
            # roughness
            cr = roughness(self.iq)
            # interface roughness
            ir = cr.interface_roughness(self.iroughness)
            
            f = np.multiply(np.multiply(self.intensity*fctr, np.mat(r).T),np.mat(ir).T)
            
            a = s + f
            
            e = np.sqrt(self.i)
            
#            plt.cla()
#            plt.plot(np.log10(e))
#            plt.plot(np.log10(abs(a)))
            
            v = np.sum(abs(mask*(np.log10(abs(a)+1e-6) - np.log10(e+1e-6))))
            
            print(v)
            
            return v
        
        re = so.minimize(varience,subr0,
                         method = 'Nelder-Mead',
                         options = {'maxiter':50})
        subr1 = re.x
        
        fctr1 = structure.film_ctr(self.var_list,
                                   self.var_table,
                                   subr1,
                                   slab_index,
                                   self.iq,
                                   self.h,
                                   self.k,
                                   key_list)
        
        # roughness
        cr = roughness(self.iq)
        # interface roughness
        ir = cr.interface_roughness(self.iroughness)
        
        return np.multiply(np.multiply(self.intensity*fctr1, np.mat(r).T),np.mat(ir).T)
    
    # return the initlatised parameters
    def return_self(self):
        return self

class vr(var_refine):
    
    # plot experiment data
    def p_rhkl(self):
        
        plt.plot(self.iq, np.log(abs(self.i) + 1e-5))
        
    # plot sqrt experiment data
    def p_shkl(self):
        
        plt.plot(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)),'o')

    # plot substrate data
    def p_sctr(self):

        plt.plot(self.iq, np.log(abs(var_refine.sub_ctr(self))))

    # plot slab data
    def p_fctr(self):

        plt.plot(self.iq, np.log(abs(var_refine.slab_ctr(self))))

    # plot sumed data -- substrate data + slab data
    def p_actr(self):

        plt.plot(self.iq, np.log(abs(var_refine.sub_ctr(self) + var_refine.slab_ctr(self))))        
    
    # return the complex data of substrate, slab. The sqrt data of experiment data
    
    # compare the fitting data and experiment data
    def p_c(self):
        
        plt.plot(self.iq, np.log(abs(var_refine.sub_ctr(self) + var_refine.slab_ctr(self))))
        plt.scatter(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)),s = 5, c = 'r')
        plt.plot(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)), c = 'r')
    
    # plot the pos_modulate result
    def p_p(self, bc, sc, cl):
        
        plt.plot(self.iq, np.log(abs(var_refine.sub_ctr(self) + var_refine.pos_modulate(self, bc, sc, cl))))
        plt.scatter(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)),s = 5, c = 'r') 
        
    # plot the surface modulate result
    def p_s(self,thrface,secface,surface):
        
        plt.plot(self.iq, np.log(abs(var_refine.sub_ctr(self) + var_refine.surface_modulate(self, thrface, secface, surface))))
        plt.scatter(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)),s = 5, c = 'r')
        plt.plot(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)), c = 'r')
        
    # plot the pos_refine result
    def p_r(self):
        
        plt.plot(self.iq, np.log(abs(var_refine.sub_ctr(self) + var_refine.pos_refine(self))))
        plt.scatter(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)),s = 5, c = 'r')
        plt.plot(self.iq, np.log(np.sqrt(abs(self.i) + 1e-5)), c = 'r')
    
    # return the fitting and sqrt experiment data
    def rdata(self):
        
        return {'iq':np.asarray(self.iq).reshape(-1),
                'shkl':np.asarray(np.sqrt(self.i) + 1e-5).reshape(-1),
                'substrate_ctr':np.asarray(var_refine.sub_ctr(self)).reshape(-1),
                'slab_ctr':np.asarray(var_refine.slab_ctr(self)).reshape(-1)}
    
    # save modulated data to folder /data
    def save(self):
        
        ctr_optimised_result = {'q':np.asarray(self.iq).reshape(-1),
                                'shkl':np.asarray(np.sqrt(self.i) + 1e-5).reshape(-1),
                                'substrate_ctr':np.asarray(var_refine.sub_ctr(self)).reshape(-1),
                                'slab_ctr':np.asarray(var_refine.slab_ctr(self)).reshape(-1)}
        
        f = open(data_path, 'wb')
        pickle.dump(ctr_optimised_result, f)
        f.close()
        