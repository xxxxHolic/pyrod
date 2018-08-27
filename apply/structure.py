# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 23:42:12 2018

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
from tool.properties import roughness
import tool.control as tc
import read.read_parameters as rp

####################### CONSTANT ############################

# constant 

####################### FUNCTIONS ###########################

def atomic_form_factor(ei_list,square_q,mode = 'factor'):
    
    """atomic form factors"""
    
    rp.check_ions(ei_list)
    
    if mode == 'factor':
    
        # atomic form factor
        af = rp.IONS_TABLE
        
        f = np.mat(np.zeros([np.size(square_q,0),len(ei_list)]))
        
        num = 0
        
        for ei in ei_list:
        
            f[:,num] = np.mat(af.at[ei,'a1']*np.exp(-af.at[ei,'b1']*square_q/4) +\
                              af.at[ei,'a2']*np.exp(-af.at[ei,'b2']*square_q/4) +\
                              af.at[ei,'a3']*np.exp(-af.at[ei,'b3']*square_q/4) +\
                              af.at[ei,'a4']*np.exp(-af.at[ei,'b4']*square_q/4) +\
                              af.at[ei,'c']).T
            num = num + 1
            
    elif mode == 'electron':
        
        # atomic form factor
        af = rp.IONS_TABLE
        
        f = np.mat(np.ones([np.size(square_q,0),len(ei_list)]))
        
        num = 0
        
        for ei in ei_list:
        
            f[:,num] = f[:,num]*(af.at[ei,'a1'] + af.at[ei,'a2'] + af.at[ei,'a3'] + af.at[ei,'a4'] + af.at[ei,'c'])
            
            num = num + 1
        
    return f

def debye_waller_factor(dw,q,u,v):
    
    """
    ______________________________________________________
    
       debye waller factors 
    """
    # debye waller factor obtained from parameters.xlsx
    # actually its thermal vibration distance of atoms
    # default setting is 0.1 A
#    dw = slab['dw']
    num = 0
    
    # debye waller factor matrix--sampling number x atoms number
    d = np.mat(np.zeros([len(q),np.size(dw,1)]))
    
    for ei in dw.columns:
        
        ds = np.exp(-np.square(2*np.pi)*(\
                    (dw.at['dw_x',ei]*u)**2+\
                    (dw.at['dw_y',ei]*v)**2+\
                    (dw.at['dw_z',ei]*q)**2)/2)
        
        d[:,num] = np.mat(ds).T
        num = num + 1
    
    return d

def substrate_ctr(var_list, iq, 
                  sub_ubr,
                  sub_index,
                  beta = 0, u = 0, v = 0, scale = 1e-4,
                  key_list = ['absorption',
                              'pos',
                              'roughness',
                              'dw',
                              'scale']):

    """
       calculate the ctr rods of substrate
    """
    
    square_q = iq**2
    # sn for sampling number
    sn = len(iq)
    
    atoms = var_list['substrate']['pos']
    ions = (np.array(var_list['substrate']['ions']).tolist())[0]
    
    p = atoms.as_matrix()
    # thermal vibration distance
    dw = var_list['substrate']['dw']
    
    # r-- crystal surface roughness
    r = np.mat(np.ones([sn,1]))
    
    "---------------------------------------------------------"
    
    # introduce the contral variable list
    # indicator to indicate the variable postion in ubr list
    
    for key in key_list:
        
        if key == 'absorption':
            # aor-- absorption indicator
            aor = tc.indicator_contral_vars(sub_ubr,sub_index,'substrate','absorption')
            beta = beta + (1-sub_ubr[aor[0]:aor[1]][0])
            
        elif key == 'roughness':
            # ror--roughness indicator
            ror = tc.indicator_contral_vars(sub_ubr,sub_index,'substrate','roughness')
            # beta_r- roughness beta distinguish with absorption beta
            beta_r = sub_ubr[ror[0]:ror[1]][0]-1
            cr = roughness(iq)
            r = np.mat(cr.robinson_roughness(beta_r)).T
            
        elif key == 'scale':
            # sor-- scale indicator
            sor = tc.indicator_contral_vars(sub_ubr,sub_index,'substrate','scale')
            scale = scale + sub_ubr[sor[0]:sor[1]][0]
            
        elif key == 'pos':
            # por-- pos indicator
            por = tc.indicator_contral_vars(sub_ubr,sub_index,'substrate','pos')
            # get pos ubr and reshape it into 3x5 matrix p.shape
            p_ubr = np.reshape(np.array(sub_ubr[por[0]:por[1]]),p.shape)
            p = np.multiply(p, p_ubr)
            
        elif key == 'dw':
            # dor-- debye waller factor indicator
            dor = tc.indicator_contral_vars(sub_ubr,sub_index,'substrate','dw')
            d_ubr = np.reshape(np.array(sub_ubr[dor[0]:dor[1]]),dw.shape)
            dw = np.multiply(dw, d_ubr)
            
    "---------------------------------------------------------"
    
    # a--atomic form factor
    a = atomic_form_factor(ions,square_q)
    
    # d--debye waller factor
    d = debye_waller_factor(dw,iq,u,v)
    
    # q--construct 3d qs space
    q = np.mat(np.zeros([sn,3]))
    
    q[:,0] = np.mat(np.ones([sn,1]))*u # qx
    q[:,1] = np.mat(np.ones([sn,1]))*v # qy
    q[:,2] = np.mat(iq).T
    
    # f--structure form factor
    f = np.multiply(np.multiply(a,d),
                    np.exp(1j*2*np.pi*(np.dot(q,p))))
    
    # s--shape function
    s = 1/(1 - np.exp(-1*1j*2*np.pi*np.dot(q,np.mat([0,0,1]).T))*\
                                                    np.exp(-beta) +\
                                                    scale)
    
    #sctr--structure form factor for crystal turncation rod
    sctr = np.multiply(f,s)
    
    return np.multiply(np.sum(sctr,1), r)

def film_ctr(varl,vart,subr,sindex, iq, u, v,key_list = ['lattice_abc','dw','pos']):
    
    "--------------------------------------------------------------------------"

    for key in key_list:
        
        # lattice abc modulation
        if key == 'lattice_abc':
            # lattice contral variable
            pos_ctable = []
            # layer loop
            for l in vart['slab_list']:
               # the indicator of the lattice constant contral variable for each layer
                pos_indicator = tc.indicator_contral_vars(subr,sindex,l,key)
                pos_contral = subr[pos_indicator[0]:pos_indicator[1]]
                pos_ctable.append(pos_contral)
                
        # debye_waller_factor modultion
        elif key == 'dw':
            # dw contral table
            dw_ctable = []
            # dw contral variable
            for l in vart['slab_list']:
                # indicator of debye waller factor for each layer
                dw_indicator = tc.indicator_contral_vars(subr,sindex,l,key)
                dw_contral = subr[dw_indicator[0]:dw_indicator[1]]
                dw_ctable.append(np.reshape(dw_contral,(3,varl[l]['atoms_num'])))
                
        # posz--modulate lattice c, modualte it indivual
            
    "--------------------------------------------------------------------------"

    m,n = vart['posz_table'].shape
    u,v = 0,0
    sn = len(iq)
    # Note! data type = complex!
    #  ComplexWarning: Casting complex values to real discards the imaginary part 
    slab_ctr = np.mat(np.zeros([sn,1],dtype = complex))
    
    # q--construct 3d qs space
    q = np.mat(np.zeros([sn,3]))
    
    q[:,0] = np.mat(np.ones([sn,1]))*u # qx
    q[:,1] = np.mat(np.ones([sn,1]))*v # qy
    q[:,2] = np.mat(iq).T
        
    square_q = np.square(iq)
    
    # loop for layer
    
    tot_layer = 0
    
    for layeri in range(m):
        
        # debye waller factor
        if 'dw' in key_list:
            dw_factor = np.multiply(vart['dw_table'][layeri],dw_ctable[layeri])
        else:
            dw_factor = vart['dw_table'][layeri]
        dw_layeri = debye_waller_factor(dw_factor,iq,u,v)
        
        # lattice constant modulation
        slabi = vart['slab_list'][layeri]
        
        if 'lattice_abc' in key_list:
            lattice_abc = pos_ctable[layeri]
        else:
            lattice_abc = [1,1,1]
            
        sub_c = varl['substrate']['lattice'].at['constant','c']

        # multivarible optimization of postion
        if 'pos' in key_list:        
            slab_c = varl[slabi]['lattice'].at['constant','c']*lattice_abc[2]*subr[layeri]
        else:
            slab_c = varl[slabi]['lattice'].at['constant','c']*lattice_abc[2]
            
        tot_layer = tot_layer + slab_c/sub_c
        
        # loop for atom
        for atomi in range(n):

            posi = [vart['posx_table'][layeri,atomi]*lattice_abc[0],
                    vart['posy_table'][layeri,atomi]*lattice_abc[1],
                    tot_layer + varl[slabi]['pos'].as_matrix()[2,atomi] - \
                    varl['substrate']['layers_num']*1]
            
            # the ocupany, contral vancany
            ocui = vart['ocu_table'][layeri,atomi]
            # ion name
            ion = vart['ion_table'][layeri][atomi]
            #
            a = ocui*atomic_form_factor([ion],square_q)
            d = dw_layeri[:,atomi]
            f = np.multiply(np.multiply(a,d),
                            np.exp(1j*2*np.pi*(np.dot(q,np.mat(posi).T))))
            
            slab_ctr = slab_ctr + f
            
    return slab_ctr

######################## CLASSS #############################

    