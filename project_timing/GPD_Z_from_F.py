#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" 
Created on Tue Nov 9th 15:58:34 2021

Functions to uses GPD parameters to estimate the return heights Z corresponding 
to the queried return frequencies F. 

@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""
import numpy as np
from scipy import interpolate

def Z_from_F_loc(scale,shape,loc,avg_exceed,sample_freqs): 
    """
    Obtain Z for F, with Z evaluating to NaN for F lower than the average 
    exceedance of the location parameter. Takes into account the lower and 
    upper supported values of the GPD distribution.
    
    Input: 
        scale:          scale parameter (scalar or vector)
        shape:          shape parameter (scalar or vector)
        loc:            location parameter
        avg_exceed:     exceedance frequeny of location parameter (POT threshold)
        sample_freqs:   sample frequencies to compute Z for
        
    Output:
        z:              heights computed for sample frequencies, relative to location parameter
    """
    try: #convert array input to float if using only 1 shape & scale
        scale=scale.item()
        shape=shape.item()
    except:
        pass
    
    z = 0*sample_freqs #initialize z
    
    if np.isscalar(shape): #if shape is scalar value
        if shape ==0: 
            z = -scale * np.log(sample_freqs/avg_exceed) + loc #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
        else:
            z = scale/shape * (np.power( (sample_freqs/avg_exceed), (shape/-1)) -1 ) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
            
    else: #if array input
        if np.any(shape)==0:
            z[shape==0] = -scale[shape==0] * np.log(sample_freqs[shape==0]/avg_exceed)
            z[shape!=0] = scale[shape!=0]/shape[shape!=0] * (np.power( (sample_freqs[shape!=0]/avg_exceed), (shape[shape!=0]/-1)) -1 )   
        else:
            z = scale/shape * (np.power( (sample_freqs/avg_exceed), (shape/-1)) -1 )
            
    #test lower bound of GPD (F=avg_exceed for z->0 with z=z0-loc (z0->loc))
    z = np.where(sample_freqs>avg_exceed,np.nan,z) #where F>avg_exceed, replace z by np.nan
    
    #test upper bound of GPD ((shape*z/scale)<-1)
    z = np.where((shape*z/scale)<-1,np.nan,z)
    return z
    
def Z_from_F_mhhw(scale,shape,loc,avg_exceed,sample_freqs,mhhw,mhhw_freq): 
    """
    Obtain Z for F, with Z for F lower than the average exceedance of the 
    location parameter evaluated using a Gumbel distribution between
    the location parameter and MHHW (with frequency of once per two days). 
    Takes into account the lower and upper supported values of the GPD distribution.
    
    Input: 
        scale:          scale parameter (scalar or vector)
        shape:          shape parameter (scalar or vector)
        loc:            location parameter
        avg_exceed:     exceedance frequeny of location parameter (POT threshold)
        sample_freqs:   sample frequencies to compute Z for
        mhhw:           z corresponding to MHHW
        mhhw_freq:      frequency corresponding to MHHW
    Output:
        z:              heights computed for sample frequencies, relative to location parameter 
                        (meaning z is negative for sample_freq>avg_exceed, i.e., in Gumbel distribution below location parameter)
    """
    try: #convert array input to float if using only 1 shape & scale
        scale=scale.item()
        shape=shape.item()
    except:
        pass
        
    z = 0*sample_freqs #initialize z
    
    if np.isscalar(shape): #if shape is scalar value
        if shape ==0: 
            z = -scale * np.log(sample_freqs/avg_exceed) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
        else:
            z = scale/shape * (np.power( (sample_freqs/avg_exceed), (shape/-1)) -1 ) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
            
    else: #if array input
        if np.any(shape)==0:
            z[shape==0] = -scale[shape==0] * np.log(sample_freqs[shape==0]/avg_exceed)
            z[shape!=0] = scale[shape!=0]/shape[shape!=0] * (np.power( (sample_freqs[shape!=0]/avg_exceed), (shape[shape!=0]/-1)) -1 )   
        else:
            z = scale/shape * (np.power( (sample_freqs/avg_exceed), (shape/-1)) -1 )

    #test upper bound of GPD ((shape*z/scale)<-1)
    z = np.where((shape*z/scale)<-1,np.nan,z)
    
    #below lower bound of GPD (F=avg_exceed for z->0 with z=z0-loc (z0->loc)), use Gumbel distribution
    z[sample_freqs>avg_exceed] = np.log(sample_freqs[sample_freqs>avg_exceed]/avg_exceed) * (mhhw-loc)/np.log(mhhw_freq/avg_exceed) #lower bound of gumbel is z=z0-loc=(mhhw-loc), i.e. (loc-mhhw) below loc if mhhw follows from long-term mean of 2-day maximum

    #below lower bound of Gumbel?
    z = np.where(sample_freqs>mhhw_freq,np.nan,z) #where F>mhhwFreq, replace z by np.nan

    return z


def Z_from_F_Sweet(scale,shape,loc,avg_exceed,sample_freqs): 
    """
    Obtain Z for F, with Z for F lower than the average exceedance of the 
    location parameter evaluated through extrapolation between .5/yr and the 
    average exceedance, up to 10/yr, following (Sweet et al., 2022).

    Input: 
        scale:          scale parameter (scalar or vector)
        shape:          shape parameter (scalar or vector)
        loc:            location parameter
        avg_exceed:     exceedance frequeny of location parameter (POT threshold)
        sample_freqs:   sample frequencies to compute Z for
       
    Output:
        z:              heights computed for sample frequencies, relative to location parameter 
                        (meaning z is negative for sample_freq>avg_exceed, i.e., in Gumbel distribution below location parameter)
    """
    
    try: #convert array input to float if using only 1 shape & scale
        scale=scale.item()
        shape=shape.item()
    except:
        pass
    
    z = 0*sample_freqs #initialize z
    
    if np.isscalar(shape): #if shape is scalar value
        if shape ==0: 
            z = -scale * np.log(sample_freqs/avg_exceed) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
        else:
            z = scale/shape * (np.power( (sample_freqs/avg_exceed), (shape/-1)) -1 ) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
            
    else: #if array input
        if np.any(shape)==0:
            z[shape==0] = -scale[shape==0] * np.log(sample_freqs[shape==0]/avg_exceed)
            z[shape!=0] = scale[shape!=0]/shape[shape!=0] * (np.power( (sample_freqs[shape!=0]/avg_exceed), (shape[shape!=0]/-1)) -1 )   
        else:
            z = scale/shape * (np.power( (sample_freqs/avg_exceed), (shape/-1)) -1 )

    #test upper bound of GPD ((shape*z/scale)<-1)
    z = np.where((shape*z/scale)<-1,np.nan,z)
    
    #below lower bound of GPD (F=avg_exceed for z->0 with z=z0-loc (z0->loc)), use extrapolation
    z = np.where(sample_freqs>avg_exceed,np.nan,z)
    
    z_0p5 = Z_from_F_loc(scale,shape,loc,avg_exceed,np.array([0.5])) #heights corresponding to 0.5/yr frequency
    z_0p2 = Z_from_F_loc(scale,shape,loc,avg_exceed,np.array([0.2])) #heights corresponding to 0.2/yr frequency
    
    
    iXtrp = ((sample_freqs>avg_exceed) & (sample_freqs<=10)) #indices of frequencies to fill with extrapolation
    
    if ((avg_exceed>0.5)&(avg_exceed<10)):
        if np.isscalar(shape): #if shape is scalar value
            z[iXtrp] = z_0p5+np.log(sample_freqs[iXtrp]/0.5) * (0-z_0p5)/np.log(avg_exceed/0.5) #log-linear similar to Gumbel
        else:
            z[iXtrp] = z_0p5[iXtrp]+np.log(sample_freqs[iXtrp]/0.5) * (0-z_0p5[iXtrp])/np.log(avg_exceed/0.5)
    
    elif ((avg_exceed>0.2)&(avg_exceed<10)):
        if np.isscalar(shape): #if shape is scalar value
            z[iXtrp] = z_0p2+np.log(sample_freqs[iXtrp]/0.2) * (0-z_0p2)/np.log(avg_exceed/0.2) #log-linear similar to Gumbel
        else:
            z[iXtrp] = z_0p2[iXtrp]+np.log(sample_freqs[iXtrp]/0.2) * (0-z_0p2[iXtrp])/np.log(avg_exceed/0.2)
    return z