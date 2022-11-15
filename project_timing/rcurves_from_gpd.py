#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" 
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl

Generates samples of return curves for query return frequencies based on GPD fit and various other input parameters.

"""
import numpy as np
from GPD_Z_from_F import Z_from_F_mhhw, Z_from_F_loc,Z_from_F_Sweet

def rcurves_from_gpd_cov(location,scale,shape,cov,avg_exceed,num_mc_samps,sample_freqs,extrap_method,mhhw=None,mhhw_freq=None):
    """
    Generate return curves based on GPD fit and covariance matrix of scale & shape samples
    
    Input:
    location, scale, shape:     best-estimate GPD parameters
    cov:                        covariance matric scale & shape parameters
    avg_exceed:                 average annual exceedance rate of threshold (location parameter)
    num_mc_samps:               number of return curve samples to generate
    sample_freqs:               return frequencies to compute return height for
    extrap_method:              'None', 'Gumbel' or 'Sweet' -> how to deal with return frequencies higher than the average annual exceedance rate (i.e., return heights below the threshold)
        required for 'Gumbel':
            mhhw:               height of MHHW
            mhhwfreq:           frequency of MHHW (365/2)
            
    Output:
    rc_ce:                      central estimate return curve using central estimate GPD parameters
    rc_samlpes:                 samples using central estimate GPD parameters and their uncertainty
    
    *resulting return heights are given relative to the inputted location parameter
    """
    #generate scale and shape samples from covariance matrix assuming multivariate normal distributions
    gp_samples = np.random.multivariate_normal([shape,scale], cov,size=num_mc_samps)
    shape_samples = gp_samples[:,0]
    scale_samples = gp_samples[:,1]
    scale_samples[scale_samples<0] = 0 #no negative scales
    
    #expand these arrays according to length of input return frequencies for matrix operations
    scale_mat    = np.matlib.repmat(scale_samples,len(sample_freqs),1) 
    shape_mat    = np.matlib.repmat(shape_samples,len(sample_freqs),1)
    sample_freq_mat = np.transpose(np.matlib.repmat(sample_freqs,num_mc_samps,1)) #expand input frequencies according to number of samples

    #generate return curves, samples & central estimate
    if extrap_method == 'None':
        rc_samples = Z_from_F_loc(scale_mat,shape_mat,location,avg_exceed,sample_freq_mat)
        rc_ce = Z_from_F_loc(scale,shape,location,avg_exceed,sample_freqs)
    elif extrap_method == 'Gumbel':
        rc_samples = Z_from_F_mhhw(scale_mat,shape_mat,location,avg_exceed,sample_freq_mat,mhhw,mhhw_freq)
        rc_ce = Z_from_F_mhhw(scale,shape,location,avg_exceed,sample_freqs,mhhw,mhhw_freq)
    elif extrap_method == 'Sweet':
        rc_samples = Z_from_F_Sweet(scale_mat,shape_mat,location,avg_exceed,sample_freq_mat)
        rc_ce = Z_from_F_Sweet(scale,shape,location,avg_exceed,sample_freqs)
    return rc_ce,rc_samples

def rcurves_from_gpd_bootstrapped(location,scale,shape,avg_exceed,sample_freqs,extrap_method,scale_samples,shape_samples,mhhw=None,mhhw_freq=None):
    """
    Generate return curves based on GPD fit and uncertainty derived from bootstrapping (e.g., as applied in the automatic threshold selection of Solari et al.)
    
    Input:
    location, scale, shape:     best-estimate GPD parameters
    scale_samples:              samples of shape parameter obtained by applying GPD to bootstrapped GESLA observations
    shape_samples:              samples of shape parameter obtained by applying GPD to bootstrapped GESLA observations
    avg_exceed:                 average annual exceedance rate of threshold (location parameter)
    num_mc_samps:               number of return curve samples to generate
    sample_freqs:               return frequencies to compute return height for
    extrap_method:              'None', 'Gumbel' or 'Sweet' -> how to deal with return frequencies higher than the average annual exceedance rate (i.e., return heights below the threshold)
        required for 'Gumbel':
            mhhw:               height of MHHW
            mhhwfreq:           frequency of MHHW (365/2)
            
    Output:
    rc_ce:                      central estimate return curve using central estimate GPD parameters
    rc_samlpes:                 samples using central estimate GPD parameters and their uncertainty
    
    *resulting return heights are given relative to the inputted location parameter
    """
    scale_samples[scale_samples<0] = 0 #no negative scales
    
    #expand these arrays according to length of input return frequencies for matrix operations
    scale_mat    = np.matlib.repmat(scale_samples,len(sample_freqs),1) 
    shape_mat    = np.matlib.repmat(shape_samples,len(sample_freqs),1)
    sample_freq_mat = np.transpose(np.matlib.repmat(sample_freqs,len(scale_samples),1)) #expand input frequencies according to number of samples

    #generate return curves, samples & central estimate
    if extrap_method == 'None':
        rc_samples = Z_from_F_loc(scale_mat,shape_mat,location,avg_exceed,sample_freq_mat)
        rc_ce = Z_from_F_loc(scale,shape,location,avg_exceed,sample_freqs)
    elif extrap_method == 'Gumbel':
        rc_samples = Z_from_F_mhhw(scale_mat,shape_mat,location,avg_exceed,sample_freq_mat,mhhw,mhhw_freq)
        rc_ce = Z_from_F_mhhw(scale,shape,location,avg_exceed,sample_freqs,mhhw,mhhw_freq)
    elif extrap_method == 'Sweet':
        rc_samples = Z_from_F_Sweet(scale_mat,shape_mat,location,avg_exceed,sample_freq_mat)
        rc_ce = Z_from_F_Sweet(scale,shape,location,avg_exceed,sample_freqs)
    return rc_ce,rc_samples