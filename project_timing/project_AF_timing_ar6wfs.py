#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" 
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl

Computes the projected timing of frequency amplifications of extreme sea levels
using the GPD distributions fitted to daily maxima derived from GESLA3 and 
projected SLR from IPCC AR6. FLOPROS estimates are used as a benchmark frequency.
"""
import numpy as np
import xarray as xr
from rcurves_from_gpd import rcurves_from_gpd_cov,rcurves_from_gpd_bootstrapped
from datetime import datetime
from get_FLOPROS_nearest_to_points import get_FLOPROS_nearest_to_points
import numpy.matlib
import os

#### SET PATHS AND CONFIGURE THE SCRIPT
gpd_file    = '/Volumes/Naamloos/PhD_Data/GESLA3/GPD_fits/gesla3_gpd_daily_max_potATS_Solari.nc' #file with GPD parameters/uncertainties
ar6_dir     = '/Volumes/Naamloos/PhD_Data/AR6_projections/ar6-GESLA3-samples-total/' #path to sea-level projections of IPCC AR6
output_dir  = '/Volumes/Naamloos/PhD_Data/GESLA3/timing_projections/' #output directory

extrap_below_gpd= 'Sweet'   # whether/how to derive return heights below the support of the GPD, one of ['Sweet','Gumbel','None']

f_ref           = 'FLOPROS' #locally varying 'FLOPROS', or single numerical value
yr_ref          = 2022      #starting year to compute timing
ssps            = ['ssp126','ssp245','ssp370','ssp585','ssp126lowconf','ssp585lowconf'] #SSPs to process

sample_freqs    = 10**np.linspace(-4,2,num=2001) #return frequencies used to generate return curves with
output_pcts     = np.arange(0,100,2.5)  #percentiles to evaluate results on
output_qnts     = output_pcts/100       #quantiles to evaluate results on

num_mc_samps = 10000 #number of Monte Carlo samples to draw from GPD parameters and sea-level projections (max. 20,000)
####

#### LOAD INPUT DATA, INITIALIZE ARRAYS
gpd_fits = xr.open_dataset(gpd_file) #load GPD fits to daily maxima at GESLA3 tide gauges

rc_dists = np.zeros((len(gpd_fits.site),len(output_qnts),len(sample_freqs)))                                            #return curve distribution quantiles
required_slr_dists = np.nan+np.zeros((len(gpd_fits.site),len(output_qnts),len(sample_freqs)))                           #required SLR distribution quantiles
timing_pboxes = np.nan+np.zeros((len(gpd_fits.site),len(ssps),len(output_qnts),len(sample_freqs)))                      #AF timing distribution quantiles
projected_slr_pboxes = np.nan+np.zeros((len(gpd_fits.site),len(ssps),len(output_qnts),len(np.arange(yr_ref,2151))))     #pbox quantiles of projected SLR

rc_ces = np.zeros((len(gpd_fits.site),len(sample_freqs)))                                                               #central-estimate return curves
required_slr_ces = np.nan+np.zeros((len(gpd_fits.site),len(sample_freqs)))                                              #central-estimate required SLR
timing_ces = np.nan+np.zeros((len(gpd_fits.site),len(ssps),len(sample_freqs)))                                          #central-estimate AF timing

if f_ref=='FLOPROS': #if f_ref=FLOPROS, find FLOPROS estimates for all GESLA3 sites
    gesla_flopros = 1/get_FLOPROS_nearest_to_points(gpd_fits.lon.values,gpd_fits.lat.values) #1/RP
####

#### COMPUTE RETURN CURVES, REQUIRED SLR and TIMING OF AFs
for st,site in enumerate(gpd_fits.site.values): #loop over GESLA3 sites
    print('Current site: '+str(st+1)+'/'+str(len(gpd_fits.site.values))+' '+site)
    
    iSite = np.where(gpd_fits.site.values==site)[0][0] #iSite=st if looping over all sites in gpd_fits
    
    #load GPD parameters (MHHW not yet incorporated at the moment)
    location = gpd_fits.sel(site=site).location.values 
    shape = gpd_fits.sel(site=site).shape.values
    scale = gpd_fits.sel(site=site).scale.values
    avg_exceed = gpd_fits.sel(site=site).avg_exceed.values
    
    
    #compute return curves (central estimates and samples)
    if 'cov' in gpd_fits.variables: #use covariance matrix of scale/shape to generate return curve samples
        cov = gpd_fits.sel(site=site).cov.values
        rc_ce,rc_samples = rcurves_from_gpd_cov(location,scale,shape,cov,avg_exceed,num_mc_samps,
                                             sample_freqs,extrap_below_gpd)
    else: #use parameter uncertainty derived from bootstrapping (as with Solari et al. (2017) method)
        rc_ce,rc_samples = rcurves_from_gpd_bootstrapped(location,scale,shape,avg_exceed,sample_freqs,
                                        extrap_below_gpd,gpd_fits.sel(site=site).scale_samples.values,
                                        gpd_fits.sel(site=site).shape_samples.values)
    
    #store distribution quantiles derived from samples, and central estimate
    rc_dists[iSite,...] = np.quantile(rc_samples,output_qnts,axis=1)
    rc_ces[iSite,...] = rc_ce
    

    #compute required SLR for AFs
    if f_ref == 'FLOPROS':
        if gesla_flopros[iSite] == np.nan : #if region nearest to gauge has no FLOPROS, results become NaN, and continue with next station
            print('Has no FLOPROS')
            continue
                
        i_fref = np.argmin(np.abs(sample_freqs-gesla_flopros[iSite])) #index of sample frequency nearest to f_ref
    else:
        i_fref = np.argmin(np.abs(sample_freqs-f_ref)) #index of sample frequency nearest to f_ref
    
    #compute required SLR by subtracting the return height of f_ref from the return curves and inverting them (see Fig. 2 of main manuscript)
    required_slr_samples = -1*(rc_samples - rc_samples[i_fref,:])
    required_slr_ce = -1*(rc_ce - rc_ce[i_fref])
    
    #store distribution quantiles derived from samples, and central estimate (quantiles of positive and negative frequency amplifications are treated separately)
    required_slr_dists[iSite,:,0:i_fref] = np.quantile(required_slr_samples[0:i_fref,:],np.flip(output_qnts),axis=1)
    required_slr_dists[iSite,:,i_fref::] = np.quantile(required_slr_samples[i_fref::,:],output_qnts,axis=1) 
    required_slr_ces[iSite,...] = required_slr_ce
    
    
    #compute timing of AFs
    for sp,ssp in enumerate(ssps): #loop over SSPs
        print('Current SSP: '+ssp)
        
        #load samples of SLR projections of each IPCC AR6 workflow 
        if 'lowconf' in ssp:
            wfs = ['1f','2f','3f','4'] #workflows to use for low-confidence projections under ssp126 and ssp585
        else:
            wfs = ['1f','2f'] #workflows to use for medium-confidence projections
   
        for wf in wfs: #for each workflow
            site_ar6_raw_samples=[] #open raw samples for current site
            site_ar6_raw_samples = xr.open_dataset(os.path.join(ar6_dir,'wf_'+wf+'/'+ssp.replace('lowconf','')+'/total-workflow.nc')).isel(locations=iSite)
            site_ar6_raw_samples.load()
            
            #to avoid loading in data for all sites at each iteration, store projection samples for this site in separate NetCDF if not already existing
            if not os.path.exists(os.path.join(ar6_dir,'wf_'+wf+'/'+ssp.replace('lowconf','')+'/per_tg')):
                os.mkdir(os.path.join(ar6_dir,'wf_'+wf+'/'+ssp.replace('lowconf','')+'/per_tg'))
                
            if not os.path.exists(os.path.join(ar6_dir,'wf_'+wf+'/'+ssp.replace('lowconf','')+'/per_tg/'+site+'.nc')):
                site_ar6_raw_samples.to_netcdf(os.path.join(ar6_dir,'wf_'+wf+'/'+ssp.replace('lowconf','')+'/per_tg/'+site+'.nc'))
              
        #load samples for site for all considered workflows
        paths = [os.path.join(ar6_dir,'wf_'+wf+'/'+ssp.replace('lowconf','')+'/per_tg/'+site+'.nc') for wf in wfs]
        site_ar6_samples = []
        site_ar6_samples = xr.open_mfdataset(paths,combine='nested',concat_dim='wf').assign_coords(wf=wfs)
        
        if np.isnan(site_ar6_samples.sea_level_change.isel(years=-1).values).all(): #if site has no projections
            print('No projections available.')
            continue
        
        #select random subset of projection samples of size num_mc_samps from the 20,000 available samples
        np.random.seed(0)
        rand_subset = np.random.randint(0,20000-1,size=num_mc_samps) 
        
        #interpolate decadal projections to annual means, reference time series to 2022 and convert projections from mm to m
        site_ar6_samples = site_ar6_samples.isel(samples=rand_subset).interpolate_na(dim='years',method='quadratic',fill_value='extrapolate').interp(years=np.arange(yr_ref,2151),method='quadratic')
        site_ar6_samples = (site_ar6_samples.sea_level_change - site_ar6_samples.sea_level_change.sel(years=yr_ref))/1000
        
        site_ar6_dists = site_ar6_samples.quantile(output_qnts,dim='samples') #compute distribution quantiles from samples
        
        #compute timing for all considered workflows
        timing_wfs = np.zeros((num_mc_samps,len(sample_freqs),len(wfs))) #initialize for current site
        for w,wf in enumerate(wfs): #loop over workflows
            site_ar6_samples_wf = site_ar6_samples.sel(wf=wf).values #load samples into memory
            
            #interpolate timing of projected SLR onto required SLR, looping over all samples (could possibly be vectorized?)
            timing_wfs[:,:,w] = np.array([np.interp(required_slr_samples[:,sample], site_ar6_samples_wf[sample,:], np.arange(yr_ref,2151)) for sample in np.arange(0,num_mc_samps)])

        #compute central-estimate timing by combining central-estimate required SLR and pbox median of projection samples (=mean of medians of considered workflows)
        timing_ce = np.interp(required_slr_ce,site_ar6_dists.sel(quantile=.5).mean(dim='wf'),np.arange(yr_ref,2151))
        timing_ce[timing_ce==2150] = np.nan #set 2150 to NaN (time series stop in 2150)
        
        #compute workflow timing distributions (computes quantiles by sorting, this is necessary because time series end in 2150)
        timing_wfs_sorted = np.sort(timing_wfs,axis=0) #sort & set 2150 to NaN
        timing_wfs_sorted[timing_wfs_sorted==2150] = np.nan
        
        #(!) this is not generic and not 100% accurate, but it works for num_mc_samps=1e4 & 2 decimal quantiles
        timing_dists_wfs = (timing_wfs_sorted[(output_qnts*num_mc_samps).astype(int),...] + timing_wfs_sorted[(output_qnts*num_mc_samps).astype(int)-1,...])/2

        #generate pbox of workflow distributions
        median_idx = np.flatnonzero(output_qnts==.5) #generate indices for which quantile are the median, >median, and <median
        above_idx = np.arange(median_idx + 1, len(output_qnts))
        below_idx = np.arange(median_idx)
        
        pbox_timing = np.full(timing_dists_wfs.shape[0:-1],np.nan) #initialize pbox array
        
    	# Use the maximum of workflows for quantiles >median
        pbox_timing[above_idx,:] = np.amax(timing_dists_wfs[above_idx,:,:], axis=-1) #NaN is latest (means >2150), so OK if amax evaluates to NaN
        
        # Use the minimum of workflows for quantiles <median
        temp = timing_dists_wfs 
        temp[np.isnan(temp)]=2150 # as np.amin will evaluate to NaN, temporarily set NaN to 2150
        
        pbox_timing[below_idx,:] = np.nanmin(temp[below_idx,:,:], axis=-1)
        
        # Use the mean of the medians of the workflows as the median for the pbox
        pbox_timing[median_idx,:] = np.mean(temp[median_idx,:,:], axis=-1)
        
        pbox_timing[pbox_timing==2150] = np.nan #set 2150 to NaN
        
        #compute and store the probability box of the projected SLR itself
        pbox_projected_slr = np.full(site_ar6_dists.isel(wf=0).shape,np.nan)         
        pbox_projected_slr[median_idx,:] = np.mean(site_ar6_dists[median_idx,:,:],axis=1).values 
        pbox_projected_slr[below_idx,:] = np.amin(site_ar6_dists[below_idx,:,:],axis=1).values
        pbox_projected_slr[above_idx,:] = np.amax(site_ar6_dists[above_idx,:,:],axis=1).values
        
        #aggregate site results in larger arrays
        timing_pboxes[iSite,sp,:,:] = pbox_timing
        projected_slr_pboxes[iSite,sp,:,:] = pbox_projected_slr
        timing_ces[iSite,sp,:] = timing_ce
        ####

       
#### PREPARE NETCDF OUTPUT
 
#store populated arrays in xarray dataset
ds = xr.Dataset(
    data_vars=dict(
        timing_pbox=(["site", "ssp","pct","freq"], timing_pboxes),
        timing_ce=(['site','ssp','freq'],timing_ces),
        projected_slr_pbox=(["site", "ssp","pct","year"], projected_slr_pboxes),
        required_slr_dist=(["site", "pct","freq"], required_slr_dists),
        required_slr_ce=(["site","freq"], required_slr_ces),
        rc_dist=(["site","pct","freq"], rc_dists),
        rc_ce=(["site","freq"], rc_ces),
        lon=(["site"], gpd_fits.lon.values),
        lat=(["site"], gpd_fits.lat.values),
    ),
    coords=dict(
        site=gpd_fits.site.values,
        ssp=ssps,  
        pct=output_pcts,
        freq=sample_freqs,
        year=np.arange(yr_ref,2151)
    )
)

del timing_pboxes
del required_slr_dists

#add FLOPROS standards if used
if f_ref=='FLOPROS':
    ds['flopros'] = 0*ds.lon+gesla_flopros
    
#add metadata
ds.attrs['esl_file'] = gpd_fits.esl_file
ds.attrs['f_ref'] = str(f_ref)
ds.attrs['yr_ref'] = str(yr_ref)
ds.attrs['num_mc_samps'] = num_mc_samps
ds.attrs['extrapolation_below_gpd'] = extrap_below_gpd
ds.attrs['history'] = 'Created '+str(datetime.now().day)+'_'+str(datetime.now().month)+'_'+str(datetime.now().year)
ds.timing_pbox.attrs['workflows'] = str(wfs)
ds.projected_slr_pbox.attrs['workflows'] = str(wfs)

ds.to_netcdf(path=output_dir+'esl_AF_timing_ar6_2150wfs_yrRef'+str(yr_ref)+'_fRef'+str(f_ref)+'_muATS_extrap'+extrap_below_gpd+'.nc',mode='w')
####