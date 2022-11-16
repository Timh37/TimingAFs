#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 14:47:27 2022

Takes the output of VMS's Solari et al. (2017) analysis, stored as individual
netcdf files per site, and stores them in a single overarching netcdf.

@author: thermans
"""

import numpy as np
import xarray as xr
import os
from gesla import GeslaDataset

#open gesla data for metadata at each site
g3dir = '/Volumes/Naamloos/PhD_Data/GESLA3/GESLA3.0_ALL/' #gesla data
g3files = os.listdir(g3dir) #find files in GESLA directory
g3 = GeslaDataset(meta_file='/Volumes/Naamloos/PhD_Data/GESLA3/GESLA3_ALL.csv',data_path=g3dir)

ats_output_dir = '/Volumes/Naamloos/PhD_Data/GESLA3/Solari_GESLA3_OUTPUT_30_04_22/'
ats_files = [i for i in os.listdir(ats_output_dir) if i.startswith('gpd_solari_')]
    
ds={}
sites = []
lon = []
lat = []
for f,file in enumerate(ats_files): #loop over each location
    print(file)
    
    new_ds = xr.open_dataset(os.path.join(ats_output_dir,file)) #open ats output for current site
    
    #drop and rename some variables
    new_ds = new_ds.drop_vars(['dec_ex','est','lower','upper','boot','xi_ci','sigma_ci'])
    new_ds = new_ds.drop_dims(['return_period','id_timeseries','quantiles'])
    
    new_ds = new_ds.rename_vars({'xi':'shape','sigma':'scale','thres':'location',
                                 'rate':'avg_exceed','prctile':'threshold_pct',
                                 'xi_boot':'shape_samples','sigma_boot':'scale_samples'})
    
    
    data = g3.file_to_pandas(file.replace('gpd_solari_','').replace('.nc',''))
    metadata = data[1]
    
    sites.append(metadata['site_name'])
    lon.append(metadata['longitude'])
    lat.append(metadata['latitude'])
    
    #add data for current site to dataset
    if f==0:
        ds = new_ds
    else:
        ds = xr.concat([ds,new_ds],dim='site')

ds = ds.squeeze(dim='est_d',drop=True)
ds['site'] = sites
ds['lon'] = ds['shape']*0 + lon
ds['lat'] = ds['shape']*0 + lat
ds['nboot'] = np.arange(0,10000)
ds['threshold_pct'] = 100*ds['threshold_pct']
ds.attrs['esl_file'] = os.path.join(ats_output_dir,file)

ds.to_netcdf('/Volumes/Naamloos/PhD_Data/GESLA3/GPD_fits/gesla3_gpd_daily_max_potATS_Solari.nc',mode='w')
