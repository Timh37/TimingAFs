#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" store_daily_maxima_GESLA3.py

Saves daily maxima extracted from GESLA3 observations to .csv files.

input:
    min_num_years: required minimum number/365.25 of daily maxima at a site
    gauge_type: 'Coastal, 'River', 'Lake','Any'
    detrend: detrend daily maxima or not
    deseasonalize: deseasonalize daily maxima or not
    subtract_annualmeans: subtract annual means from daily maxima or not

    pathing
    
output:
    sites_data: dictionary organized by site name containing daily maxima and site metadata 
    
Created on Fri Mar 4 14:26 2022
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""
import numpy as np
import os
import pickle
from datetime import datetime
from gesla import GeslaDataset
from daily_max_from_GESLA3 import daily_max_from_GESLA3

#prepare loading GESLA3 data
g3dir = '/Volumes/Naamloos/PhD_Data/GESLA3/GESLA3.0_ALL/' #gesla data
g3files = os.listdir(g3dir) #find files in GESLA directory
g3 = GeslaDataset(meta_file='/Volumes/Naamloos/PhD_Data/GESLA3/GESLA3_ALL.csv',data_path=g3dir) #load metadata

#set preprocessing steps
min_num_years = 30 
gauge_type = 'Coastal'
detrend = True
deseasonalize = True
subtract_annualmeans = False

#generate output directory name
daily_max_output_dir = os.path.join('/Volumes/Naamloos/PhD_Data/GESLA3/daily_maxima',
                          str(datetime.now().day)+'_'+str(datetime.now().month)+'_'+str(datetime.now().year)+'_'+
                          '_'.join(list(np.array(['detrended','deseason','noAnnualMean'])[np.array([detrend,deseasonalize,subtract_annualmeans])]))
                          +'_'+gauge_type.lower()+'_minyr'+str(min_num_years))
        
sites_data = {} #initialize dictionary

for filename in g3files: #loop over station files to extract daily maxima
    sites_data = daily_max_from_GESLA3(g3,sites_data,filename,gauge_type,min_num_years,
                                                       detrend=detrend,deseasonalize=deseasonalize,subtract_annualmeans=subtract_annualmeans)

#save dictionary to pickle file
open(os.path.join(daily_max_output_dir,'sites_data.txt'), 'a').close() #this makes the file if not already existing
with open(os.path.join(daily_max_output_dir,'sites_data.txt'), 'wb') as myFile:    
          pickle.dump(sites_data,myFile) #store sites_data dictionary to pickle

for s,site_name in enumerate(sites_data): #write daily maxima for each site to csv file
    print('Storing daily maxima for: '+site_name)
    if not os.path.isdir(daily_max_output_dir):
        os.mkdir(daily_max_output_dir)
    site_data = sites_data[site_name]
    
    site_data['daily_max'].to_csv(os.path.join(daily_max_output_dir,'daily_max_'+site_data['metadata']['file_name']+'.csv'),index=True)