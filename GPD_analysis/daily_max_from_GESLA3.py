#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" daily_max_from_GESLA3.py

Extracts daily maxima from GESLA3 raw data and adds them to a dictionary of multiple sites data. 
For sites with the same name or closer than 3 km together, the site with the longest record is used.

input:
    g3object: GESLA3 object used to read data from
    sites_data: dictionary to fill with data from sites
    g3_file: name of GESLA3 file to use
    min_num_years: required minimum number/365.25 of daily maxima at a site
    gauge_type: 'Coastal, 'River', 'Lake','Any'
    detrend: detrend daily maxima or not
    deseasonalize: deseasonalize daily maxima or not
    subtract_annualmeans: subtract annual means from daily maxima or not
    
output:
    sites_data: dictionary appended with daily maxima of g3file
    
Created on Fri Mar 4 14:26 2022
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""
import numpy as np
import pandas as pd
import gesla
from compute_distances_to_point import compute_distances_to_point

def daily_max_from_GESLA3(g3object,sites_data,g3file,gauge_type,min_num_years,
                                   detrend=None,deseasonalize=None,subtract_annualmeans=None):
    #load data
    data = g3object.file_to_pandas(g3file)
    rawdata = data[0]
    metadata = data[1]
    
    print('Processing station: ' + metadata.site_name)
    
    #check metadata
    if metadata.overall_record_quality!= 'No obvious issues':
        print('Overall record quality: '+metadata.overall_record_quality)
        return sites_data
    
    if metadata.number_of_years < min_num_years:
        print('Insufficient total years of observations')
        return sites_data
    
    if ((metadata.gauge_type != gauge_type) & (gauge_type != 'Any')):
        print('Wrong gauge type')     
        return sites_data
    
    rawdata = rawdata[rawdata['use_flag']==1] #only use data with use flag 1 (qf=0-2) (see GESLA3 documentation)
    
    #convert data to hourly means
    rawdata['date_time']=rawdata.index
    data_hourly = rawdata.resample('H',on='date_time').mean() #regardless of # obs in hour
    data_hourly = data_hourly[np.isfinite(data_hourly['sea_level'])] #remove hours without any data
    data_hourly['date_time']=data_hourly.index #add date time index as separate column
    
    daily_max = rawdata.resample('D',on='date_time').max() #daily maxima
    
    #only consider days with 12 or more hourly means of obs
    hourly_pday = data_hourly.resample('D',on='date_time').count()
    daily_max = daily_max[hourly_pday['sea_level']>=12]
    
    if len(daily_max) < min_num_years*365.25:
        print('Insufficient daily maxima of observations')
        return sites_data
    
    #if site with same name included, pick the one with the longest record
    if metadata.site_name in sites_data:
        if len(sites_data[metadata.site_name]['daily_max']) > len(daily_max):
            print('Site already included with longer record')
            return sites_data #else, overwritten at end of function, as dictionary key = site name
        
    #if sites within 3 km included, pick the one with the longest record
    lons = np.array([sites_data[k]['metadata']['longitude'] for k in sites_data]) #coordinates of all included sites
    lats = np.array([sites_data[k]['metadata']['latitude'] for k in sites_data])
    
    distances = compute_distances_to_point(metadata.longitude,metadata.latitude,lons,lats) #distances from current site to all included sites

    if np.any(distances < 3): #distances lower than 3 km
        for d,distance in enumerate(distances[distances < 3]):
            i = np.where(distances<3)[0][d] #index in sites_data dictionary
            
            if len(sites_data[list(sites_data)[i]]['daily_max']) > len(daily_max):
                print('Site in close proximity to included site with more observations')
                return sites_data
            else:
                del sites_data[list(sites_data)[i]]
    
    #detrend, deseasonalize and/or subtract annual means:
    if detrend: 
        x =  daily_max.index.values.astype(np.int64) // 10 ** 9 #convert to seconds timestamp
        y =  daily_max['sea_level'].values.astype('float64')
        lrcoefs = np.polyfit(x,y,1)
        trend = np.polyval(lrcoefs,x)
        
        daily_max['sea_level'] = daily_max['sea_level'] - trend + lrcoefs[-1]

    if deseasonalize:
        daily_max['sea_level'] = daily_max['sea_level'] - daily_max.groupby(daily_max.index.month).transform('mean')['sea_level'].astype('float64') + np.mean(daily_max.groupby(daily_max.index.month).mean()['sea_level'])
        
    if subtract_annualmeans:
        daily_max['sea_level'] = daily_max['sea_level'] - daily_max.groupby(daily_max.index.year).transform('mean')['sea_level'].astype('float64')
    
    sites_data[metadata.site_name] = {'daily_max':daily_max['sea_level'],'metadata':metadata} #append current site to dictionary
    return sites_data
