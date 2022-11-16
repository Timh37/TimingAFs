#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" 
Plot central estimate and histograms of timing of flopros AF10 and 100 for SSP126 and SSP370
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import xarray as xr
import os
import cartopy.crs as ccrs
import matplotlib
import cartopy
import cmocean
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import PercentFormatter
from matplotlib import colors as c

plt.close('all')

in_dir = '/Volumes/Naamloos/PhD_Data/GESLA3/timing_projections' #path to timing projections

#load timing projections
esl_ds_flopros = xr.open_dataset(os.path.join(in_dir,'esl_AF_timing_ar6_2150wfs_yrRef2022_fRefFLOPROS_muATS_extrapSweet.nc'))
esl_ds_flopros = esl_ds_flopros.drop_isel(site=np.where(np.isnan(esl_ds_flopros.flopros))[0])
esl_ds_flopros = esl_ds_flopros.drop_isel(site=np.where(np.isnan(esl_ds_flopros.projected_slr_pbox.sel(year=2150).sel(ssp='ssp370',pct=50)))[0]) #for 3 TGs there are no projections

#sample frequencies matrix
freq_mat = np.matlib.repmat(esl_ds_flopros.freq.values,len(esl_ds_flopros.site),1)

#figure
fig=plt.figure(figsize=(7,8.8))   
gs = fig.add_gridspec(17,2)
gs.update(left=.07,wspace=.55,right=.935,top=.96,hspace=.7,bottom=.06)

#make discrete colormap
cmdict = cmocean.tools.get_dict(cmocean.cm.amp_r)
cmap=matplotlib.colors.LinearSegmentedColormap('amp_discrete',cmdict,16).reversed()

#configure map insets for densely covered regions
map_insets_bbox = [(.2, -.3,.7,.53),(.65, -.3,.7,.53),(-.25, -.3,.7,.53)]
map_insets_coords = [[-11,29,40,65],[104,144,15,40],[-99,-59,25,50]]

#set up subplots
subplot_ids = np.array([[0,9],[4,13]])  
subplot_labels = np.array([['ae','bf'],['cg','dh']]) 
 
for iaf,af in enumerate([10,100]): #loop over AFs
    for sp,ssp in enumerate(['ssp126','ssp370']): #loop over scenarios
        
        #locate and fetch sample frequency equal to af*flopros frequency
        fut_flopros_mat = np.transpose(np.matlib.repmat(esl_ds_flopros.flopros.values*af,len(esl_ds_flopros.freq),1))
        idx = np.argmin(np.abs(fut_flopros_mat-freq_mat),axis=1)
        
        #plot central estimate
        timing = esl_ds_flopros.timing_ce.sel(ssp=ssp).values[tuple(np.array([np.arange(0,len(esl_ds_flopros.site)),idx]))]
        
        isRising = np.any(esl_ds_flopros.timing_ce.sel(ssp=ssp).values>2022,axis=1) #distinguish between relative sea-level rise and fall
        
        #make map
        ax = plt.subplot(gs[0+9*sp:3+9*sp,iaf], projection=ccrs.Robinson(central_longitude=0)) #maps
        ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='darkgrey')
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='dimgrey')
        ax.scatter(esl_ds_flopros.lon[((np.isfinite(timing))&(isRising))],esl_ds_flopros.lat[((np.isfinite(timing))&(isRising))],c=timing[((np.isfinite(timing))&(isRising))],cmap=cmap,vmin=2022,vmax=2150,s=7,transform=ccrs.PlateCarree(),zorder=3)
        ax.scatter(esl_ds_flopros.lon[((~np.isfinite(timing))&(isRising))],esl_ds_flopros.lat[((~np.isfinite(timing))&(isRising))],c='None',edgecolor='cyan',s=7,transform=ccrs.PlateCarree(),zorder=1)
        ax.scatter(esl_ds_flopros.lon[~isRising],esl_ds_flopros.lat[~isRising],c='None',edgecolor='yellow',s=7,transform=ccrs.PlateCarree(),zorder=2,marker='D')
        
        ax.set_title('('+subplot_labels[0,iaf][sp]+') AF$_{FLOPROS}$ = '+str(af))
        ax.set_extent([-170, 170, -75, 75], crs=ccrs.PlateCarree())
        
        #make insets
        for i in np.arange(3):
            axins = inset_axes(ax,width="100%", height="100%",bbox_to_anchor=map_insets_bbox[i],bbox_transform=ax.transAxes,#width=2, height=2, #bbox_to_anchor=(0.9, 0.9),#loc="upper right", 
                               axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                               axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree()))
            axins.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='darkgrey')
            axins.add_feature(cartopy.feature.LAND, zorder=0,facecolor='dimgrey')
            sc=axins.scatter(esl_ds_flopros.lon[((np.isfinite(timing))&(isRising))],esl_ds_flopros.lat[((np.isfinite(timing))&(isRising))],c=timing[((np.isfinite(timing))&(isRising))],cmap=cmap,vmin=2022,vmax=2150,s=7,transform=ccrs.PlateCarree(),zorder=3)
            axins.scatter(esl_ds_flopros.lon[((~np.isfinite(timing))&(isRising))],esl_ds_flopros.lat[((~np.isfinite(timing))&(isRising))],c='None',edgecolor='cyan',s=7,transform=ccrs.PlateCarree(),zorder=1)
            axins.scatter(esl_ds_flopros.lon[~isRising],esl_ds_flopros.lat[~isRising],c='None',edgecolor='yellow',s=7,transform=ccrs.PlateCarree(),zorder=2,marker='D')
            axins.set_extent(map_insets_coords[i],crs=ccrs.PlateCarree())
        
        #draw colorbar
        if iaf==0:
            if sp==1:
                cax=inset_axes(ax,width="100%", height="100%",bbox_to_anchor=(0, -2.32,.78*1.05,.18),bbox_transform=ax.transAxes)
                cb=fig.colorbar(sc, cax=cax,orientation='horizontal',extend='max')
                ax.annotate('Timing [yr]',xy=(.345,.034),xycoords='figure fraction',rotation=0)
                
                cax=inset_axes(ax,width="100%", height="100%",bbox_to_anchor=(1.55, -2.32,.78*1.05,.18),bbox_transform=ax.transAxes)
                cb=fig.colorbar(sc, cax=cax,orientation='horizontal',extend='max')
                ax.annotate('Timing [yr]',xy=(.345+1*.525,.034),xycoords='figure fraction',rotation=0)
        
        
        #make cumulative mesh plots
        ax = plt.subplot(gs[4+9*sp:7+9*sp,iaf])
        
        timing = esl_ds_flopros.timing_pbox.sel(ssp=ssp).values.transpose(0,2,1)
        timing_nfold= np.zeros((timing.shape[0],timing.shape[2])) #initialize timing of AFs at each location
        isRising = np.any(timing>2022,axis=1) #determine whether relative sea-level rise or fall
        
        qntsNotRising = np.nan*timing_nfold #becomes yellow
        laterThan2150 = np.nan*timing_nfold #becomes cyab
        
        for s,site  in enumerate(esl_ds_flopros.site.values): #loop over sites, could be vectorized for speed  
            timing_nfold[s,:] = timing[(s,idx[s])]  #get timing of flopros*af for current site
        
        for pc,pct in enumerate(esl_ds_flopros.pct.values): #for each output percentile, determine  
            qntsNotRising[np.sum(isRising,axis=0)[pc]::,pc]=1 #if quantiles have relative sea-level fall
            laterThan2150[np.sum(np.isfinite(timing_nfold),axis=0)[pc]::,pc]=1 #if quantiles have less relative sea-level rise than required
            
        timing_nfold_sorted = np.sort(timing_nfold,axis=0) #sort timing across sites to get histogram for each percentile
        p = 1. * np.arange(len(esl_ds_flopros.site)) / (len(esl_ds_flopros.site) - 1) #cumulative fraction of TG sites
        
        #do the plotting
        pm = ax.pcolormesh(np.arange(0,100,2.5)[0:20],p,timing_nfold_sorted[:,0:20],cmap=cmap,vmin=2022,vmax=2150)
        pm = ax.pcolormesh(np.arange(0,100,2.5)[21::],p,timing_nfold_sorted[:,21::],cmap=cmap,vmin=2022,vmax=2150)
        ax.pcolormesh(np.arange(0,100,2.5)[21::],p,laterThan2150[:,21::],cmap=c.ListedColormap(['cyan']),vmin=2022,vmax=2150,zorder=0)
        ax.pcolormesh(np.arange(0,100,2.5)[0:20],p,laterThan2150[:,0:20],cmap=c.ListedColormap(['cyan']),vmin=2022,vmax=2150,zorder=0)
        ax.pcolormesh(np.arange(0,100,2.5)[21::],p,qntsNotRising[:,21::],cmap=c.ListedColormap(['yellow']),vmin=2022,vmax=2150,zorder=1)
        ax.pcolormesh(np.arange(0,100,2.5)[0:20],p,qntsNotRising[:,0:20],cmap=c.ListedColormap(['yellow']),vmin=2022,vmax=2150,zorder=1)
        ax.plot([48.75,48.75],[0,1],c='black',linewidth=.75)
        ax.plot([51.25,51.25],[0,1],c='black',linewidth=.75)
        
        cs = ax.contour(np.arange(0,100,2.5)[0:20],p,timing_nfold_sorted[:,0:20],levels=[2052,2072],colors=['white','white'],linestyles=['dashed','dashed'])
        cs = ax.contour(np.arange(0,100,2.5)[21::],p,timing_nfold_sorted[:,21::],levels=[2052,2072],colors=['white','white'],linestyles=['dashed','dashed'])
        
        if af==100:
            if sp==0:
                ax.text(12,.5,'2072',color='white')
                ax.text(28,.1,'2052',color='white')
            else:
                ax.text(17,.55,'2072',color='white')
                ax.text(28,.16,'2052',color='white')
        else:
            if sp==0:
                ax.text(12,.65,'2072',color='white')
                ax.text(28,.29,'2052',color='white')
            else:
                ax.text(17,.69,'2072',color='white')
                ax.text(28,.31,'2052',color='white')
                
        ax.set_xlabel('Percentile at each TG [-]')
        ax.set_ylabel('Fraction of TGs [-]')
        ax.set_facecolor('white')
        ax.grid(color='darkgrey',linestyle='--', linewidth=.75)
        ax.set_xticks([5,20,35,50,65,80,95])
        ax.set_yticks([0,.2,.4,.6,.8,1])
        ax.set_title('('+subplot_labels[1,iaf][sp]+') AF$_{FLOPROS}$ = '+str(af))
        ax.set_ylim([0,1])
        ax.set_xlim([2.5,97.5])
        
        ax.yaxis.set_major_formatter(PercentFormatter(1))

        #add central-estimate cdf on the right hand side of the plot
        a = np.sort(esl_ds_flopros.timing_ce.sel(ssp=ssp).values[tuple(np.array([np.arange(0,len(esl_ds_flopros.site)),idx]))])
        isRising = np.any(esl_ds_flopros.timing_ce.sel(ssp=ssp).values>2022,axis=1)
        
        b = np.nan*a #fraction with relative sea-level fall
        b[-np.sum(~isRising)::]=1
        
        axins = inset_axes(ax,width="100%", height="100%",bbox_to_anchor=(1.065, 0.059,.05,1),bbox_transform=ax.transAxes)
        pm = axins.pcolormesh([1,2],p,np.transpose(np.matlib.repmat(a,2,1)),cmap=cmap,vmin=2022,vmax=2150)
        axins.pcolormesh([1,2],p,np.transpose(np.matlib.repmat(b,2,1)),cmap=c.ListedColormap(['yellow']))
        cs = axins.contour([1,2],p,np.transpose(np.matlib.repmat(a,2,1)),levels=[2052,2072],colors=['white','white'],linestyles=['dashed','dashed'])

        axins.set_facecolor('cyan') #fraction with relative sea-level rise less than required
        ax.set_ylim([0,1])
        axins.set_xlim([1,1.5])
        axins.set_xticks([1.25],labels=['ce'])
        axins.set_yticks([0,.2,.4,.6,.8,1],labels=[])
        axins.grid(color='darkgrey',linestyle='--', linewidth=.75)

        for tick in axins.yaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
        
        

ax.annotate('SSP1-2.6',xy=(.02,.985),xycoords='figure fraction',rotation=0,color='blue',fontweight='bold')
ax.annotate('SSP3-7.0',xy=(.02,.503),xycoords='figure fraction',rotation=0,color='orange',fontweight='bold')

#plt.savefig('/Users/timhermans/Documents/PhD/Phase5_timing_esl/Paper/Figures/Fig4_timingPbox_flopros_AF10_100_v9.png',format='png',dpi=300)
