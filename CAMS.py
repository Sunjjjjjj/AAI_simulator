# -*- coding: utf-8 -*-
"""
Intercomparison among satellite data

@author: sunj
"""


######### load python packages
import sys, os
import shutil
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 
import glob
import tables
import math
from scipy import ndimage
from mpl_toolkits.basemap import Basemap
from scipy import spatial
from scipy import interpolate
from scipy import ndimage
import matplotlib.mlab as mlab
import matplotlib.patches as patches
from scipy.interpolate import griddata






def gridCAMS(year, month, day, ROI, res): 
    
    print '**** Reading CAMS %02i-%02i-%04i' %(day,month,year)      
    path = '/nobackup/users/sunj/CAMS/nrt_surface/%4i/%02i/%02i/' %(year, month, day)
    filelist = glob.glob( path + '*.nc')
    
    data = netCDF4.Dataset(filelist[0],'r')
    
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]-360.    
    aod550 = data.variables['aod550'][0,:,:]    # total AOD
    ssaod550 = data.variables['ssaod550'][:]    # sea salt AOD 
    duaod550 = data.variables['duaod550'][:]    # dust AOD
    omaod550 = data.variables['omaod550'][:]    # organic matter AOD 
    bcaod550 = data.variables['bcaod550'][:]    # black carbon AOD 
    suaod550 = data.variables['suaod550'][:]    #  sulphate AOD 
    
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
    
    
    lon,lat = np.meshgrid(lon,lat)
    lonarr = np.asarray(lon.reshape(-1))
    latarr = np.asarray(lat.reshape(-1))
    aod550arr = np.asarray(aod550.reshape(-1))
    x,y = np.meshgrid(lonnew,latnew)
    aod550new = griddata((lonarr, latarr), aod550arr, (x, y), method = 'linear')    
    
    return latnew, lonnew, aod550new 




def daytsCAMS(year, month, days, aerlat, aerlon, jetlag): 

    aodts = []
    ts = [] 
    

    
    for i, iday in enumerate(days):
    
        path = '/nobackup/users/sunj/CAMS/nrt_surface/%4i/%02i/%02i/' %(year, month, iday)
        filelist = glob.glob( path + '*.nc')
        
        data = netCDF4.Dataset(filelist[0],'r')
        
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]-360.    
#        time = data.variables['time'][:]
        
        aod550 = data.variables['aod550'][0,:,:]    # total AOD
        ssaod550 = data.variables['ssaod550'][:]    # sea salt AOD 
        duaod550 = data.variables['duaod550'][:]    # dust AOD
        omaod550 = data.variables['omaod550'][:]    # organic matter AOD 
        bcaod550 = data.variables['bcaod550'][:]    # black carbon AOD 
        suaod550 = data.variables['suaod550'][:]    #  sulphate AOD 


        lon,lat = np.meshgrid(lon,lat)
        lonarr = np.asarray(lon.reshape(-1))
        latarr = np.asarray(lat.reshape(-1))
        aod550arr = np.asarray(aod550.reshape(-1))

#        latgol = latgol + list(np.asarray(lat).reshape(-1))
#        longol = longol + list(np.asarray(lon).reshape(-1))
#        aodgol = aodgol + list(np.asarray(aod550).reshape(-1))
#        tsgol = tsgol +list(np.asarray([iday+0.5]).reshape(-1))
#    
#    latarr = np.array(latgol)
#    lonarr = np.array(longol)    
#    aodarr = np.array(aodgol) 
#    tsarr = np.array(tsgol)
        
        
        tree = spatial.KDTree(zip(latarr, lonarr))
        locs = tree.query([aerlat,aerlon], k = 1)
    
    
        aodts.append(aod550arr[locs[1]])
        ts.append(iday+(12+jetlag)/24.)    
        
        data.close()
    
    return ts, aodts
    
    
    

def CAMS4cgf(year, month, day, ROI, plumemsk, jetlag, res): 
    
    print '**** Reading CAMS for Config.in %02i-%02i-%04i' %(day,month,year)    
    path = '/nobackup/users/sunj/CAMS/nrt_surface/%4i/%02i/%02i/' %(year, month, day)
    filelist = glob.glob( path + '*.nc')
    
    data = netCDF4.Dataset(filelist[0],'r')
    
    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]-360.    
    aod550 = data.variables['aod550'][0,:,:]    # total AOD
    ssaod550 = data.variables['ssaod550'][:]    # sea salt AOD 
    duaod550 = data.variables['duaod550'][:]    # dust AOD
    omaod550 = data.variables['omaod550'][:]    # organic matter AOD 
    bcaod550 = data.variables['bcaod550'][:]    # black carbon AOD 
    suaod550 = data.variables['suaod550'][:]    #  sulphate AOD 
    
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
    

    
    lon,lat = np.meshgrid(lon,lat)
    lonarr = np.asarray(lon.reshape(-1))
    latarr = np.asarray(lat.reshape(-1))
    aodarr = np.asarray(aod550.reshape(-1))
    x,y = np.meshgrid(lonnew,latnew)
    aodnew = griddata((lonarr, latarr), aodarr, (x, y), method = 'linear')    


#    # applying MODIS mask 
#    latomi, lonomi, aodomi0 = gridOMAERO(year, month, iday, ROI, jetlag, 'AOD')    
#    latmod, lonmod, aodmod0 = gridMODIS(year, month, day, ROI)
#    
#    mask = np.logical_or(np.isnan(aodomi0), np.isnan(aodmod0))
#    mask = np.logical_or(mask, aodmod0 < crival)
#    mask = np.logical_or(mask, aodnew < crival)
#    mask = np.logical_or(mask, aodmod0 < 0.)

    latm = np.ma.masked_array(y,plumemsk)  
    lonm = np.ma.masked_array(x,plumemsk) 
    aodm = np.ma.masked_array(aodnew,plumemsk)      

    latrm = np.array(latm[latm.mask == 0])    
    lonrm = np.array(lonm[lonm.mask == 0])    
    aodrm = np.array(aodm[aodm.mask == 0])

        
    
    
    return latrm, lonrm, aodrm     