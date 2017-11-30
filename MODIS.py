# -*- coding: utf-8 -*-
"""
Read OMAERO data 
- rawOMAERO: return orbit data, organized by list
- gridMODIS: return gridded data over ROI, organized by np.2darray
- daytsMODIS: return time series 


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
from scipy import ndimage
from mpl_toolkits.basemap import Basemap
import pandas as pd
from pandas import Series, DataFrame, Panel
from scipy import spatial
from datetime import datetime
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata


def gridMODIS(year, month, day, ROI, res, data): 
    
    print '**** Reading MODIS %02i-%02i-%4i' % (day, month, year)
    moddir = '/nobackup/users/sunj/MODIS/AQUA/MYD04/%4i/%02i/%02i/' % (year, month, day)

    aodmod = []
    latmod = []
    lonmod = [] 
    
    filelist = glob.glob( moddir + '*' + data + '*.hdf')
#    print filelist
    for io in filelist[:]:  
#        print io
        data = netCDF4.Dataset(io,'r')    

        lat = data.variables['Latitude'][:]
        lon = data.variables['Longitude'][:]
        aod = data.variables['Optical_Depth_Land_And_Ocean'][:]
        sza = data.variables['Solar_Zenith'][:]
        vza = data.variables['Sensor_Zenith'][:]
        saa = data.variables['Solar_Azimuth'][:]
        vaa = data.variables['Sensor_Azimuth'][:]
        sca = data.variables['Scattering_Angle'][:]
        
        raa = saa + 180 - vaa 
        
        scacal = np.arccos(-np.cos(sza/180.*np.pi)*np.cos(vza/180.*np.pi)+np.sin(sza/180.*np.pi)*np.sin(vza/180.*np.pi)*np.cos(raa/180.*np.pi)) / np.pi * 180
#        print 'scacal - sca :', scacal - sca
      
        
        latmod = latmod + list(np.asarray(lat).reshape(-1))
        lonmod = lonmod + list(np.asarray(lon).reshape(-1))        
        aodmod = aodmod + list(np.asarray(aod).reshape(-1))
        
        data.close() 
    
    latarr = np.array(latmod)
    lonarr = np.array(lonmod)   
    aodarr = np.array(aodmod) 
  

    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
    
    x,y = np.meshgrid(lonnew,latnew)
    aodnew = griddata((lonarr, latarr), aodarr , (x, y), method = 'linear')
    aodnew[np.isnan(aodnew)] = -32767.

        
    return latnew, lonnew, aodnew
    
    
    
    

def daytsMODIS(year, month, day, aerlat, aerlon, jetlag):
    aodmod = [] 
    tsmod = [] 
    for iday in range(19,32):
        print '**** Reading MODIS %02i-%02i-%4i' % (iday, month, year)
        moddir = '/nobackup/users/sunj/MODIS/AQUA/MYD04/%4i/%02i/%02i/' % (year, month, iday)
        modfile = glob.glob( moddir + '*.txt')[0]
        moddata = pd.read_csv(modfile, delim_whitespace=True, header = 7, \
                            names = ['date','time','longitude','latitude','AOD','Angstrom',\
                            'fine_mode_frac','sza','vza','scattering_angle','sunglint_angle'])
        lat = np.array(moddata['latitude'])
        lon = np.array(moddata['longitude'])
        ts = np.array(moddata['time'])
        aod = np.array(moddata['AOD'])
        
        tree = spatial.KDTree(zip(lat, lon))
        locs = tree.query([aerlat,aerlon], k = 1)
        
        aodmod.append(aod[locs[1]])
        
        
        temp = datetime.strptime(ts[locs[1]],'%H:%M:%S.00')
        tsjulian = iday + ( temp.hour + temp.minute/60. + temp.second/3600. + jetlag) / 24. 
        
        tsmod.append(tsjulian)
    
    return tsmod, aodmod 

    



def MODIS4cfg(year, month, day, ROI, plumemsk, res,data): 
    
    print '**** Reading MODIS %02i-%02i-%4i' % (day, month, year)
    moddir = '/nobackup/users/sunj/MODIS/AQUA/MYD04/%4i/%02i/%02i/' % (year, month, day)

    aodmod = []
    latmod = []
    lonmod = [] 
    angmod = []
    asymod = []
    ssamod = []     
    
    filelist = glob.glob( moddir + '*' + data + '*.hdf')
    for io in filelist[:]:  
#        print io
        data = netCDF4.Dataset(io,'r')    

        lat = data.variables['Latitude'][:]
        lon = data.variables['Longitude'][:]
        aod = data.variables['Optical_Depth_Land_And_Ocean'][:]
        ang = data.variables['Angstrom_Exponent_1_Ocean'][:]
        asy = data.variables['Asymmetry_Factor_Best_Ocean'][1,:,:]              # @ 550 nm 
#        ssa = data.variables['Deep_Blue_Spectral_Single_Scattering_Albedo_Land'][:].mean(axis=0)    
        sza = data.variables['Solar_Zenith'][:]
        vza = data.variables['Sensor_Zenith'][:]
        saa = data.variables['Solar_Azimuth'][:]
        vaa = data.variables['Sensor_Azimuth'][:]
    
        latmod = latmod + list(np.asarray(lat).reshape(-1))
        lonmod = lonmod + list(np.asarray(lon).reshape(-1))        
        aodmod = aodmod + list(np.asarray(aod).reshape(-1))
        angmod = angmod + list(np.asarray(ang).reshape(-1))
        asymod = asymod + list(np.asarray(asy).reshape(-1))
#        ssamod = ssamod + list(np.asarray(ssa).reshape(-1))
        
        data.close()
        
    
    latarr = np.array(latmod)
    lonarr = np.array(lonmod)   
    aodarr = np.array(aodmod) 
    angarr = np.array(angmod) 
    asyarr = np.array(asymod)
#    ssaarr = np.array(ssamod)
  
  
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
    
    x,y = np.meshgrid(lonnew,latnew)
    aodnew = griddata((lonarr, latarr), aodarr , (x, y), method = 'linear')
    angnew = griddata((lonarr, latarr), angarr , (x, y), method = 'linear')
    asynew = griddata((lonarr, latarr), asyarr , (x, y), method = 'linear')
#    ssanew = griddata((lonarr, latarr), ssaarr , (x, y), method = 'linear')
    
    
    latm = np.ma.masked_array(y,plumemsk)  
    lonm = np.ma.masked_array(x,plumemsk) 
    aodm = np.ma.masked_array(aodnew,plumemsk)    
    angm = np.ma.masked_array(angnew,plumemsk)    
    asym = np.ma.masked_array(asynew,plumemsk)    
#    ssam = np.ma.masked_array(ssanew,plumemsk)    

    latrm = np.array(latm[latm.mask == 0])    
    lonrm = np.array(lonm[lonm.mask == 0])  
    aodrm = np.array(aodm[aodm.mask == 0])
    angrm = np.array(angm[angm.mask == 0])
    asyrm = np.array(asym[asym.mask == 0])
#    ssarm = np.array(ssam[ssam.mask == 0])

        
    return latrm, lonrm, aodrm, angrm, asyrm