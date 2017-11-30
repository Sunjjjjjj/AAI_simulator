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
import seaborn as sns
from scipy import ndimage
from mpl_toolkits.basemap import Basemap
from scipy import spatial
from scipy import interpolate
from scipy import ndimage
import matplotlib.mlab as mlab
import matplotlib.patches as patches
from scipy.interpolate import griddata
#from CAMS import gridCAMS, daytsCAMS
from datetime import datetime





year = 2017
month = 1
days = list(np.arange(26,32,1))
days.remove(28)
days = [26]
aerlat = -33.46
aerlon = -70.66
time0 = '1900-01-01 00:00:0.0'


ROI = [-40, -5, -105, -70]


#for i, iday in enumerate(days):
    
#    latc, lonc, aod550c = gridCAMS(year, month, iday, ROI)


path = '/nobackup/users/sunj/CAMS/nrt_surface/' 
#    path = '/usr/people/sunj/' 
filelist = glob.glob( path + '*.nc')

data = netCDF4.Dataset(filelist[0],'r')

lat = data.variables['latitude'][:]
lon = data.variables['longitude'][:]-360.
pres = data.variables['level'][:]
time = data.variables['time'][:]

OMd = data.variables['aermr07'][:]
OMw = data.variables['aermr08'][:]
BCd = data.variables['aermr09'][:]
BCw = data.variables['aermr10'][:]

data.close()

plt.figure(figsize = (4,4))
plt.plot(OMd[4:8,:,:,:].mean(3).mean(2).mean(0),pres,label = 'Organic matter hydropho')
plt.plot(OMw[4:8,:,:,:].mean(3).mean(2).mean(0),pres,label = 'Organic matter hydrophi')
plt.plot(BCd[4:8,:,:,:].mean(3).mean(2).mean(0),pres,label = 'Black carbon hydropho')
plt.plot(BCw[4:8,:,:,:].mean(3).mean(2).mean(0),pres,label = 'Black carbon hydrophi')
plt.plot(OMd[4:8,:,:,:].mean(3).mean(2).mean(0) + OMw[4:8,:,:,:].mean(3).mean(2).mean(0) + \
         BCd[4:8,:,:,:].mean(3).mean(2).mean(0) + BCw[4:8,:,:,:].mean(3).mean(2).mean(0) ,pres,'k.-',label = 'OM + BC')
plt.xlabel('Mixing ratio [kg/kg]')
plt.ylabel('Pressure [hPa]')    
plt.gca().invert_yaxis()
plt.legend(frameon = False)


#    aod550 = data.variables['aod550'][:]    # total AOD
#    omaod550 = data.variables['omaod550'][:]    # organic matter AOD 
#    bcaod550 = data.variables['bcaod550'][:]    # black carbon AOD 
    
#    for j in np.arange(0,len(time)): 
#        print 'Date:',15+j/2
#        print aod550[j,:,:].mean() 
#        print omaod550[j,:,:].mean() 
#        print bcaod550[j,:,:].mean() 

#    ssaod550 = data.variables['ssaod550'][:]    # sea salt AOD 
#    duaod550 = data.variables['duaod550'][:]    # dust AOD
#    omaod550 = data.variables['omaod550'][:]    # organic matter AOD 
#    bcaod550 = data.variables['bcaod550'][:]    # black carbon AOD 
#    suaod550 = data.variables['suaod550'][:]    #  sulphate AOD 
#    
#    latnew = np.arange(ROI[0], ROI[1], 0.1)
#    lonnew = np.arange(ROI[2], ROI[3], 0.1)
#    
#    
#    lon,lat = np.meshgrid(lon,lat)
#    lonarr = np.asarray(lon.reshape(-1))
#    latarr = np.asarray(lat.reshape(-1))
#    aod550arr = np.asarray(aod550.reshape(-1))
#    x,y = np.meshgrid(lonnew,latnew)
#    aod550new = griddata((lonarr, latarr), aod550arr, (x, y), method = 'linear')    
#
#    fig = plt.figure(figsize=(4.5,4))
#    cbax = fig.add_axes([0.85,0.1,0.05,0.8])  
#    ax = fig.add_axes([0.0,0.1,0.8,0.8])
#    map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], lat_0 = 0, lon_0 = 0, projection='cyl',resolution='c')
#    map.drawcoastlines(color='black',linewidth=0.6)
#    map.drawcountries(color = 'black',linewidth=0.4)
#    map.drawmeridians(np.arange(-115,-45,30),labels = [0,1,1,0])
#    map.drawparallels(np.arange(-65,10,30), labels = [0,1,1,0])
#    cb = plt.pcolor(lonc,latc,aod550c,cmap='pink_r', vmin=0, vmax=2)   
#    
#    ttl = plt.title('%4i%02i%02i' %(year, month, iday)) 
#    ttl.set_position([0.2,0.9])
#    plt.xlabel('Longitude')
#    plt.ylabel('Latitude')
#    cbar = plt.colorbar(cb,cax = cbax, ticks = np.arange(0,2.1,1))
#    cbar.set_label('AOD', rotation =270)
#
#    plt.savefig('../Figures/CAMS_AOD_%02i.png' % (iday), dpi = 300)


