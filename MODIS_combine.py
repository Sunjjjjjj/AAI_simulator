# -*- coding: utf-8 -*-
"""
Intercomparison among satellite data

@author: sunj
"""


import sys, os
import shutil
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset 
import glob
import tables
from scipy import ndimage
from mpl_toolkits.basemap import Basemap
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn.metrics import mean_squared_error
from scipy.optimize import leastsq
from shutil import copyfile
import subprocess 
from scipy import optimize
from scipy.interpolate import griddata
from math import sqrt
#from OMAERO import gridOMAERO, OMAERO4cfg
#from MODIS import gridMODIS, MODIS4cfg
from pyhdf.SD import SD, SDC 
from scipy.misc import bytescale

#plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':2,'font.size':22,'axes.labelsize':22,'axes.titlesize':22,'xtick.labelsize':18,'ytick.labelsize':18,'legend.fontsize':22})
plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':2,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':14})


year = 2017
month = 1
days = list(np.arange(26,30,1))
days.remove(28)
day = 26
aerlat = -33.46
aerlon = -70.66

print 'Calculate time zone'
if abs(aerlon)%15 < 7.5:
    jetlag = abs(aerlon)//15
else:
    jetlag = abs(aerlon)//15+1
    
if aerlon<0:
    jetlag = - jetlag 
print jetlag 


ROI = [-40,-12.5,-97.5,-70]

crival = 2

res = 0.5

    
print '**** Reading MODIS %02i-%02i-%4i' % (day, month, year)
moddir = '/nobackup/users/sunj/MODIS/AQUA/MYD021KM/%4i/%02i/%02i/' % (year, month, day)
coordir = '/nobackup/users/sunj/MODIS/AQUA/MYD03/%4i/%02i/%02i/' % (year, month, day)

filelist = glob.glob( moddir + '*.hdf')
coorlist = glob.glob( coordir + '*.hdf')

red = [] 
modlat = [] 
modlon = [] 

plt.figure(figsize=(4.5,4), frameon = False) 
for i in range(len(filelist)):  
    data = SD(filelist[i], SDC.READ)
    coor = SD(coorlist[i], SDC.READ)
    
    selected_sds = data.select('EV_250_Aggr1km_RefSB')
    selected_sds_attributes = selected_sds.attributes()
    
    for key, value in selected_sds_attributes.iteritems():
        if key == 'reflectance_scales':
    		reflectance_scales_250_Aggr1km_RefSB = np.asarray(value)		
        if key == 'reflectance_offsets':
    		reflectance_offsets_250_Aggr1km_RefSB = np.asarray(value)	
    
    sds_data_250_Aggr1km_RefSB = selected_sds.get()
    


    selected_sds = data.select('EV_500_Aggr1km_RefSB')
    selected_sds_attributes = selected_sds.attributes()
    
    for key, value in selected_sds_attributes.iteritems():
        if key == 'reflectance_scales':
    		reflectance_scales_500_Aggr1km_RefSB = np.asarray(value)	
        if key == 'reflectance_offsets':
    		reflectance_offsets_500_Aggr1km_RefSB = np.asarray(value)	
    
    sds_data_500_Aggr1km_RefSB = selected_sds.get()
    
    
    
    selected_sds = coor.select('Latitude')
    myd03_lat = selected_sds.get()
    
    selected_sds = coor.select('Longitude')
    myd03_long = selected_sds.get()

    data_shape = sds_data_250_Aggr1km_RefSB.shape
    
    along_track = data_shape[1]
    cross_trak = data_shape[2]
    
    
    
    z = np.zeros((along_track, cross_trak,3))
    
    for i in np.arange(along_track):
        for j in np.arange(cross_trak): 
            z[i,j,0] = ( sds_data_250_Aggr1km_RefSB[0,i,j] - \
            reflectance_offsets_250_Aggr1km_RefSB[0] ) * \
            reflectance_scales_250_Aggr1km_RefSB[0] 
    
    for i in np.arange(along_track):
        for j in np.arange(cross_trak): 
            z[i,j,1] = ( sds_data_500_Aggr1km_RefSB[1,i,j] - \
            reflectance_offsets_500_Aggr1km_RefSB[1] ) * \
            reflectance_scales_500_Aggr1km_RefSB[1]  
    
    for i in np.arange(along_track):
        for j in np.arange(cross_trak): 
            z[i,j,2] = ( sds_data_500_Aggr1km_RefSB[0,i,j] - \
            reflectance_offsets_500_Aggr1km_RefSB[0] ) * \
            reflectance_scales_500_Aggr1km_RefSB[0] 
    
    z[ z > 1 ] = 1.0
    z[ z < 0 ] = 0.0   
    
    lat_min = myd03_lat[0,0]
    lat_max = myd03_lat[along_track-1,cross_trak-1]
    
    lat_0 = lat_min + (lat_max - lat_min) / 2.
    
    long_min = min(myd03_long[0,0],myd03_long[along_track-1,cross_trak-1])
    long_max = max(myd03_long[0,0],myd03_long[along_track-1,cross_trak-1])
    
    lon_0 = long_min + (long_max - long_min) / 2.
    
    #----------------------------------------------------------------------------------------#
    # Orthographic Map Projection
    
    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    
    ax.patch.set_facecolor((0.75,0.75,0.75))
    
    m1 = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
    
    xpt0, ypt0 = m1(lon_0,lat_0) 
    
    xpt1, ypt1 = m1(myd03_long[0,0],myd03_lat[0,0]) 
    xpt2, ypt2 = m1(myd03_long[0,cross_trak-1],myd03_lat[0,cross_trak-1]) 
    xpt3, ypt3 = m1(myd03_long[along_track-1,cross_trak-1], \
                    myd03_lat[along_track-1,cross_trak-1])
    xpt4, ypt4 = m1(myd03_long[along_track-1,0],myd03_lat[along_track-1,0])
    
    llx = min(xpt1,xpt2,xpt3,xpt4) - xpt0  # lower left
    lly = min(ypt1,ypt2,ypt3,ypt4) - ypt0
    
    urx = max(xpt1,xpt2,xpt3,xpt4) - xpt0  # upper right
    ury = max(ypt1,ypt2,ypt3,ypt4) - ypt0
    
    m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',\
        llcrnrx=llx,llcrnry=lly,urcrnrx=urx,urcrnry=ury)
    
    x_igrid, y_igrid = m(myd03_long,myd03_lat)
    
    x_igrid = x_igrid - xpt0
    y_igrid = y_igrid - ypt0
    
    z_igrid_01 = np.zeros((along_track, cross_trak))
    z_igrid_02 = np.zeros((along_track, cross_trak))
    z_igrid_03 = np.zeros((along_track, cross_trak))
    
    for i in np.arange(2030):
        for j in np.arange(1354): 
            z_igrid_01[i,j] = z[i,j,0]
            z_igrid_02[i,j] = z[i,j,1]
            z_igrid_03[i,j] = z[i,j,2]
    
    x1_igrid = x_igrid.ravel()
    y1_igrid = y_igrid.ravel()
    z_igrid_01 = z_igrid_01.ravel()
    z_igrid_02 = z_igrid_02.ravel()
    z_igrid_03 = z_igrid_03.ravel()
    
    xy1_igrid = np.vstack((x1_igrid, y1_igrid)).T
    xi, yi = np.mgrid[llx:urx:1000j, lly:ury:1000j]
    
    z_01 = griddata(xy1_igrid, z_igrid_01, (xi, yi), method='cubic')
    z_02 = griddata(xy1_igrid, z_igrid_02, (xi, yi), method='cubic')
    z_03 = griddata(xy1_igrid, z_igrid_03, (xi, yi), method='cubic')
    
    rgb_projected = np.zeros((1000, 1000,3))
    for i in np.arange(1000):
        for j in np.arange(1000): 
            rgb_projected[i,j,0] = z_01[i,j]
            rgb_projected[i,j,1] = z_02[i,j]
            rgb_projected[i,j,2] = z_03[i,j]
    
    #rgb_projected[ z > 1 ] = 1.0
    #rgb_projected[ z < 0 ] = 0.0
    whereAreNaNs = np.isnan(rgb_projected);
    rgb_projected[whereAreNaNs] = 0.75;
    
    img = m.imshow(np.rot90(np.fliplr(rgb_projected)), origin='lower')
    
    m.drawcoastlines()
    
    m.drawparallels(np.arange(-90.,120.,5.), color='k', labels=[True,False,False,False])
    m.drawmeridians(np.arange(0.,420.,5.), color='k', labels=[False,False,False,True])
    
    ax.set_xlabel("", fontsize=10)
    ax.set_ylabel("", fontsize=10)       
