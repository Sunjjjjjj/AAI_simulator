# -*- coding: utf-8 -*-
"""
Read Metop A&B GOME2 data 
- rawGOME2: return Metop A&B data, organized by list
- gridGOME2: return gridded data over ROI, organized by np.2darray

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
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata


def rawGOME2(year, month, day): 

    print '**** Reading GOME2 %02i-%02i-%04i' %(day,month,year) 
    
    pathA = '/nobackup/users/sunj/GOME-2/MetopA/%4i/%02i/%02i/' %(year, month, day)
    filelistA = glob.glob( pathA + '*.hdf5')
    pathB = '/nobackup/users/sunj/GOME-2/MetopB/%4i/%02i/%02i/' %(year, month, day)
    filelistB = glob.glob( pathB + '*.hdf5')


    latA = []
    lonA = []
    aaiA = []     
    for io in filelistA[:]:  

        dataA = netCDF4.Dataset(io,'r')    
        data = dataA.groups['DATA']
        geo = dataA.groups['GEOLOCATION']

        lat = geo.variables['LatitudeCenter'][:]
        lon = geo.variables['LongitudeCenter'][:]
        aai = data.variables['AAI'][:]
        sza = geo.variables['SolarZenithAngle'][:]
        sunglint = data.variables['SunGlintFlag'][:]
        
        
        mask = np.logical_or(sza[0:-1,0:-1] < 0. , sza[0:-1,0:-1] > 75.)  
        mask = np.logical_or(mask, sunglint[0:-1,0:-1] < 0.)                # sunglint pixel
        mask = np.logical_or(mask, aai[0:-1,0:-1] < 0.)                     # scattering aerosol pixel  
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel 
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
        
        
        if mask.all() != True: 
            lon = np.ma.masked_array(lon[0:-1,0:-1],mask)
            lat = np.ma.masked_array(lat[0:-1,0:-1],mask)
            aai = np.ma.masked_array(aai[0:-1,0:-1],mask)  
            
            latA.append(lat)  
            lonA.append(lon)
            aaiA.append(aai)
            

    latB = []
    lonB = [] 
    aaiB = [] 
    for io in filelistB[:]: 

        dataB = netCDF4.Dataset(io,'r')
        data = dataB.groups['DATA']
        geo = dataB.groups['GEOLOCATION']


        lat = geo.variables['LatitudeCenter'][:]
        lon = geo.variables['LongitudeCenter'][:]    
        aai = data.variables['AAI'][:]
        sza = geo.variables['SolarZenithAngle'][:]
        sunglint = data.variables['SunGlintFlag'][:]
        
        
        mask = np.logical_or(sza[0:-1,0:-1] < 0. , sza[0:-1,0:-1] > 75.)  
        mask = np.logical_or(mask, sunglint[0:-1,0:-1] < 0.)                # sunglint pixel
        mask = np.logical_or(mask, aai[0:-1,0:-1] < 0.)                     # scattering aerosol pixel  
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel 
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
        
        
        if mask.all() != True: 
            lon = np.ma.masked_array(lon[0:-1,0:-1],mask)
            lat = np.ma.masked_array(lat[0:-1,0:-1],mask)
            aai = np.ma.masked_array(aai[0:-1,0:-1],mask)  
            
            latB.append(lat)
            lonB.append(lon)
            aaiB.append(aai)
            
    return latA, lonA, aaiA, latB, lonB, aaiB 
            
            








#    # GOME
def gridGOME2(year, month, day, ROI, res):
    print '**** Reading GOME2 %02i-%02i-%04i' %(day,month,year) 
    
    pathA = '/nobackup/users/sunj/GOME-2/MetopA/%4i/%02i/%02i/' %(year, month, day)
    filelistA = glob.glob( pathA + '*.hdf5')
    pathB = '/nobackup/users/sunj/GOME-2/MetopB/%4i/%02i/%02i/' %(year, month, day)
    filelistB = glob.glob( pathB + '*.hdf5')
  
    
 
    aaiA = []
    latA = []
    lonA = [] 
    maskA = []   
    raaA = [] 
    scaA = []
    
    for io in filelistA[:]:  

        dataA = netCDF4.Dataset(io,'r')    
        data = dataA.groups['DATA']
        geo = dataA.groups['GEOLOCATION']

        lat = geo.variables['LatitudeCenter'][:]
        lon = geo.variables['LongitudeCenter'][:]
        aai = data.variables['AAI'][:]
        sza = geo.variables['SolarZenithAngle'][:]
        sunglint = data.variables['SunGlintFlag'][:]  
        saa = geo.variables['SolarAzimuthAngle'][:]
        vaa = geo.variables['LineOfSightAzimuthAngle'][:]
        sca = geo.variables['ScatteringAngle'][:]
        
        raa = saa + 180.- vaa 
        
        mask = np.logical_or(sza[0:-1,0:-1] < 0. , sza[0:-1,0:-1] > 75.)  
        mask = np.logical_or(mask, lat[0:-1,0:-1]<ROI[0])
        mask = np.logical_or(mask, lat[0:-1,0:-1]>ROI[1])
        mask = np.logical_or(mask, lon[0:-1,0:-1]<ROI[2])
        mask = np.logical_or(mask, lon[0:-1,0:-1]>ROI[3])
        mask = np.logical_or(mask, sunglint[0:-1,0:-1] < 0.)                # sunglint pixel
        mask = np.logical_or(mask, aai[0:-1,0:-1] < 0.)                     # scattering aerosol pixel  
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel 
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
#        mask = np.logical_or(mask, aai[0:-1,0:-1]<5.)
        
        latA = latA + list(np.asarray(lat[0:-1,0:-1]).reshape(-1))
        lonA = lonA + list(np.asarray(lon[0:-1,0:-1]).reshape(-1))
        aaiA = aaiA + list(np.asarray(aai[0:-1,0:-1]).reshape(-1))
        raaA = raaA + list(np.asarray(raa[0:-1,0:-1]).reshape(-1))
        scaA = scaA + list(np.asarray(sca[0:-1,0:-1]).reshape(-1))
        maskA = maskA + list(np.asarray(mask).reshape(-1))
        
        dataA.close()

    latarrA = np.array(latA)
    lonarrA = np.array(lonA)    
    aaiarrA = np.array(aaiA) 
    raaarrA = np.array(raaA)
    scaarrA = np.array(scaA)
    maskarrA = np.array(maskA)       
    
    aaiarrA = np.ma.masked_array(aaiarrA,maskarrA)

    latvldA = latarrA[aaiarrA.mask == 0] 
    lonvldA = lonarrA[aaiarrA.mask == 0] 
    aaivldA = aaiarrA[aaiarrA.mask == 0]     
    raavldA = raaarrA[aaiarrA.mask == 0]     
    scavldA = scaarrA[aaiarrA.mask == 0]     


 
    aaiB = []
    latB = []
    lonB = [] 
    maskB = []
    raaB = [] 
    scaB = [] 
    
    for io in filelistB[:]:  

        dataB = netCDF4.Dataset(io,'r')    
        data = dataB.groups['DATA']
        geo = dataB.groups['GEOLOCATION']

        lat = geo.variables['LatitudeCenter'][:]
        lon = geo.variables['LongitudeCenter'][:]
        aai = data.variables['AAI'][:]
        sza = geo.variables['SolarZenithAngle'][:]
        sunglint = data.variables['SunGlintFlag'][:]  
        saa = geo.variables['SolarAzimuthAngle'][:]
        vaa = geo.variables['LineOfSightAzimuthAngle'][:]
        sca = geo.variables['ScatteringAngle'][:]
        
        raa = saa + 180.- vaa 
        
        mask = np.logical_or(sza[0:-1,0:-1] < 0. , sza[0:-1,0:-1] > 75.)  
        mask = np.logical_or(mask, lat[0:-1,0:-1]<ROI[0])
        mask = np.logical_or(mask, lat[0:-1,0:-1]>ROI[1])
        mask = np.logical_or(mask, lon[0:-1,0:-1]<ROI[2])
        mask = np.logical_or(mask, lon[0:-1,0:-1]>ROI[3])
        mask = np.logical_or(mask, sunglint[0:-1,0:-1] < 0.)                # sunglint pixel
        mask = np.logical_or(mask, aai[0:-1,0:-1] < 0.)                     # scattering aerosol pixel  
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel 
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
#        mask = np.logical_or(mask, aai[0:-1,0:-1]<5.)
        
        latB = latB + list(np.asarray(lat[0:-1,0:-1]).reshape(-1))
        lonB = lonB + list(np.asarray(lon[0:-1,0:-1]).reshape(-1))
        aaiB = aaiB + list(np.asarray(aai[0:-1,0:-1]).reshape(-1))
        raaB = raaB + list(np.asarray(raa[0:-1,0:-1]).reshape(-1))
        scaB = scaB + list(np.asarray(sca[0:-1,0:-1]).reshape(-1))
        maskB = maskB + list(np.asarray(mask).reshape(-1))

    latarrB = np.array(latB)
    lonarrB = np.array(lonB)    
    aaiarrB = np.array(aaiB) 
    raaarrB = np.array(raaB) 
    scaarrB = np.array(scaB) 
    maskarrB = np.array(maskB)       
    
    aaiarrB = np.ma.masked_array(aaiarrB,maskarrB)

    latvldB = latarrB[aaiarrB.mask == 0] 
    lonvldB = lonarrB[aaiarrB.mask == 0] 
    aaivldB = aaiarrB[aaiarrB.mask == 0]        
    raavldB = raaarrB[aaiarrB.mask == 0]        
    scavldB = scaarrB[aaiarrB.mask == 0]
    
    latvld = np.concatenate((latvldA, latvldB))
    lonvld = np.concatenate((lonvldA, lonvldB))
    aaivld = np.concatenate((aaivldA, aaivldB))
    raavld = np.concatenate((raavldA, raavldB))
    scavld = np.concatenate((scavldA, scavldB))
    
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
    
    x,y = np.meshgrid(lonnew,latnew)
    aainew = griddata((lonvld, latvld), aaivld, (x, y), method = 'linear')
    aainew[np.isnan(aainew)] = -32767.
    
#    raavld[raavld < 0] = 360 + raavld[raavld<0]
#    plt.figure()
#    plt.subplot(1,2,1)
#    plt.scatter(raavldA, aaivldA, c ='k')
#    plt.scatter(raavldB, aaivldB, c ='r')
#    plt.subplot(1,2,2)
#    plt.scatter(scavldA, aaivldA, c = 'k')    
#    plt.scatter(scavldB, aaivldB, c = 'r')    
    
    return latnew, lonnew, aainew 





def GOME24cfg(year, month, day, ROI, jetlag, plumemsk, crival, res):     
    
    print '**** Reading GOME2 %02i-%02i-%04i' %(day,month,year) 
    
    pathA = '/nobackup/users/sunj/GOME-2/MetopA/%4i/%02i/%02i/' %(year, month, day)
    filelistA = glob.glob( pathA + '*.hdf5')
    pathB = '/nobackup/users/sunj/GOME-2/MetopB/%4i/%02i/%02i/' %(year, month, day)
    filelistB = glob.glob( pathB + '*.hdf5')
  
    
 
    aaiA = []
    latA = []
    lonA = [] 
    maskA = []   
    raaA = [] 
    scaA = []
    
    for io in filelistA[:]:  

        dataA = netCDF4.Dataset(io,'r')    
        data = dataA.groups['DATA']
        geo = dataA.groups['GEOLOCATION']

        lat = geo.variables['LatitudeCenter'][:]
        lon = geo.variables['LongitudeCenter'][:]
        aai = data.variables['AAI'][:]
        sunglint = data.variables['SunGlintFlag'][:]  
        
        saa = geo.variables['SolarAzimuthAngle'][:]
        sza = geo.variables['SolarZenithAngle'][:]
        vaa = geo.variables['LineOfSightAzimuthAngle'][:]
        vza = geo.variables['LineOfSightZenithAngle'][:]
        sca = geo.variables['ScatteringAngle'][:]
        raa = geo.variables['RelAzimuthAngle'][:]
        raa1 = saa + 180.- vaa 
        
        print (raa - raa1).mean()          
        
        mask = np.logical_or(sza[0:-1,0:-1] < 0. , sza[0:-1,0:-1] > 75.)  
        mask = np.logical_or(mask, lat[0:-1,0:-1]<ROI[0])
        mask = np.logical_or(mask, lat[0:-1,0:-1]>ROI[1])
        mask = np.logical_or(mask, lon[0:-1,0:-1]<ROI[2])
        mask = np.logical_or(mask, lon[0:-1,0:-1]>ROI[3])
        mask = np.logical_or(mask, sunglint[0:-1,0:-1] < 0.)                # sunglint pixel
        mask = np.logical_or(mask, aai[0:-1,0:-1] < 0.)                     # scattering aerosol pixel  
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel 
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
#        mask = np.logical_or(mask, aai[0:-1,0:-1]<5.)
        
        latA = latA + list(np.asarray(lat[0:-1,0:-1]).reshape(-1))
        lonA = lonA + list(np.asarray(lon[0:-1,0:-1]).reshape(-1))
        aaiA = aaiA + list(np.asarray(aai[0:-1,0:-1]).reshape(-1))
        raaA = raaA + list(np.asarray(raa[0:-1,0:-1]).reshape(-1))
        scaA = scaA + list(np.asarray(sca[0:-1,0:-1]).reshape(-1))
        maskA = maskA + list(np.asarray(mask).reshape(-1))
        
        dataA.close()

    latarrA = np.array(latA)
    lonarrA = np.array(lonA)    
    aaiarrA = np.array(aaiA) 
    raaarrA = np.array(raaA)
    scaarrA = np.array(scaA)
    maskarrA = np.array(maskA)       
    
    aaiarrA = np.ma.masked_array(aaiarrA,maskarrA)

    latvldA = latarrA[aaiarrA.mask == 0] 
    lonvldA = lonarrA[aaiarrA.mask == 0] 
    aaivldA = aaiarrA[aaiarrA.mask == 0]     
    raavldA = raaarrA[aaiarrA.mask == 0]     
    scavldA = scaarrA[aaiarrA.mask == 0]     


 
    aaiB = []
    latB = []
    lonB = [] 
    maskB = []
    raaB = [] 
    scaB = [] 
    
    for io in filelistB[:]:  

        dataB = netCDF4.Dataset(io,'r')    
        data = dataB.groups['DATA']
        geo = dataB.groups['GEOLOCATION']

        lat = geo.variables['LatitudeCenter'][:]
        lon = geo.variables['LongitudeCenter'][:]
        aai = data.variables['AAI'][:]
        sza = geo.variables['SolarZenithAngle'][:]
        sunglint = data.variables['SunGlintFlag'][:]  
        saa = geo.variables['SolarAzimuthAngle'][:]
        vaa = geo.variables['LineOfSightAzimuthAngle'][:]
        sca = geo.variables['ScatteringAngle'][:]
        
        raa = saa + 180.- vaa 
        
        mask = np.logical_or(sza[0:-1,0:-1] < 0. , sza[0:-1,0:-1] > 75.)  
        mask = np.logical_or(mask, lat[0:-1,0:-1]<ROI[0])
        mask = np.logical_or(mask, lat[0:-1,0:-1]>ROI[1])
        mask = np.logical_or(mask, lon[0:-1,0:-1]<ROI[2])
        mask = np.logical_or(mask, lon[0:-1,0:-1]>ROI[3])
        mask = np.logical_or(mask, sunglint[0:-1,0:-1] < 0.)                # sunglint pixel
        mask = np.logical_or(mask, aai[0:-1,0:-1] < 0.)                     # scattering aerosol pixel  
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel 
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
#        mask = np.logical_or(mask, aai[0:-1,0:-1]<5.)
        
        latB = latB + list(np.asarray(lat[0:-1,0:-1]).reshape(-1))
        lonB = lonB + list(np.asarray(lon[0:-1,0:-1]).reshape(-1))
        aaiB = aaiB + list(np.asarray(aai[0:-1,0:-1]).reshape(-1))
        raaB = raaB + list(np.asarray(raa[0:-1,0:-1]).reshape(-1))
        scaB = scaB + list(np.asarray(sca[0:-1,0:-1]).reshape(-1))
        maskB = maskB + list(np.asarray(mask).reshape(-1))

    latarrB = np.array(latB)
    lonarrB = np.array(lonB)    
    aaiarrB = np.array(aaiB) 
    raaarrB = np.array(raaB) 
    scaarrB = np.array(scaB) 
    maskarrB = np.array(maskB)       
    
    aaiarrB = np.ma.masked_array(aaiarrB,maskarrB)

    latvldB = latarrB[aaiarrB.mask == 0] 
    lonvldB = lonarrB[aaiarrB.mask == 0] 
    aaivldB = aaiarrB[aaiarrB.mask == 0]        
    raavldB = raaarrB[aaiarrB.mask == 0]        
    scavldB = scaarrB[aaiarrB.mask == 0]
    
    latvld = np.concatenate((latvldA, latvldB))
    lonvld = np.concatenate((lonvldA, lonvldB))
    aaivld = np.concatenate((aaivldA, aaivldB))
    raavld = np.concatenate((raavldA, raavldB))
    scavld = np.concatenate((scavldA, scavldB))
    
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
    
    x,y = np.meshgrid(lonnew,latnew)
    aainew = griddata((lonvld, latvld), aaivld, (x, y), method = 'linear')
    aainew[np.isnan(aainew)] = -32767.
    
#    raavld[raavld < 0] = 360 + raavld[raavld<0]
#    plt.figure()
#    plt.subplot(1,2,1)
#    plt.scatter(raavldA, aaivldA, c ='k')
#    plt.scatter(raavldB, aaivldB, c ='r')
#    plt.subplot(1,2,2)
#    plt.scatter(scavldA, aaivldA, c = 'k')    
#    plt.scatter(scavldB, aaivldB, c = 'r')    
    
    return latnew, lonnew, aainew 

#year = 2017
#month = 10
#day = 17
#
## region of interest 
#ROI = [45,65,30,60]
#
#res = 0.5
#
#GOME24cfg(year, month, day, ROI, jetlag, plumemsk, 2, res)