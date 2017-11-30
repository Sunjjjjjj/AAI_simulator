# -*- coding: utf-8 -*-
"""
Read OMAERO data 
- rawOMAERO: return orbit data, organized by list
- gridOMAERO: return gridded data over ROI, organized by np.2darray
- daytsOMAERO: return time series 

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
import pandas as pd
from scipy import spatial
import matplotlib.mlab as mlab
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata


def rawOMAERO(year, month, day, jetlag, parameter): 
    # read OMAERO orbit data, global coverage 


    aerdir = '/nobackup/users/sunj/AERONET/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340
    aerfile = glob.glob( aerdir + '*.lev15')[0]
    aerdata = pd.read_csv(aerfile, sep=",", header = 4)
    
    ts1 = aerdata['Julian_Day']
    alpha = aerdata['440-675Angstrom']
    
    
  
    
    
    print '**** Reading OMAEROv003 %02i-%02i-%04i' %(day,month,year)  
    paragol = []
    latgol = []
    longol = []

    
    path = '/nobackup/users/sunj/OMI/OMAERO/%4i/%02i/%02i/' %(year, month, day)
    filelist = glob.glob( path + '*.he5')
    
    
    for io in filelist[:]: 

        orbdata = netCDF4.Dataset(io,'r')
   
        data = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Data Fields']
        geo = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Geolocation Fields']
        add = orbdata.groups['HDFEOS'].groups['ADDITIONAL'].groups['FILE_ATTRIBUTES']

        lat = geo.variables['Latitude'][:]
        lon = geo.variables['Longitude'][:]    
        meatime = geo.variables['Time'][:] 
        sza = geo.variables['SolarZenithAngle'][:]
        vza = geo.variables['ViewingZenithAngle'][:]
        saa = geo.variables['SolarAzimuthAngle'][:]
        vaa = geo.variables['ViewingAzimuthAngle'][:]
        gpqf = geo.variables['GroundPixelQualityFlags'][:]          
        row = data.variables['XTrackQualityFlags'][:]    
        dt = add.TAI93At0zOfGranule
        wvl = add.OPF_Wavelengths


#        secday = meatime - dt + jetlag * 3600.                             # - time correction + jetlat second of the day         
#        tstemp = day + secday / 24. / 3600. 
#        alpha1 = np.interp(tstemp, tsaer, alpha)
#        
#        alpha1 = np.array([alpha1] * lat.shape[1]).transpose()  


#        alpha2 = np.interp(ts2,ts1,alpha)

        if parameter == 'AOD': 
            para= data.variables['AerosolOpticalThicknessMW'][:,:,0]*1e-3           # aod [lat, lon, wvl] 
        elif parameter == 'AAI':
            para = data.variables['AerosolIndexUV'][:]*1e-2
        
        
        # print 'Calculating flag'
        szarad = sza/180.*np.pi 
        vzarad = vza/180.*np.pi
        
        raz = np.fabs(vaa - saa)
        raz = raz.data            
        raa = 180 - np.where(raz>180,360-raz,raz)
        raarad = raa/180.*np.pi
        
        scatangle = np.arccos(np.cos(szarad)*np.cos(vzarad) - np.sin(szarad)*np.sin(vzarad)*np.cos(np.pi-raarad))       
        waterflag = np.where(np.bitwise_and(gpqf,15) == 1, False, True)        
        sunglint = np.logical_and(waterflag,np.where(np.rad2deg(scatangle)<30., True, False))
        eclipse = np.where(np.bitwise_and(gpqf,32),True,False) 


        mask = np.logical_or(sunglint[0:-1,0:-1],eclipse[0:-1,0:-1])        # sunglint pixel, sun eclipse pixel                      
        mask = np.logical_or(mask, row[0:-1,0:-1]>0)                        # cross track row anomaly 
        mask = np.logical_or(mask, sza[0:-1,0:-1]>60.)
#        mask = np.logical_or(mask, vza[0:-1,0:-1]<35.)
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
        mask = np.logical_or(mask, lat[0:-1,0:-1]>lat[1:,0:-1])             # decreasing oribit pixel 

        if mask.all() == False: 
            lon = np.ma.masked_array(lon[0:-1,0:-1],mask)
            lat = np.ma.masked_array(lat[0:-1,0:-1],mask)
            para = np.ma.masked_array(para[0:-1,0:-1],mask)  
          

    
            longol.append(lon)
            latgol.append(lat)
            paragol.append(para)
            
        
    return latgol, longol, paragol


#lat, lon, para = rawOMAERO(2017, 1, 26, -5, 'AAI')
#
#plt.figure()
#
#for i in np.arange(0,len(lat)):
#    para[i][para[i].mask == 1]
#    plt.pcolor(np.array(lon[i]),np.array(lat[i]),np.array(para[i]),cmap = 'jet', vmin= 0, vmax = 10)
  
  

def gridOMAERO(year, month, day, ROI, jetlag, parameter, res): 
    # read OMAERO data, grid over ROI 




    aerdir = '/nobackup/users/sunj/AERONET/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340
    aerfile = glob.glob( aerdir + '*.lev15')[0]
    aerdata = pd.read_csv(aerfile, sep=",", header = 4)
    
    tsaer = aerdata['Julian_Day']
    alpha = aerdata['440-675Angstrom']
#    alpha2 = aerdata['440-870Angstrom']
#    plt.figure()
#    plt.plot(tsaer, alpha,'k-')
#    plt.plot(tsaer, alpha2, 'b-')
#    plt.title('Angstrom exponent')

    print '**** Reading OMI %02i-%02i-%04i' %(day,month,year)     
    paragol = []
    latgol = []
    longol = [] 
    maskgol = []   
    
    path = '/nobackup/users/sunj/OMI/OMAERO/%4i/%02i/%02i/' %(year, month, day)
    filelist = glob.glob( path + '*.he5')
    
    for io in filelist[:]: 

        orbdata = netCDF4.Dataset(io,'r')
   
        data = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Data Fields']
        geo = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Geolocation Fields']
        add = orbdata.groups['HDFEOS'].groups['ADDITIONAL'].groups['FILE_ATTRIBUTES']

        lat = geo.variables['Latitude'][:]
        lon = geo.variables['Longitude'][:]       
        meatime = geo.variables['Time'][:] 
        sza = geo.variables['SolarZenithAngle'][:]
        vza = geo.variables['ViewingZenithAngle'][:]
        saa = geo.variables['SolarAzimuthAngle'][:]
        vaa = geo.variables['ViewingAzimuthAngle'][:]
        gpqf = geo.variables['GroundPixelQualityFlags'][:] 
        dt = add.TAI93At0zOfGranule
        
        secday = meatime - dt + jetlag * 3600.                             # - time correction + jetlat second of the day         
        tstemp = day + secday / 24. / 3600. 
        alpha1 = np.interp(tstemp, tsaer, alpha)
        
#        plt.plot(tstemp,alpha1,'r.')
#        print alpha1.mean()
        
        alpha1 = np.array([alpha1] * lat.shape[1]).transpose()  
       

#        aod = data.variables['AerosolOpticalThicknessMW'][:]*1e-3       # aod [lat, lon, wvl] 
#        aai = data.variables['AerosolIndexUV'][:]*1e-2
        row = data.variables['XTrackQualityFlags'][:]    

        if parameter == 'AOD': 
            aod = data.variables['AerosolOpticalThicknessMW'][:]*1e-3           # aod [lat, lon, wvl] 
#            aod = data.variables['AerosolOpticalThicknessPassedThresholdMean'][:]*1e-3           # aod [lat, lon, wvl] 
#            aodstd = data.variables['AerosolOpticalThicknessPassedThresholdStd'][:]*1e-3           # aod [lat, lon, wvl] 
#            para = aod[:,:,-5] * (550./442)** (-alpha1)
#            para = aod[:,:,0] * (550./342.5)** (-alpha1)
#            para = aod[:,:,:].mean(2)
            para = aod[:,:,-1]          # 483.5 nm 
#            print aod[:,:,0].max()
#            print para.max()
#            print para1.max()
        elif parameter == 'AAI':
            para = data.variables['AerosolIndexUV'][:]*1e-2        
                
        # print 'Calculating flag'
        szarad = sza/180.*np.pi 
        vzarad = vza/180.*np.pi
        
        raz = np.fabs(vaa - saa)
        raz = raz.data            
        raa = 180 - np.where(raz>180,360-raz,raz)
        raarad = raa/180.*np.pi
        
        scatangle = np.arccos(np.cos(szarad)*np.cos(vzarad) - np.sin(szarad)*np.sin(vzarad)*np.cos(np.pi-raarad))       
        waterflag = np.where(np.bitwise_and(gpqf,15) == 1, False, True)        
        sunglint = np.logical_and(waterflag,np.where(np.rad2deg(scatangle)<30., True, False))
        eclipse = np.where(np.bitwise_and(gpqf,32),True,False) 


        
        mask = np.logical_or(lat[0:-1,0:-1]>ROI[1], lat[0:-1,0:-1]<ROI[0])
        mask = np.logical_or(mask, lon[0:-1,0:-1]<ROI[2])
        mask = np.logical_or(mask, lon[0:-1,0:-1]>ROI[3])
        mask = np.logical_or(mask, eclipse[0:-1,0:-1])        # sunglint pixel, sun eclipse pixel 
        mask = np.logical_or(mask,sunglint[0:-1,0:-1] )
        mask = np.logical_or(mask, row[0:-1,0:-1]>0)                        # cross track row anomaly 
        mask = np.logical_or(mask, sza[0:-1,0:-1]>75.)
#        mask = np.logical_or(mask, vza[0:-1,0:-1]<35.)
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
        mask = np.logical_or(mask, lat[0:-1,0:-1]>lat[1:,0:-1])             # decreasing oribit pixel 
#        mask = np.logical_or(mask, aai[0:-1,0:-1]<5.)

        latgol = latgol + list(np.asarray(lat[0:-1,0:-1]).reshape(-1))
        longol = longol + list(np.asarray(lon[0:-1,0:-1]).reshape(-1))
        paragol = paragol + list(np.asarray(para[0:-1,0:-1]).reshape(-1))
        maskgol = maskgol + list(np.asarray(mask).reshape(-1))

    latarr = np.array(latgol)
    lonarr = np.array(longol)    
    paraarr = np.array(paragol) 
    maskarr = np.array(maskgol)
    
#    print 'Applying mask'
    paraarr = np.ma.masked_array(paraarr,maskarr)
    
    latvld = np.array(latarr[paraarr.mask == 0]) 
    lonvld = np.array(lonarr[paraarr.mask == 0]) 
    paravld = paraarr[paraarr.mask == 0]
    
    
#    print 'Grid over ROI'
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
     
    x,y = np.meshgrid(lonnew,latnew)
    paranew = griddata((lonvld, latvld), paravld, (x, y), method = 'linear')
#    aainew = griddata(lonvld, latvld, aaivld, lonnew, latnew, interp = 'linear')
#    paranew[np.isnan(paranew)] = -32767.
    
    return latnew, lonnew, paranew




def daytsOMAERO(year, month, days, aerlat, aerlon, jetlag): 

#    print 'Calculate time zone'
#    if abs(aerlon)%15 < 7.5:MODIS_OMI_statistics.py
#        jetlag = abs(aerlon)//15
#    else:
#        jetlag = abs(aerlon)//15+1
#        
#    if aerlon<0:
#        jetlag = - jetlag 
#    print jetlag 

    
    aodts1 = []
    aodts2 = []
    ts = []
    
    for iday in days[:]:
        print '**** Reading OMI %02i-%02i-%04i' %(iday, month, year) 
        aodgol1 = []
        aodgol2 = []
        latgol = []
        longol = [] 
        tsgol = [] 
        
        
        path = '/nobackup/users/sunj/OMI/OMAERO/%4i/%02i/%02i/' %(year, month, iday)
        filelist = glob.glob( path + '*.he5')
    
    
        for io in filelist[:]: 
            orbdata = netCDF4.Dataset(io,'r')
            
            data = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Data Fields']
            geo = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Geolocation Fields']
            add = orbdata.groups['HDFEOS'].groups['ADDITIONAL'].groups['FILE_ATTRIBUTES']
            
            lat = geo.variables['Latitude'][:]
            lon = geo.variables['Longitude'][:] 
            meatime = geo.variables['Time'][:]                                 # measurement time
            aod = data.variables['AerosolOpticalThicknessMW'][:]*1e-3          # aod [lat, lon, wvl] 
            dt = add.TAI93At0zOfGranule            
            wvl = add.OPF_Wavelengths
            
    
    #        print '**** Transferring TAI time to Julian Day'
            secday = meatime - dt + jetlag * 3600.                             # - time correction + jetlat second of the day         
            tstemp = iday + secday / 24. / 3600. 
            tsomi = np.array([tstemp] * lat.shape[1]).transpose()  
            
    #        print '**** Applying Angstrom exponent to calculate AOD@550nm'
    #        alpha = np.log(aod[:,:,0]/aod[:,:,-1]) / np.log(483.5/342.5)
    #        aod550 = aod[:,:,0] * (550./342.5)**(-alpha) 
            aod342 = aod[:,:,0]
            aod442 = aod[:,:,-5]
            
    #        print '**** Collecting all orbits for one specific day'
            latgol = latgol + list(np.asarray(lat).reshape(-1))
            longol = longol + list(np.asarray(lon).reshape(-1))
            aodgol1 = aodgol1 + list(np.asarray(aod342).reshape(-1))
            aodgol2 = aodgol2 + list(np.asarray(aod442).reshape(-1))
            tsgol = tsgol +list(np.asarray(tsomi).reshape(-1))
            
            orbdata.close()
    #    print '**** Collecting all orbits for one specific day'
        latarr = np.array(latgol)
        lonarr = np.array(longol)    
        aodarr1 = np.array(aodgol1) 
        aodarr2 = np.array(aodgol2) 
        tsarr = np.array(tsgol)
        
    #    print '**** Searching the nearest pixel near Santiago_Beauchef AERONET station'
        tree = spatial.KDTree(zip(latarr, lonarr))
        locs = tree.query([aerlat,aerlon], k = 1)
        
        aodts1.append(aodarr1[locs[1]])
        aodts2.append(aodarr2[locs[1]])
        ts.append(tsarr[locs[1]])
    #    print latarr[locs[1]], lonarr[locs[1]]
        
    
    #print '**** Masking data'
    aodts1 = np.array(aodts1)
    aodts2 = np.array(aodts2)
    mask = np.zeros(aodts1.shape)
    mask[aodts1<0]=1
    aod342 = np.ma.masked_array(aodts1,mask)
    aod442 = np.ma.masked_array(aodts2,mask)    
    
    return ts, aod342, aod442






def OMAERO4cfg(year, month, day, ROI, jetlag, plumemsk, crival, res):     
    
    aerdir = '/nobackup/users/sunj/AERONET/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340
    aerfile = glob.glob( aerdir + '*.lev15')[0]
    aerdata = pd.read_csv(aerfile, sep=",", header = 4)
    
    tsaer = aerdata['Julian_Day']
    alpha = aerdata['440-675Angstrom']

    print '**** Reading OMI for Config.in %02i-%02i-%04i' %(day,month,year)     
    aaigol = []    
    aodgol = []
    aodstdgol = []
    latgol = []
    longol = [] 
    szagol = []
    vzagol = [] 
    saagol = []
    vaagol = []
    Hsgol = []
    Psgol = []
    As1gol = [] 
    As2gol = [] 
    ssagol = [] 
    
    maskgol = []   
    
    path = '/nobackup/users/sunj/OMI/OMAERO/%4i/%02i/%02i/' %(year, month, day)
    filelist = glob.glob( path + '*.he5')
    
    for io in filelist[:]: 

        orbdata = netCDF4.Dataset(io,'r')
   
        data = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Data Fields']
        geo = orbdata.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Geolocation Fields']
        add = orbdata.groups['HDFEOS'].groups['ADDITIONAL'].groups['FILE_ATTRIBUTES']

        lat = geo.variables['Latitude'][:]
        lon = geo.variables['Longitude'][:]       
        meatime = geo.variables['Time'][:] 
        sza = geo.variables['SolarZenithAngle'][:]
        vza = geo.variables['ViewingZenithAngle'][:]
        saa = geo.variables['SolarAzimuthAngle'][:]
        vaa = geo.variables['ViewingAzimuthAngle'][:]
        Hs = geo.variables['TerrainHeight'][:]
        Ps = data.variables['TerrainPressure'][:]
        As1 = data.variables['TerrainReflectivity'][:,:,1]*1e-3
        As2 = data.variables['TerrainReflectivity'][:,:,2]*1e-3
#        As = np.ccolumn_stack([As1, As2])
        cf = data.variables['EffectiveCloudFraction'][:]*1e-2
        ssa = data.variables['SingleScatteringAlbedoPassedThresholdMean'][:,:,0]*1e-3  
        gpqf = geo.variables['GroundPixelQualityFlags'][:] 
        dt = add.TAI93At0zOfGranule
        
        
        secday = meatime - dt + jetlag * 3600.                             # - time correction + jetlat second of the day         
        tstemp = day + secday / 24. / 3600. 
        alpha1 = np.interp(tstemp, tsaer, alpha)
        alpha1 = np.array([alpha1] * lat.shape[1]).transpose()  
#        print 'Angstrom exponent OMI:', alpha1

        row = data.variables['XTrackQualityFlags'][:]    

        aai = data.variables['AerosolIndexUV'][:]*1e-2
        aod442 = data.variables['AerosolOpticalThicknessMW'][:,:,-5]*1e-3           # aod [lat, lon, wvl] 
        aodstd = data.variables['AerosolOpticalThicknessMWPrecision'][:]*1e-3           # aod [lat, lon, wvl] 
        aod = aod442 * (550./442)** (-alpha1)
        
                
        # print 'Calculating flag'
        szarad = sza/180.*np.pi 
        vzarad = vza/180.*np.pi
        
        raz = np.fabs(vaa - saa)
        raz = raz.data            
        raa = 180 - np.where(raz>180,360-raz,raz)
        raarad = raa/180.*np.pi
        
        scatangle = np.arccos(np.cos(szarad)*np.cos(vzarad) - np.sin(szarad)*np.sin(vzarad)*np.cos(np.pi-raarad))       
        waterflag = np.where(np.bitwise_and(gpqf,15) == 1, False, True)        
        sunglint = np.logical_and(waterflag,np.where(np.rad2deg(scatangle)<30., True, False))
        eclipse = np.where(np.bitwise_and(gpqf,32),True,False) 


        mask = np.logical_or(lat[0:-1,0:-1]>ROI[1], lat[0:-1,0:-1]<ROI[0])
        mask = np.logical_or(mask, lon[0:-1,0:-1]<ROI[2])
        mask = np.logical_or(mask, lon[0:-1,0:-1]>ROI[3])
        mask = np.logical_or(mask, sunglint[0:-1,0:-1]) # sunglint pixel, sun eclipse pixel
        mask = np.logical_or(mask, eclipse[0:-1,0:-1])
        mask = np.logical_or(mask, row[0:-1,0:-1]>0)                        # cross track row anomaly 
        mask = np.logical_or(mask, sza[0:-1,0:-1]>75.)
#        mask = np.logical_or(mask, vza[0:-1,0:-1]<35.)
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[1:,0:-1])<-100)      # dateline pixel
        mask = np.logical_or(mask, (lon[0:-1,0:-1]*lon[0:-1,1:])<-100)
        mask = np.logical_or(mask, lat[0:-1,0:-1]>lat[1:,0:-1])             # decreasing oribit pixel 
#        mask = np.logical_or(mask, cf[0:-1,0:-1]>0.35)
#        mask = np.logical_or(mask, As1[0:-1,0:-1]>0.1)
#        mask = np.logical_or(mask, As2[0:-1,0:-1]>0.1)
        

        latgol = latgol + list(np.asarray(lat[0:-1,0:-1]).reshape(-1))
        longol = longol + list(np.asarray(lon[0:-1,0:-1]).reshape(-1))
        aaigol = aaigol + list(np.asarray(aai[0:-1,0:-1]).reshape(-1))
        aodgol = aodgol + list(np.asarray(aod[0:-1,0:-1]).reshape(-1))
        aodstdgol = aodstdgol + list(np.asarray(aodstd[0:-1,0:-1]).reshape(-1))
        szagol = szagol + list(np.asarray(sza[0:-1,0:-1]).reshape(-1))
        saagol = saagol + list(np.asarray(saa[0:-1,0:-1]).reshape(-1))
        vzagol = vzagol + list(np.asarray(vza[0:-1,0:-1]).reshape(-1))
        vaagol = vaagol + list(np.asarray(vaa[0:-1,0:-1]).reshape(-1))
        Hsgol = Hsgol + list(np.asarray(Hs[0:-1,0:-1]).reshape(-1))
        Psgol = Psgol + list(np.asarray(Ps[0:-1,0:-1]).reshape(-1))
        As1gol = As1gol + list(np.asarray(As1[0:-1,0:-1]).reshape(-1))
        As2gol = As2gol + list(np.asarray(As2[0:-1,0:-1]).reshape(-1))
        ssagol = ssagol + list(np.asarray(ssa[0:-1,0:-1]).reshape(-1))
        
        maskgol = maskgol + list(np.asarray(mask).reshape(-1))
        
    
    latarr = np.array(latgol)
    lonarr = np.array(longol)    
    aaiarr = np.array(aaigol) 
    aodarr = np.array(aodgol) 
    aodstdarr = np.array(aodstdgol)
    szaarr = np.array(szagol)
    saaarr = np.array(saagol)    
    vzaarr = np.array(vzagol) 
    vaaarr = np.array(vaagol)
    Psarr = np.array(Psgol)    
    Hsarr = np.array(Hsgol) 
    As1arr = np.array(As1gol)
    As2arr = np.array(As2gol)
    ssaarr = np.array(ssagol)

    


    maskarr = np.array(maskgol)
#    print 'Applying mask'
    aodarr = np.ma.masked_array(aodarr,maskarr)
    
    latvld = np.array(latarr[aodarr.mask == 0]) 
    lonvld = np.array(lonarr[aodarr.mask == 0]) 
    aaivld = np.array(aaiarr[aodarr.mask == 0])
    aodvld = np.array(aodarr[aodarr.mask == 0])
    aodstdvld = np.array(aodstdarr[aodarr.mask == 0])
    szavld = np.array(szaarr[aodarr.mask == 0]) 
    saavld = np.array(saaarr[aodarr.mask == 0]) 
    vzavld = np.array(vzaarr[aodarr.mask == 0])
    vaavld = np.array(vaaarr[aodarr.mask == 0]) 
    Psvld = np.array(Psarr[aodarr.mask == 0]) 
    Hsvld = np.array(Hsarr[aodarr.mask == 0]) 
    As1vld = np.array(As1arr[aodarr.mask == 0]) 
    As2vld = np.array(As2arr[aodarr.mask == 0]) 
    ssavld = np.array(ssaarr[aodarr.mask == 0]) 
    
   

   
#    print 'Grid over ROI'
    latnew = np.arange(ROI[0], ROI[1], res)
    lonnew = np.arange(ROI[2], ROI[3], res)
     
     
    x,y = np.meshgrid(lonnew,latnew)
    aainew = griddata((lonvld, latvld), aaivld, (x, y), method = 'linear')
    aodnew = griddata((lonvld, latvld), aodvld, (x, y), method = 'linear')
    aodstdnew = griddata((lonvld, latvld), aodstdvld, (x, y), method = 'linear')
    szanew = griddata((lonvld, latvld), szavld, (x, y), method = 'linear')
    saanew = griddata((lonvld, latvld), saavld, (x, y), method = 'linear')
    vzanew = griddata((lonvld, latvld), vzavld, (x, y), method = 'linear')
    vaanew = griddata((lonvld, latvld), vaavld, (x, y), method = 'linear')
    Psnew = griddata((lonvld, latvld), Psvld, (x, y), method = 'linear')
    Hsnew = griddata((lonvld, latvld), Hsvld, (x, y), method = 'linear')
    As1new = griddata((lonvld, latvld), As1vld, (x, y), method = 'linear')
    As2new = griddata((lonvld, latvld), As2vld, (x, y), method = 'linear')
    ssanew = griddata((lonvld, latvld), ssavld, (x, y), method = 'linear')
    
    
    # set NaN to negative vales
    aainew[np.isnan(aainew)] = -32767.
    aodnew[np.isnan(aodnew)] = -32767.
    plumemsk = np.logical_or(plumemsk, aodnew < 0.)
    plumemsk = np.logical_or(plumemsk, aainew < crival)
#    Psnew[np.isnan(szanew)] = -32767.
#    Psnew[np.isnan(vzanew)] = -32767.
#    Psnew[np.isnan(saanew)] = -32767.
#    Psnew[np.isnan(vaanew)] = -32767.
#    Psnew[np.isnan(Psnew)] = -32767.
#    Psnew[np.isnan(As1new)] = -32767.
#    Psnew[np.isnan(As2new)] = -32767.
#    plumemsk = np.logical_or(plumemsk, Psnew < 0.)


    latm = np.ma.masked_array(y,plumemsk)  
    lonm = np.ma.masked_array(x,plumemsk) 
    aaim = np.ma.masked_array(aainew,plumemsk)    
    aodm = np.ma.masked_array(aodnew,plumemsk)    
    aodstdm = np.ma.masked_array(aodstdnew,plumemsk)    
    szam = np.ma.masked_array(szanew,plumemsk)    
    saam = np.ma.masked_array(saanew,plumemsk)    
    vzam = np.ma.masked_array(vzanew,plumemsk)    
    vaam = np.ma.masked_array(vaanew,plumemsk)    
    Psm = np.ma.masked_array(Psnew,plumemsk)    
    Hsm = np.ma.masked_array(Hsnew,plumemsk)    
    As1m = np.ma.masked_array(As1new,plumemsk)   
    As2m = np.ma.masked_array(As2new,plumemsk)  
    ssam = np.ma.masked_array(ssanew,plumemsk)       


    latrm = np.array(latm[latm.mask == 0])    
    lonrm = np.array(lonm[lonm.mask == 0])  
    aairm = np.array(aaim[aaim.mask == 0])
    aodrm = np.array(aodm[aodm.mask == 0])
    aodstdrm = np.array(aodstdm[aodstdm.mask == 0])
    szarm = np.array(szam[szam.mask == 0])
    saarm = np.array(saam[saam.mask == 0])
    vzarm = np.array(vzam[vzam.mask == 0])
    vaarm = np.array(vaam[vaam.mask == 0])
    Psrm = np.array(Psm[Psm.mask == 0])
    Hsrm = np.array(Hsm[Hsm.mask == 0])
    As1rm = np.array(As1m[As1m.mask == 0])
    As2rm = np.array(As2m[As2m.mask == 0])
    ssarm = np.array(ssam[ssam.mask == 0])
    
    
    Asrm = np.column_stack([As1rm, As2rm])
    
    return latrm, lonrm, aairm, aodrm, aodstdrm, szarm, saarm, vzarm, vaarm, ssarm, Psrm, plumemsk, Asrm
    
