# -*- coding: utf-8 -*-
"""
Run PYDISAMAR
apply OMI level2 data 

@author: sunj
"""


######### load python packages
import sys, os
import shutil
import time
import numpy as np
#from AERONET import AERONET_dailymean
#import numpy.ma as ma
#import matplotlib.pyplot as plt
#import netCDF4 
#import glob
#import tables
#import math
#from scipy import ndimage
#from GOME2 import rawGOME2, gridGOME2 
#from OMAERO import rawOMAERO, gridOMAERO, daytsOMAERO, OMAERO4cfg
#from MODIS import gridMODIS, daytsMODIS, MODIS4cfg 
#from CAMS import *  




######### directory
maindir = '/usr/people/sunj/Documents/DISAMAR/'
pydir  = maindir + 'pydisamar/src/'
exedir = maindir + 'disamar/Disamar.exe'
cfgdir = maindir + 'cfgfile/'
outputdir = maindir + 'pydisamar/output/'


#outputdir = maindir + 'pydisamar/output/'
#nobackupdir = '/nobackup/users/sunj/'
#
#exedir = maindir + 'disamar/createExpCoefFiles/'
#expindir = exedir + 'inputFiles/'
#expoutdir = exedir + 'expCoefFiles/'
#os.chdir(exedir)
#
#expcoefdir = maindir + 'disamar/expCoefFiles/'
#template = expindir + 'smoke_template.in'




#readTM5dir = '/nobackup/users/sunj/aai_sim_sunj/src/'
#TM5datadir = '/nobackup/users/sunj/TM5/AP3-INSITU-AAI_LT/'
#OMIdir = '/nobackup/users/sunj/OMI/'
#AER_AI_lut = '/data/project/tropomi-l2/data/LUT_AAI___/S5P_OPER_LUT_AAI____00000000T000000_99999999T999999_20150605T114510.nc'


######### load packages
sys.path.insert( 0, pydir)
import rt_cfg, rt_run

#sys.path.insert( 0, readTM5dir)
#import TM5_Class_sunj
#sys.path.insert( 0, maindir + 'pydisamar/')
#import AAI



######## prepare to run DISAMAR
def rundisamar(expname,year,month,day,cfgtemplate, simulationOnly, Gaussdp,atminterval,aerlev, ps, SZA, SAA, VZA, VAA, wvl, AOD, As, aerosolType, **other):
    
 
#    if len(AOD.shape) == 1: 
#        AOD = AOD.reshape([AOD.shape[0],1])    
    
    
    
    ######## modify configuration file 
    cfg = rt_cfg.RT_configuration(filename = cfgdir + cfgtemplate)    # your base config file
    
    
    
    
    
    
    # SECTION PRESSURE_TEMPERATURE
    cfg['GENERAL','overall','simulationOnly'].setvalue([simulationOnly])

    if aerlev == 0: 
#        for ilev in np.arange(0,AOD.shape[-1]-2):
        cfg['GENERAL','specifyFitting','numIntervalFit'].setvalue([2])         # +2 because python counts from 0
    else: 
        cfg['GENERAL','specifyFitting','numIntervalFit'].setvalue([aerlev])
        
        
    # SECTION GEOMETRY
    cfg['GEOMETRY','geometry','solar_zenith_angle_sim'].setvalue([SZA])
    cfg['GEOMETRY','geometry','solar_zenith_angle_retr'].setvalue([SZA])
    cfg['GEOMETRY','geometry','solar_azimuth_angle_sim'].setvalue([SAA])
    cfg['GEOMETRY','geometry','solar_azimuth_angle_retr'].setvalue([SAA])
    cfg['GEOMETRY','geometry','instrument_nadir_angle_sim'].setvalue([VZA])
    cfg['GEOMETRY','geometry','instrument_nadir_angle_retr'].setvalue([VZA])
    cfg['GEOMETRY','geometry','instrument_azimuth_angle_sim'].setvalue([VAA])
    cfg['GEOMETRY','geometry','instrument_azimuth_angle_retr'].setvalue([VAA])
    
    # SECTION SURFACE 
    cfg['SURFACE','pressure','surfPressureSim'].setvalue([ps])
    cfg['SURFACE','pressure','surfPressureRetr'].setvalue([ps])    
    ######## wavelength dependent
    if isinstance(wvl,list) == False:
        wvl = [wvl]
        As  = [As]
        
    cfg['SURFACE','surfaceType','surfaceTypeSim'].setvalue(['wavelDependent'])
    cfg['SURFACE','surfaceType','surfaceTypeRetr'].setvalue(['wavelDependent'])
    
    cfg['SURFACE','wavelDependentSim','wavelSurfAlbedo'].setvalue(wvl)
    cfg['SURFACE','wavelDependentSim','surfAlbedo'].setvalue(As)     
    cfg['SURFACE','wavelDependentSim','wavelSurfEmission'].setvalue(wvl)
    cfg['SURFACE','wavelDependentSim','surfEmission'].setvalue(As)     

    cfg['SURFACE','wavelDependentRetr','wavelSurfAlbedo'].setvalue(wvl)
    cfg['SURFACE','wavelDependentRetr','surfAlbedo'].setvalue(As)     
    cfg['SURFACE','wavelDependentRetr','wavelSurfEmission'].setvalue(wvl)
    cfg['SURFACE','wavelDependentRetr','surfEmission'].setvalue(As)  
    

    
    
    # SECTION ATMOSPHERIC_INTERVALS
    cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','topPressureSim'].setvalue(list(atminterval))
    cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','topPressureRetr'].setvalue(list(atminterval))
    item = cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','APvarianceTopPressure']
    item.setvalue((len(atminterval)-1)*[25e4]+[1e-4])


    # SECTION RADIATIVE TRANSFER 
    cfg['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltSim'].setvalue([Gaussdp]*len(atminterval))
    cfg['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltRetr'].setvalue([Gaussdp]*len(atminterval))
    cfg['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltLaySim'].setvalue([Gaussdp])
    cfg['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltLayRetr'].setvalue([Gaussdp])
    cfg['RADIATIVE_TRANSFER','RTM_Sim_Retr','nstreamsSim'].setvalue([28])
    cfg['RADIATIVE_TRANSFER','RTM_Sim_Retr','nstreamsRetr'].setvalue([28])


    # SECTION O3
    cfg['O3','climatology','latitude'].setvalue([25.0])
    cfg['O3','climatology','month'].setvalue([month])
    cfg['O3','climatology','day_of_month'].setvalue([day])



    # SECTION AEROSOL
    ######## simulation
    if aerosolType == 'none': 
        # other = {ANG , SSA , g}
    
        cfg['AEROSOL','aerosolType','aerosolTypeSim'].setvalue(['none'])
        cfg['AEROSOL','aerosolType','aerosolTypeRetr'].setvalue(['none'])

    if aerosolType == 'HG': 
        # other = {ANG , SSA , g}
    
        cfg['AEROSOL','aerosolType','aerosolTypeSim'].setvalue(['HGscattering'])
        cfg['AEROSOL','aerosolType','aerosolTypeRetr'].setvalue(['HGscattering'])
        
        ANG = other['ANG']
        SSA = other['SSA']
        g = other['g']
        
#        if aerlev == 0: 
            
        AODlev = [] 
        ANGlev = []
        SSAlev = []
        glev = [] 
        
        for ilev in np.arange(2,AOD.shape[-1]):
            AODlev.append([ilev,AOD[ilev-1]])
            ANGlev.append([ilev,ANG[ilev-1]])
            SSAlev.append([ilev,SSA[ilev-1]])
            glev.append([ilev,g[ilev-1]])
            
        item = cfg['AEROSOL','HGscatteringSim','opticalThickness']
        item.set_rawvalue(AODlev)
        item = cfg['AEROSOL','HGscatteringSim','angstromCoefficient']
        item.set_rawvalue(ANGlev)
        item = cfg['AEROSOL','HGscatteringSim','singleScatteringAlbedo']
        item.set_rawvalue(SSAlev)
        item = cfg['AEROSOL','HGscatteringSim','HGparameter_g']
        item.set_rawvalue(glev)
    

#            item = cfg['AEROSOL','HGscatteringRetr','opticalThickness']
#            item.set_rawvalue(AODlev)
#            item = cfg['AEROSOL','HGscatteringRetr','angstromCoefficient']
#            item.set_rawvalue(ANGlev)
#            item = cfg['AEROSOL','HGscatteringRetr','singleScatteringAlbedo']
#            item.set_rawvalue(SSAlev)
#            item = cfg['AEROSOL','HGscatteringRetr','HGparameter_g']
#            item.set_rawvalue(glev)

        
#        else: 
#            print 'aerosol level:', aerlev
#            AODlev = [aerlev, AOD]
#            ANGlev = [aerlev, ANG]
#            SSAlev = [aerlev, SSA]
#            glev = [aerlev, g] 
#            
#            cfg['AEROSOL','HGscatteringSim','opticalThickness'].setvalue([aerlev, AOD])
#            cfg['AEROSOL','HGscatteringSim','angstromCoefficient'].setvalue([aerlev, ANG])
#            cfg['AEROSOL','HGscatteringSim','singleScatteringAlbedo'].setvalue([aerlev, SSA])
#            cfg['AEROSOL','HGscatteringSim','HGparameter_g'].setvalue([aerlev, g])
            

#            cfg['AEROSOL','HGscatteringRetr','opticalThickness'].setvalue([aerlev, AOD])
#            cfg['AEROSOL','HGscatteringRetr','angstromCoefficient'].setvalue([aerlev, ANG])
#            cfg['AEROSOL','HGscatteringRetr','singleScatteringAlbedo'].setvalue([aerlev, SSA])
#            cfg['AEROSOL','HGscatteringRetr','HGparameter_g'].setvalue([aerlev, g])

       
        
    if aerosolType == 'Mie':
        # other = {expcoef}    
        expcoefFile = other['expcoef']
        
        cfg['AEROSOL','aerosolType','aerosolTypeSim'].setvalue(['MieScattering'])
#        cfg['AEROSOL','aerosolType','aerosolTypeRetr'].setvalue(['MieScattering'])
        
#        if aerlev == 0: 
        AODlev = [] 
    
        for ilev in np.arange(2,AOD.shape[-1]):
            AODlev.append([ilev,AOD[ilev-1],'expCoefFiles/' + expcoefFile])
#                AODlev.append([ilev+2,AOD[ilev],'expCoefFiles/strong_abs_320-400.dat'])
            
        item = cfg['AEROSOL','MieScatteringSim','aerosolOpticalThickness']
        item.set_rawvalue(AODlev)            

#        else:
    #        for ilev in np.arange(1,AOD.shape[0]):
#            item = cfg['AEROSOL','MieScatteringSim','aerosolOpticalThickness'] 
    #        item.set_rawvalue([[intervalnr,AOD,'expCoefFiles/strong_abs_320-400.dat']])
    
    #        item.set_rawvalue([[intervalnr, AOD ,'expCoefFiles/strong_abs_320-400.dat'],\
    #                           [intervalnr+1, AOD ,'expCoefFiles/strong_abs_320-400.dat']])
                    
    #        item.set_rawvalue([[aerlev, AOD ,'expCoefFiles/strong_abs_320-400.dat']])
#            item.set_rawvalue([[aerlev, AOD ,'expCoefFiles/' + expcoefFile]])
#            item.set_rawvalue([[aerlev, AOD ,'expCoefFiles/strong_abs_320-400.dat']])

                            

    
    cfg.setfile(cfgdir + expname + '.in')
    cfg.write()
 
    ######### run DISAMAR
    if os.path.exists(outputdir + expname) == True: 
        pass 
    
    if os.path.exists(outputdir + expname) == False: 
        os.mkdir(outputdir + expname)
    outputname = outputdir + expname + '/disamar.h5'
    
    run = rt_run.rt_run( cfg=cfg , disamar=exedir , output=outputname , spectrum=None , debug=False , quiet=False , tempbase=None)
    run() 
    
    return outputname 


    



######### read parameters from TM5, for now use slab average 
#def tm54cfg(year, month, day, atminterval):
#    model_name = 'tm5'
#    model_nlayers = 34
#    
#    ######## get TM5 outputs 2 deg by 3 deg 
#    TM5 = TM5_Class_sunj.TM5_Class(year=year, month=month, day=day, model='tm5')
#    
#    # read TM5 output
#    aod550 = TM5.ec2aod(TM5.ec5504D)
#    aod350 = TM5.ec2aod(TM5.ec3504D)      
#    ssa550 = TM5.ssa5504D
#    g550 = TM5.g5504D    
#    ps = TM5.ps
#    A550 = TM5.compute_Angexp(350,550)
#    tm5lat = TM5.lat
#    tm5lon = TM5.lon
#
#    
#    # in case of zero denominator
#    aod550[aod550==0] += 1e-8
#    aod350[aod350==0] += 1e-8
#    ssa550[ssa550==0] += 1e-8
#    
#    
#    # initialization
#    aod550int = np.zeros([7,len(atminterval),aod550.shape[2],aod550.shape[3]])
#    aod350int = np.zeros(aod550int.shape)
#    ssa550int = np.zeros(aod550int.shape)
#    g550int = np.zeros(aod550int.shape)
#    A550int = np.zeros(aod550int.shape)
#
#
#    ######## use volume weighed AOD, volume and AOD weighed SSA, volume, AOD and SSA weighed g 
#    for imode in range(7):
#        for iatm, ipres in enumerate(atminterval):
#            for i in range(aod550.shape[2]):
#                for j in range(aod550.shape[3]):
#                    
#                    # if DISAMAR interval is outside TM5, use the value at the lowest/ highest level
#                    if ipres > TM5.preslay[:,i,j].max():
#                        idx1 = 0
#                        idx2 = 0
#                        w1 = w2 = 0.5
#        
#                    if ipres < TM5.preslay[:,i,j].min():
#                        idx1 = -1 
#                        idx2 = -1
#                        w1 = w2 = 0.5
#                        
#                    # if DISAMAR interval is inside TM5, use the pressure /volume calculate weight 
#                    if (ipres > TM5.preslay[:,i,j].min()) and (ipres < TM5.preslay[:,i,j].max()) :
#                        idx1 = np.where(TM5.preslay[:,i,j] > ipres)[0][-1]
#                        idx2 = np.where(TM5.preslay[:,i,j] < ipres)[0][0]  
#                        p1 = TM5.preslay[idx1,i,j]
#                        p2 = TM5.preslay[idx2,i,j]    
#
#                        w1 = 1./(p1-ipres) / (1./(p1-ipres) + 1./(ipres-p2))
#                        w2 = 1./(ipres-p2) / (1./(p1-ipres) + 1./(ipres-p2))
#                        
#                    # volume weighed AOD 
#                    aod550int[imode,iatm,i,j] = w1 * aod550[imode,idx1,i,j] + w2 * aod550[imode,idx2,i,j]
#                    aod350int[imode,iatm,i,j] = w1 * aod350[imode,idx1,i,j] + w2 * aod350[imode,idx2,i,j]
#                    
#                    # volume and AOD weighed ssa
#                    ssa550int[imode,iatm,i,j] = (w1 * aod550[imode,idx1,i,j] * ssa550[imode,idx1,i,j] + \
#                                             w2 * aod550[imode,idx2,i,j] * ssa550[imode,idx2,i,j]) / aod550int[imode,iatm,i,j] 
#                                     
#                    
#                    # volume, AOD and ssa weighed g 
#                    g550int[imode,iatm,i,j] = (w1 * aod550[imode,idx1,i,j] * ssa550[imode,idx1,i,j] * g550[imode,idx1,i,j] + \
#                                           w2 * aod550[imode,idx2,i,j] * ssa550[imode,idx2,i,j] * g550[imode,idx2,i,j]) \
#                                           / ssa550int[imode,iatm,i,j] / aod550int[imode,iatm,i,j]                    
#                    
#                    # weighed AOD calculate Angstrom exponent 
#                    A550int[imode,iatm,i,j] = np.log(aod550int[imode,iatm,i,j]/aod350int[imode,iatm,i,j]) / np.log(350./550.)
#    
#
#
#
#    ######## mode combination 
#    # select modes: accs, accui, coas, coai 
#    aod5504m = np.delete(aod550int,[0,1,4],0)
#    aod3504m = np.delete(aod350int,[0,1,4],0)
#    ssa5504m = np.delete(ssa550int,[0,1,4],0)
#    g5504m = np.delete(g550int,[0,1,4],0)
#
#    
#    # combine mode 
#    # sum up 4 modes into total AOD
#    aod550tot = np.sum(aod5504m,0)
#    aod350tot = np.sum(aod3504m,0)
#    
#    # total AOD weighed ssa 
#    ssa550tot = np.sum(aod5504m*ssa5504m,0) / aod550tot
#    
#    # total AOD, ssa weighed g 
#    g550tot = np.sum(aod5504m*ssa5504m*g5504m,0) / (aod550tot*ssa550tot)
#    
#    # total AOD calculate Angstrom exponent 
#    A550tot = np.log(aod550tot / aod350tot) / np.log(350./550.)
#
#
#    return tm5lat, tm5lon, aod550tot, ssa550tot, g550tot, A550tot, ps
#
#    
#    
#def omi4cfg(year, month, day, orbit):
#    ######## get OMI L2 data 0.25 deg by 0.25 deg OMAERUV /OMEARO  
#    
#    # saa, vaa are avaialble in OMAERO
#    OMIL2dir = '/nobackup/users/sunj/OMI/OMAERO/%4i/%02i/%02i/' %(year, month, day)
#    OMIL2orb = glob.glob( OMIL2dir + '*.he5')
#
#    # saa, vaa, O3 is avaialble in OMDOAO3
#    OMIL2dir3 = '/nobackup/users/sunj/OMI/OMDOAO3/%4i%02i%02i/' %(year, month, day)
#    OMIL2orb3 = glob.glob( OMIL2dir3 + '*.he5')
#    
##    # saa, vaa, O3 is avaialble in OMTO3
##    OMIL2dir4 = '/nobackup/users/sunj/OMI/OMTO3/'
##    OMIL2orb4 = glob.glob( OMIL2dir4 + '*.he5')
#
#
#    
#    for io in OMIL2orb[:]: 
#
#        if io.find(str(orbit)) != -1:  # find the orbit number
#            OMIL2 = netCDF4.Dataset(io,'r')    
#    
#            OMIL2data = OMIL2.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Data Fields']
#            OMIL2geo = OMIL2.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Geolocation Fields']
#    
#            lat = OMIL2geo.variables['Latitude'][:]
#            lon = OMIL2geo.variables['Longitude'][:]       
#            sza = OMIL2geo.variables['SolarZenithAngle'][:]
#            vza = OMIL2geo.variables['ViewingZenithAngle'][:]
#            saa = OMIL2geo.variables['SolarAzimuthAngle'][:]
#            vaa = OMIL2geo.variables['ViewingAzimuthAngle'][:]
#            surfh = OMIL2geo.variables['TerrainHeight'][:]
#            alb = OMIL2data.variables['TerrainReflectivity'][:]*1e-3                     # alb [lat, lon, wvl]     
#            aodomi = OMIL2data.variables['AerosolOpticalThicknessMW'][:]*1e-3       # aod [lat, lon, wvl] 
#            aaiomi = OMIL2data.variables['AerosolIndexUV'][:]*1e-2
##            wvl = OMIL2data.variables['Wavelength'][:]                        # wavelength [354, 388, 500]
#            
#    
#
##    for io in OMIL2orb2[:]:
##        if io.find(str(orbit)) != -1: 
##            OMIL2 = netCDF4.Dataset(io,'r')    
##            OMIL2geo2 = OMIL2.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountAerosol'].groups['Geolocation Fields']
##            saa = OMIL2geo2.variables['SolarAzimuthAngle'][:]
##            vaa = OMIL2geo2.variables['ViewingAzimuthAngle'][:]
#
#
#    for io in OMIL2orb3[:]:
#        if io.find(str(orbit)) != -1: 
#            OMIL2 = netCDF4.Dataset(io,'r')    
#            OMIL2data2 = OMIL2.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountO3'].groups['Data Fields']
#            OMIL2geo2 = OMIL2.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountO3'].groups['Geolocation Fields']
#            vcd_o3 = OMIL2data2.variables['ColumnAmountO3'][:]
#
#    
#    return lat,lon, sza, vza, saa, vaa, alb, aodomi, aaiomi, surfh, vcd_o3
#
#
######### interpolate TM5 onto OMI L2 orbit
#def interp2omi(lat, lon, tm5lat,tm5lon, aod550tot, ssa550tot, g550tot, A550tot, ps):     
#              
#    model_layer = aod550tot.shape[0]
#    
#    aod550_grid = resampletm5(lat,lon,tm5lat, tm5lon, aod550tot,model_layer) 
#    ssa550_grid = resampletm5(lat,lon,tm5lat,tm5lon,ssa550tot,model_layer)  
#    g550_grid = resampletm5(lat,lon,tm5lat, tm5lon, g550tot,model_layer)  
#    A550_grid = resampletm5(lat,lon,tm5lat, tm5lon, A550tot,model_layer)
#    ps_grid = resampletm5(lat,lon,tm5lat,tm5lon,ps,1)
#    
#
#    ######### for checking
#    mask = np.logical_or(np.isnan(aaiomi), aaiomi.data<0)
#    mask = np.logical_or(mask, aodomi.data[:,:,-1]<0)
#    mask = np.logical_or(mask, sza>75.)
#    
#    lon = np.ma.masked_array(lon,mask)
#    lat = np.ma.masked_array(lat,mask)
#    
#    aodomi[:,:,-1] = np.ma.masked_array(aodomi.data[:,:,-1],mask)
#    aod550_grid[0,:,:] = np.ma.masked_array(aod550_grid[0,:,:],mask)
#
#    
#
#            
#    return aod550_grid, ssa550_grid, g550_grid, A550_grid, ps_grid
#
#
#
#def resampletm5(satlat, satlon, tm5lat, tm5lon, data, nlayers, order=1, prefilter=False):
#    # compute the coordinates
#    coords_lat = np.interp( satlat.flatten(), tm5lat, np.arange(tm5lat.shape[0]) )
#    coords_lon = np.interp( satlon.flatten(), tm5lon, np.arange(tm5lon.shape[0]) )	
#    coords = np.array( [coords_lat[:], coords_lon[:]] )
#
#    # for 3D interpolation
#    if nlayers > 1: 
#        res = np.zeros( (nlayers, lat.shape[0], lat.shape[1]) )        
#        for ilayer in range(nlayers):
#            res[ilayer,:,:] = np.reshape( ndimage.map_coordinates(data[ilayer,:,:], coords, order=order, prefilter=prefilter), (1, res.shape[1], res.shape[2]) )
#
#    # for 2D interpolation 
#    if nlayers == 1: 
#        res = np.zeros( ( lat.shape[0], lat.shape[1]) )        
#        for ilayer in range(nlayers):
#            res[:,:] = np.reshape( ndimage.map_coordinates(data[:,:], coords, order=order, prefilter=prefilter), (res.shape[0], res.shape[1]) )
#    return res
#    
#


