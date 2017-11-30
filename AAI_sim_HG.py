# -*- coding: utf-8 -*-
"""
Run PYDISAMAR AAI simulation HG scheme.

@author: sunj
"""


######### load python packages
import os
import shutil
import time
import numpy as np
import numpy.ma as ma
import tables
from scipy.interpolate import griddata 
from DISAMAR import rundisamar
from CALIPSO import proCALIOP 
from OMI import gridOMAERO, OMAERO4cfg
from MODIS import gridMODIS, MODIS4cfg
from GOME2 import gridGOME2
from AERONET import AERONET_dailymean
from timezone import timezone

import matplotlib.pyplot as plt
maindir = '/usr/people/sunj/Documents/DISAMAR'
outputdir = maindir + '/pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'

maindir = '/usr/people/sunj/Documents/DISAMAR/disamar'
exedir = maindir + '/createExpCoefFiles'
expindir = exedir + '/inputFiles/'
expoutdir = exedir + '/expCoefFiles/'
scaoutdir = exedir +'/scaMatFiles/'
expcoefdir = maindir + '/expCoefFiles/'

LUTdir = 'AAI_LUT/'



# initialization 
#ssaval = np.arange(0.7,0.96,0.05)             # single scattering albedo
#ssaval = [0.75]
#
#Angval = np.arange(0.5,2.1,0.5)               # angstrom exponent
#Angval = [2.0]
#
#asyval = np.arange(0.5,1.0,0.1)                       # asymmetry factor
#asyval = [0.9]

# initialization 
# construct aerosol profile 
L = -0.0065         # lapse rate, unit: K/m
R = 8.31432         # gas constant 
g0 = 9.80665        # gravitational accelaration constant
M = 0.0289644       # air molar mass, unit: kg/mol
T0 = 15.+ 273.15     # surface temperature, unit: K 
P0 = 1013.     # surface pressure, unit: hPa
z0 = 0.     




aaiomi = []
aaidsm = []
aaiAAI = [] 

start=time.time() 

years = range(2017,2018)
months = range(1,2)
days = range(26,30)
#days.remove(19)
#days.remove(21)

days = [27]

aerlat = -33.46
aerlon = -70.66

# calculate time difference
jetlag = timezone(aerlon)
wvl = np.arange(340,677,2)


ROI = [-40,-5,-105,-70]

crival = 1
res = 0.5

nlay = 4
t0 = time.time()


# Find the minimum RMS
#alhval = [5,3,2,4] 
LUTfile = 'LUT_ALH_ALT_%4i%02i%02i-%4i%02i%02i.h5' % (years[0], months[0], days[0], years[-1], months[-1], days[-1]) 
LUT = tables.open_file(nobackupdir + LUTdir + LUTfile, 'r').root
rms = LUT.Root_mean_suquared_error[:]
aaimean = LUT.DSM_AAI_mean[:]
aaimedian = LUT.DSM_AAI_median[:]
alh = LUT.Aerosol_layer_height[:]
alt = LUT.Aerosol_layer_thickness[:]


v = np.arange(0, np.max(rms)+0.5, 0.01)
for dd, iday in enumerate(days):
    plt.figure(figsize=(5,4), frameon = False)
    plt.contourf(alt,alh, rms[dd,:,:], v,cmap = 'rainbow')
    plt.xlabel('ALT')
    plt.ylabel('ALH')
    plt.xticks(alt)
    plt.yticks(alh)
    cbar = plt.colorbar()
    cbar.set_label('RMS')
    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    plt.title('2017-01-%02i' %(iday))
    
#    plt.savefig(figdir + 'AAI_sim_Mie_stats_%4i%02i%02i_rms.png' % (iyear, imon, iday), dpi = 300, transparent = False) 



#alhval = np.zeros(len(days))
#altval = np.zeros(len(days))
#
#
#for dd, iday in enumerate(days):
#    idx = np.argmin(rms[dd,:,:])
#    ih = idx / len(alh)
#    it = idx % len(alt)
#    alhval[dd] = alh[ih]
#    altval[dd] = alt[it] 
#    print alh[ih], alt[it]
#    
##    print rms[dd,ih,it]
#
#print alhval, altval 
#    
#
#
#hh, tt = np.meshgrid(alh, alt)
#hh = hh.reshape(1,hh.size)[0]
#tt = tt.reshape(1,tt.size)[0]
#rms =rms.reshape(1, rms.size)[0]
#
#A = np.vstack([hh,tt,np.ones(rms.size)]).T
#a,b,c = np.linalg.lstsq(A, rms)[0]
#
#
#
#testname = 'AAI_HG_based_on_Mie/'
#nobackupdir += testname 


#for iyear in years:
#    for imon in months:
#        for dd, iday in enumerate(days): 
#            
#            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res)
##            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
##            latomi, lonomi, aaiomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AAI', res) 
##            latomi, lonomi, aodomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AOD', res) 
#            plumemsk = np.logical_or(aodmod0 < 0., aodmod0 < 0.)
##            plumemsk = np.logical_or(aodomi0 < 0., aodmod0 < 0.)
##            plumemsk = np.logical_or(plumemsk, aaiomi0 < 0.)
##            plumemsk = np.logical_or(plumemsk, aaiomi0 < crival)     
##            plumemsk = np.logical_or(plumemsk, aaigome0 < crival) 
#            
##            aaiomi = np.ma.masked_array(aaiomi0,plumemsk)   
##            aaigome = np.ma.masked_array(aaigome0,plumemsk)   
##            aodomi = np.ma.masked_array(aodomi0,plumemsk)    
##            aodmod = np.ma.masked_array(aodmod0,plumemsk)
#            
#            
#            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
#            lat, lon, aodm, Ang, asy = MODIS4cfg(iyear, imon, iday, ROI, plumemskomi, res)
#            aod550 = AERONET_dailymean(iyear,imon,iday, 'aod')
#            Angval = [AERONET_dailymean(iyear,imon,iday, 'AE')]
#
#            aodscl = aod550 / aodm.mean()
#            aodms = aodm * aodscl   
#            
#
#
#            
##            t1 = time.time()
##            print '#### time for reading data:', t1 - t0, 's'
#            
#            for imode in ['bimodal']:
#                expfilename = 'AAI_sim_Mie_%s_%04i%02i%02i' % (imode, iyear, imon, iday)
#                expfile = expcoefdir + expfilename + '.dat'
#                with open(expfile,'r') as f: 
#                    content = f.readlines()
#                
#                ssamie = [] 
#                for il in content: 
#                    if il.find('ssa') != -1:
#                        ssamie.append(float(il[5:11]))
#                ssamie = np.array(ssamie)
#                ssaval = [ssamie[np.where(wvl == 550)[0]]]
#                
##                
##                
##                expfilename = 'AAI_sim_Mie_%s_%04i%02i%02i' % (imode, iyear, imon, iday)
##                expfile = expoutdir + expfilename + '.dat'
##                with open(expfile,'r') as f: 
##                    content = f.readlines()
##                
#                gmie = [] 
#                for ig, il in enumerate(content): 
#                    if il.find('alpha1') != -1:
#                        temp1 = content[ig+2]
#                        temp2 = float(temp1.split('  ')[2]) / 3.
#                        gmie.append(temp2)  
#                asyval = [gmie[np.where(wvl == 550)[0]]]
#                
#                print ssaval, Angval, asyval 
#                
#                
#                
#        
#                for ip in range(len(lat)):
#    #            for ip in range(26,27):
#                    print ip, '/', len(lat)
##                    DSMpres, DSMnormpro, DSMpresg, DSMnormprom  = proCALIOP(iyear, imon, iday, ROI, lat[ip], lon[ip], nlay)
##                        
##                    
##                    if DSMpres[0] > 1011.:
##                        diff = DSMpres[0] - 1009.
##                        DSMpres = DSMpres - diff  
##                    
##                    DSMinterval = [1011.] + list(DSMpresg) + [0.3]
#                    
#    #                for il in range(0,len(DSMinterval)-1):
#    #                    print il 
#    #                    if DSMinterval[il+1] / DSMinterval[il] > 0.999:
#    #                        DSMinterval[il+1] = DSMinterval[il+1] - 2.
#                    
#                    
##                    DSMpro = np.array([1.0e-4]+list(DSMnormprom*aodms[ip])+[1.0e-4])
#                    ialh = alhval[dd]
#                    ialt = altval[dd]
#                    z = np.array([40] + [10] + [ialh + ialt / 2.] + [ialh - ialt / 2.])  # assume aerosol layer thickness is 1 km 
#                    z = z[::-1]
#                    DSMpres = P0 * (1+L/T0 * (z*1e3-z0)) ** (-g0*M/R/L)     # unit: hPa
#                    DSMpro = np.array([1.0e-4]+[aodms[ip]]+[1.0e-4])                    
#                    
##                    print DSMpres 
#                    aerlev = np.where(DSMpro == DSMpro.max())[0][0]+1
#                    if aerlev == 1: 
#                        aerlev = 2 
#    
#    
#                     
##                    t2 = time.time()
#    #                print '#### time for preparing profile:', t2 - t1, 's'
#                    
#                    for issa in ssaval: 
#                        for iAng in Angval:
#                            for iasy in asyval:
##                                t3 = time.time()
#                                
#                                
#                                expname = 'AAI_sim_HG_%4i%02i%02i' %(iyear, imon, iday)
#                                ######## Test HG: all test
#                                outputname = rundisamar(expname = expname, year=iyear, month=imon, day=iday, cfgtemplate = 'Config_AAI_354_388_100hPa_SA_HG.in',\
#                                            simulationOnly = 0, Gaussdp = 8, atminterval = DSMpres, aerlev = aerlev,  ps = 1013., SZA = sza[ip], SAA = saa[ip], \
#                                            VZA = vza[ip], VAA = vaa[ip], wvl = [354., 388.], AOD = DSMpro, As = list(alb[ip,:]), \
#                                            aerosolType = 'HG', ANG = np.ones(len(DSMpres))*iAng, SSA = np.ones(len(DSMpres))*issa, g = np.ones(len(DSMpres))*iasy)
#                        
#                                src = '/AAI_sim_HG_lat_%2.2f_lon_%3.2f_saa-%1.2f_Ang-%1.1f_asy-%1.1f.h5' % (lat[ip], lon[ip], issa, iAng, iasy)
#                                print '**** output file:',src 
#                                
#                                
#                                if not os.path.exists(nobackupdir + expname):
#                                    os.makedirs(nobackupdir + expname)   
#                                
#                                os.rename(outputdir + expname + '/disamar.h5', outputdir + expname + src)
#         
#                                if os.path.isfile(nobackupdir + expname + src):
#                                    os.remove(nobackupdir + expname + src)
#        
#        
#                                shutil.move(outputdir + expname + src, nobackupdir + expname + '/')       
#                                
##                                t4 = time.time()
##                                print '#### time for running DISAMAR per simulation:', t4 - t3, 's'
#            
#        
#                                ######### call AAI module 
#                                disamarhdf = tables.open_file(nobackupdir + expname + src)
#                                
#                                disamarout = disamarhdf.root                    
#                                aairetr = disamarout.parameters.AAI[0]
#                                disamarhdf.close()
#                                
#                                print '**** AAI retrieved: %1.4f'% aairetr 
#                    print '**** OMI AAI: %1.4f' % aai[ip]
##                t5 = time.time()
##                print '#### total time for running DISAMAR per day:', t5 - t1, 's'
#            
#t6 = time.time() 
#print '#### total time for running DISAMAR:', (t6 - t0) / 3600., 'h'
