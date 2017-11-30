# -*- coding: utf-8 -*-
"""
Run AAI simulation Mie scheme

@author: sunj
"""


######### load python packages
import os
import shutil
import time
import numpy as np
import numpy.ma as ma
import tables 
from DISAMAR import rundisamar
from OMI import gridOMAERO, OMAERO4cfg
from MODIS import gridMODIS, MODIS4cfg
from GOME2 import gridGOME2
from AERONET import AERONET_day, AERONET_dailymean
from timezone import timezone
from DSM_Mie_func import create_Expcoef


# directory ===================================================================
maindir = '/usr/people/sunj/Documents/DISAMAR/'
outputdir = maindir + 'pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'

exedir = maindir + 'disamar/createExpCoefFiles/'
expindir = exedir + 'inputFiles/'
expoutdir = exedir + 'expCoefFiles/'
os.chdir(exedir)

expcoefdir = maindir + 'disamar/expCoefFiles/'
template = expindir + 'smoke_template.in'



# simulation initialization ===================================================
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

start = time.time() 

years = range(2017,2018)
months = range(10,11)
days = range(17,18)
#days.remove(28)
#days = [27]

aerlat = -33.46
aerlon = -70.66

# calculate time difference
jetlag = timezone(aerlon)

ROI = [45, 65, 0, 50]      # portugal 
ROIaer = [50, 65, 5, 45]
#ROI = [-40,-5,-105,-70]

crival = 2.
res = 0.5

nlay = 4

## testing aerosl layer height 1km, 2km, 3km, 4km, 5km
#alhval = np.arange(3,6.1,0.5)
#altval = np.arange(0.5,2.1,0.5)
#factorval = np.arange(1,4,1.)

alhval = np.arange(3,4)
altval = np.arange(2,3)
factorval = np.arange(1,4,1.)
aodval = np.arange(1,3,0.5)

t0 = time.time() 
for iyear in years:
    for imon in months:
        for iday in days: 
           
#            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res)
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
            latomi, lonomi, aaiomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AAI',  res) 
#            plumemsk = np.logical_or(aodmod0 < 0., aodmod0 < 0.)
            
            plumemsk = np.logical_or(np.isnan(aaiomi0), np.isnan(aaigome0))
            plumemsk = np.logical_or(plumemsk, aaigome0 < crival)
            plumemsk = np.logical_or(plumemsk, aaiomi0 < crival)

            
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
#            lat, lon, aodm, Ang, asy = MODIS4cfg(iyear, imon, iday, ROI, plumemskomi, res)
        
        
            aerosite = AERONET_day(iyear, imon, iday, ROIaer)
            

#                aodms = aodm * 1. 
                
            for ifactor in factorval: 
                fine, coarse, bimodal = create_Expcoef(iyear, imon, iday, aerosite, ifactor)    
       
                for iaod in aodval: 
                    aodms = iaod * aai / aai.mean()             
                    
                    t1 = time.time()
                    print '#### time for reading data:', t1 - t0, 's'
                    
                    for ialh in alhval: 
                        for ialt in altval: 
                            if ialh > ialt:
                                for ip in range(len(lat)): 
                                    print ip, '/', len(lat)
                    
                                    z = np.array([40] + [ialh + ialt / 2.] + [ialh - ialt / 2.])    # assume aerosol layer thickness is 1 km 
                                    z = z[::-1]
                                    DSMpres = P0 * (1+L/T0 * (z*1e3-z0)) ** (-g0*M/R/L)     # unit: hPa
                                    DSMpro = np.array([1.0e-4]+[aodms[ip] / ialt]+[1.0e-4])
                                    
                                    
                                    aerlev = np.where(DSMpro == DSMpro.max())[0][0]+1
                                    if aerlev == 1: 
                                        aerlev = 2 
                                    
                                     
                                    for expname in [bimodal]:
                                    ###### DSM Mie         
                                        outputname = rundisamar(expname = expname, year=iyear, month=imon, day=iday, cfgtemplate = 'Config_AAI_354_388_100hPa_SA_Mie.in',\
                                                    simulationOnly = 0, Gaussdp = 20, atminterval = DSMpres, aerlev = aerlev,  ps = 1013., SZA = sza[ip], SAA = saa[ip], \
                                                    VZA = vza[ip], VAA = vaa[ip], wvl = [354., 388.], AOD = DSMpro, As = list(alb[ip,:]), \
                                                    aerosolType = 'Mie', expcoef = expname+'.dat')         
                        
                        
                                
                                        src = '/AAI_sim_Mie_lat_%2.2f_lon_%3.2f.h5' % (lat[ip], lon[ip])
                                        
                                        print src
#                                        testname = '/alh-%02i_alt-%1.2f_ni-%1.2f' % (ialh, ialt,ifactor)
                                        testname = '/aod-%1.2f_ni-%1.2f' % (iaod, ifactor)
                                        
                                        if not os.path.exists(nobackupdir + expname + testname):
                                            os.makedirs(nobackupdir + expname + testname)   
                                        
                                        os.rename(outputdir + expname + '/disamar.h5', outputdir + expname + src)
                         
                                        if os.path.isfile(nobackupdir + expname + testname + src):
                                            os.remove(nobackupdir + expname + testname + src)
                        
                        
                                        shutil.move(outputdir + expname + src, nobackupdir + expname + testname)       
                        
                                                          
                                        t2 = time.time()  
                                        print 'Time period of running DISAMAR per simulation:',t2-t1,'s'
                                        
                                        ######### call AAI module 
                                        disamarhdf = tables.open_file(nobackupdir + expname + testname + src)
                                        
                                        disamarout = disamarhdf.root                    
                                        aairetr = disamarout.parameters.AAI[0]
                                        disamarhdf.close()
                                        
                                        print expname
                                        print '**** AAI retrieved: %1.4f'% aairetr   
        
 
   
t3 = time.time()
print 'Time period of running DISAMAR:',t3-t0,'s'    
    
    
