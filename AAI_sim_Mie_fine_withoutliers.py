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
from AERONET import AERONET_dailymean
from timezone import timezone
from DSM_Mie_func import create_Expcoef

maindir = '/usr/people/sunj/Documents/DISAMAR/'
outputdir = maindir + 'pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'
LUTdir = 'AAI_LUT/'


exedir = maindir + 'disamar/createExpCoefFiles/'
expindir = exedir + 'inputFiles/'
expoutdir = exedir + 'expCoefFiles/'
os.chdir(exedir)

expcoefdir = maindir + 'disamar/expCoefFiles/'
template = expindir + 'smoke_template.in'



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


years = range(2017,2018)
months = range(1,2)
days = range(26,31)
days.remove(28)
#days = [27]

aerlat = -33.46
aerlon = -70.66

# calculate time difference
jetlag = timezone(aerlon)

ROI = [-40,-5,-105,-70]

crival = 1.
res = 0.5

nlay = 4


alhest = [] 
altest = [] 
factorest = []



# prepare sensitivity range 

LUTfile = 'coarse_LUT_%4i%02i%02i-%4i%02i%02i.h5' % (years[0], months[0], days[0], years[-1], months[-1], days[-1])     
LUT = tables.open_file(nobackupdir + LUTdir + LUTfile, 'r').root
rmse = LUT.Root_mean_suquared_error[:]
alh = LUT.Aerosol_layer_height[:]
alt = LUT.Aerosol_layer_thickness[:]
factor = LUT.Refractive_index_factor[:]


for il in range(len(rmse)): 
    ix, iy, iz = np.where(rmse[il,:,:,:] == np.nanmin(rmse[il,:,:,:]))
    factorest.append(factor[ix[0]])
    alhest.append(alh[iy[0]])
    altest.append(alt[iz[0]])



t0 = time.time() 
for iyear in years:
    for imon in months:
        for dd, iday in enumerate(days):
            
            factorval = np.array([factorest[dd] - 0.2, factorest[dd] - 0.1, factorest[dd], factorest[dd] + 0.1, factorest[dd] + 0.2])
            alhval = np.array([alhest[dd] - 0.5, alhest[dd] - 0.25, alhest[dd], alhest[dd] + 0.25, alhest[dd] + 0.5]) 
            altval = np.array([altest[dd] - 0.25,  altest[dd], altest[dd] + 0.25])
            
            print 'factor range: ', factorval 
            print 'alh range: ', alhval 
            print 'alt range: ', altval 
            
           
            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res, '3K')
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
            
            plumemsk = np.logical_or(aaigome0 < crival, aodmod0 < 0.)
            
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
            lat, lon, aodm, Ang, asy = MODIS4cfg(iyear, imon, iday, ROI, plumemskomi, res, '3K')
        
            aod550 = AERONET_dailymean(iyear,imon,iday, 'aod')
            
            aodms = aodm * 1. 
            
            for ifactor in factorval: 
                fine, coarse, bimodal = create_Expcoef(iyear, imon, iday, ifactor)    
   
                
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
                                
                                 
                                t2 = time.time()
                                
                                
                                
                                for expname in [bimodal]:
                                    outputname = rundisamar(expname = expname, year=iyear, month=imon, day=iday, cfgtemplate = 'Config_AAI_354_388_100hPa_SA_Mie.in',\
                                                simulationOnly = 0, Gaussdp = 20, atminterval = DSMpres, aerlev = aerlev,  ps = 1013., SZA = sza[ip], SAA = saa[ip], \
                                                VZA = vza[ip], VAA = vaa[ip], wvl = [354., 388.], AOD = DSMpro, As = list(alb[ip,:]), \
                                                aerosolType = 'Mie', expcoef = expname+'.dat')         
                    
                    
                            
                                    src = '/AAI_sim_Mie_lat_%2.2f_lon_%3.2f.h5' % (lat[ip], lon[ip])
                                    
                                    print src
                                    testname = '/fine_LUT_alh-%1.2f_alt-%1.2f_ni-%1.2f' % (ialh, ialt,ifactor)
                                    
                                    if not os.path.exists(nobackupdir + expname + testname):
                                        os.makedirs(nobackupdir + expname + testname)   
                                    
                                    os.rename(outputdir + expname + '/disamar.h5', outputdir + expname + src)
                     
                                    if os.path.isfile(nobackupdir + expname + testname + src):
                                        os.remove(nobackupdir + expname + testname + src)
                    
                    
                                    shutil.move(outputdir + expname + src, nobackupdir + expname + testname)       
                    
                                                      
                                    t3 = time.time()  
                                    print 'Time period of running DISAMAR per simulation:',t2-t1,'s'
                                    
                    
                                    ######### call AAI module 
                                    disamarhdf = tables.open_file(nobackupdir + expname + testname + src)
                                    
                                    disamarout = disamarhdf.root                    
                                    aairetr = disamarout.parameters.AAI[0]
                                    disamarhdf.close()
                                    
                                    print expname
                                    print '**** AAI retrieved: %1.4f'% aairetr   
        
 
   
t4 = time.time()
print 'Time period of running DISAMAR:',t4-t0,'s'    
    
    
