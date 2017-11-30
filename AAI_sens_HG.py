# -*- coding: utf-8 -*-
"""
Run PYDISAMAR AAI sensitivity HG scheme

@author: sunj
"""


######### load python packages
import os
import shutil
import time
import numpy as np
from DISAMAR import rundisamar


maindir = '/usr/people/sunj/Documents/DISAMAR'
outputdir = maindir + '/pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'

# initialization 
ssaval = np.arange(0.7,0.96,0.05)             # single scattering albedo
ssaval = [0.85]

Angval = np.arange(0,2.1,0.5)               # angstrom exponent
Angval = [1.]

asyval = np.array([0.]+list(np.arange(0.5,1.0,0.1)))                       # asymmetry factor
asyval = [0]

szaval = np.arange(0.,76.,15)
#szaval = [30.]

saaval = np.arange(0.,181.,45)
#saaval = [0.]

vzaval = np.arange(0.,76.,15)
#vzaval = [0.]

vaaval = np.arange(0.,181.,45)
#vaaval = [180.]

aodval = np.arange(0.5,2.1,0.5)
aodval = [1.]

# construct aerosol profile 
L = -0.0065         # lapse rate, unit: K/m
R = 8.31432         # gas constant 
g0 = 9.80665        # gravitational accelaration constant
M = 0.0289644       # air molar mass, unit: kg/mol
T0 = 15.+ 273.15     # surface temperature, unit: K 
P0 = 1013.     # surface pressure, unit: hPa
z0 = 0. 



alhval = np.arange(2.5,8.6,2)
alhval = [4.5]
    
altval = np.arange(0.5,2.1,0.5)
altval = [1.]

dpval = np.array([8,12,20,40,80,100,150,200,300])
dpval = [8]
    
# the following parameters using plume mean value
Ps_avg = 1013.
As_avg = [0.05]*2            # @ 354 and 388 nm 


aaiomi = []
aaidsm = []
aaiAAI = [] 

start=time.time() 

for issa in ssaval: 
    for iAng in Angval:
        for iasy in asyval:
            for isza in szaval:
                for isaa in saaval:
                    for ivza in vzaval: 
                        for ivaa in vaaval:
                            for iaod in aodval:                                 
                                for i, ialh in enumerate(alhval): 
                                    for ialt in altval: 
                                        
                                        z = np.array([40] + [10] + [ialh + ialt / 2.] + [ialh - ialt / 2.])
                                        z = z[::-1]
                                        pres = P0 * (1+L/T0 * (z*1e3-z0)) ** (-g0*M/R/L)     # unit: hPa
                                        
                                        aodpro = np.ones(len(pres))*1e-4
                                        aodpro[1] = iaod / ialt


                                        for idp in dpval:
                                            
#                                            expname = 'AAI_sens_HG_mphy_4' 
                                            expname = 'AAI_sens_HG_geom_4' 
                                            
                                            ######## Test HG: all test
                                            outputname = rundisamar(expname = expname, year=2017, month=1, day=27, cfgtemplate = 'Config_AAI_354_388_100hPa_SA_HG.in',\
                                                        simulationOnly = 0, Gaussdp = idp, atminterval = pres, aerlev = 0,  ps = Ps_avg, SZA = isza , SAA = isaa, \
                                                        VZA = ivza, VAA = ivaa, wvl = [354., 388.], AOD = aodpro, As = list(As_avg), \
                                                        aerosolType = 'HG', ANG = np.ones(len(pres))*iAng, SSA = np.ones(len(pres))*issa, g = np.ones(len(pres))*iasy)
                                    
                                            src = '/AAI_sens_saa-%1.2f_Ang-%1.2f_asy-%1.2f_sza-%03i_saa-%03i_vza-%03i_vaa-%03i_aod-%1.2f_alh-%1.2f_alt-%1.2f_dp-%02i.h5' \
                                                                              %(issa,iAng,iasy,isza,isaa,ivza,ivaa,iaod,ialh,ialt,idp)
                                            print src 
                                                                              
        
                                            if not os.path.exists(nobackupdir + expname):
                                                os.makedirs(nobackupdir + expname)   
                                            
                                            os.rename(outputdir + expname + '/disamar.h5', outputdir + expname + src)
                     
                                            if os.path.isfile(nobackupdir + expname + src):
                                                os.remove(nobackupdir + expname + src)            
                    
                                            shutil.move(outputdir + expname + src, nobackupdir + expname + '/')    
        

        
    
 
end = time.time()    
print 'Time period of running DISAMAR:',end-start,'s'
   
    
    
    
