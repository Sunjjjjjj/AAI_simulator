# -*- coding: utf-8 -*-
"""
Run PYDISAMAR AAI sensitivity for Mie scheme

@author: sunj
"""


# load python packages ========================================================
import os
import shutil
import time
import numpy as np
from DISAMAR import rundisamar

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


# initialization ==============================================================
# construct aerosol profile 
L = -0.0065         # lapse rate, unit: K/m
R = 8.31432         # gas constant 
g0 = 9.80665        # gravitational accelaration constant
M = 0.0289644       # air molar mass, unit: kg/mol
T0 = 15.+ 273.15     # surface temperature, unit: K 
P0 = 1013.     # surface pressure, unit: hPa
z0 = 0. 


# smoke model 
# size distribution ===========================================================
rgval = np.arange(0.1,0.41,0.05)
rgval = [0.14] 

sigma = 1.5

# refractive index ============================================================
reval = np.arange(1.300, 1.501, 0.05)
#reval = np.arange(1.400, 1.601, 0.05)
reval = [1.5]

imval = np.arange(0.04, 0.11, 0.02)
#imval = np.arange(0.01, 0.09, 0.01)
imval = [0.06]

Re550 = 1.505
Im550 = 0.024

# geometry ====================================================================
szaval = np.arange(0.,76.,15)
szaval = [30.]

saaval = np.arange(0.,181.,45)
saaval = [0.]

vzaval = np.arange(0.,76.,15)
vzaval = [0.]

vaaval = np.arange(0.,181.,45)
vaaval = [180.]


# macro-physics ===============================================================
aodval = np.arange(0.5,2.1,0.5)
#aodval = np.arange(0.2,0.9,0.1)
#aodval = [0.9,1.0]
aodval = [1]

alhval = np.arange(2.5,8.6,2)
#alhval = np.arange(3,5.1,0.25)
#alhval = [5.25,5.5,5.75,6.0]
alhval = [4.5]

altval = np.arange(0.5,2.1,0.5)
#altval = np.arange(0.25,2.1,0.25)
#altval = [0.75]
altval = [1]

# Gaussian division points (sub-layers)========================================
dpval = np.array([8,12,20,40,80,100,150,200,300])
dpval = [20]


# surface condition ===========================================================
Ps_avg = 1013.
As_avg = [0.05]*2            # @ 354 and 388 nm 


aaiomi = []
aaidsm = []
aaiAAI = [] 

# sumulation ==================================================================
start=time.time() 

for irg in rgval: 
    for ire in reval:
        for iim in imval: 
            for isza in szaval:
                for isaa in saaval:
                    for ivza in vzaval: 
                        for ivaa in vaaval:
                            for iaod in aodval:                                 
                                for i, ialh in enumerate(alhval): 
                                    for ialt in altval: 
                                        for idp in dpval:
                                        
#"""
#**** construct aerosol profile
#"""
                                            z = np.array([40] + [10] + [ialh + ialt / 2.] + [ialh - ialt / 2.])
                                            z = z[::-1]
                                            pres = P0 * (1+L/T0 * (z*1e3-z0)) ** (-g0*M/R/L)    # unit: hPa
        
                                            aodpro = np.ones(len(pres))*1e-4
                                            aodpro[1] = iaod / ialt                             # normalize AOD          
                                            
    
#"""
#**** specify expname (output directory) and expfilename (output fire)
#"""  
#                                            expname = 'AAI_sens_Mie_uncertainty_4' 
#                                            expname = 'AAI_sens_Mie_mphy_7'
#                                            expname = 'AAI_sens_Mie_geom_5'
                                            expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (irg, ire, iim)
                                            
                                            
#"""
#**** construct expansion coefficient and scattering coefficient file
#"""
                                        
                                            with open(template,'r') as f: 
                                                content = f.read()
                                                
                                            content = content.replace('XnameX',expfilename)
                                            content = content.replace('XrgX','%1.2f' %(irg))
                                            content = content.replace('XsigmaX','%1.2f' %(sigma))   
                                            content = content.replace('Xre354X','%1.4f' %(ire))
                                            content = content.replace('Xim354X','%1.4f' %(iim))
                                            content = content.replace('Xre388X','%1.4f' %(ire + 0.002))
                                            content = content.replace('Xim388X','%1.4f' %(iim - 0.002))        
                                            content = content.replace('Xre550X','%1.4f' %(ire + 0.01))
                                            content = content.replace('Xim550X','%1.4f' %(iim - 0.01))  
                                            
                                            with open(expindir + expfilename+ '.in','w') as f:
                                                f.write(content)                        
                                            
                                            if os.path.isfile(exedir + 'Config.in'):
                                                os.remove(exedir + 'Config.in')
                                            os.symlink(expindir + expfilename + '.in', exedir + 'Config.in')
                                            os.system(exedir + "createExpCoefFile.exe") 
                                            
                                            shutil.copy2(expoutdir + expfilename+'.dat', expcoefdir + expfilename+'.dat')
                                            
                                            
#"""
#**** run DISAMAR
#"""
                                            outputname = rundisamar(expname = expname, year = 2017, month=1, day=27, cfgtemplate = 'Config_AAI_354_388_100hPa_SA_Mie.in',\
                                                        simulationOnly = 0, Gaussdp = idp, atminterval = pres, aerlev = 0,  ps = Ps_avg, SZA = isza, SAA = isaa, \
                                                        VZA = ivza, VAA = ivaa, wvl = [354., 388.], AOD = aodpro, As = list(As_avg), \
                                                        aerosolType = 'Mie', expcoef = expfilename+'.dat')         
                    
                                    
                                            src = '/AAI_sens_rg-%1.2f_re-%1.4f_im-%1.4f_sza-%03i_saa-%03i_vza-%03i_vaa-%03i_aod-%1.2f_alh-%1.2f_alt-%1.2f_dp-%02i.h5' \
                                                                                %(irg,ire,iim,isza,isaa,ivza,ivaa,iaod,ialh,ialt,idp)
                                            
                                            print src
                                            if not os.path.exists(nobackupdir + expname):
                                                os.makedirs(nobackupdir + expname)   
                                            
                                            os.rename(outputdir + expname + '/disamar.h5', outputdir + expname + src)
                     
                                            if os.path.isfile(nobackupdir + expname + src):
                                                os.remove(nobackupdir + expname + src)            
                    
                                            shutil.move(outputdir + expname + src, nobackupdir + expname + '/')    
                                              
end = time.time()    
print 'Time period of running DISAMAR:',end-start,'s'
                    

    
 
   
    
    
    
