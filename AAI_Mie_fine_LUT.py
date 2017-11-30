"""
Plot AAI simulation:
- spatial pattern 
- statitics of the plume 

@author: sunj
"""


######### load python packages
import sys, os
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec
import glob
import seaborn as sns
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
from sklearn.metrics import mean_squared_error
from scipy import stats
import tables
from scipy.interpolate import griddata 
from OMI import gridOMAERO, OMAERO4cfg
from MODIS import gridMODIS, MODIS4cfg
from GOME2 import gridGOME2
from AERONET import AERONET_dailymean
from scipy import spatial 
from timezone import timezone
from scipy import signal
import subprocess
import shutil
from scipy.cluster.vq import vq, kmeans, whiten
from scipy.stats import iqr


plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)


maindir = '/usr/people/sunj/Documents/DISAMAR'
outputdir = maindir + '/pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'
LUTdir = 'AAI_LUT/'

casename = 'chile201701/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'+ casename



# the following parameters using plume mean value

aaiomi = []
aaidsm = []


start=time.time() 

years = range(2017,2018)
months = range(1,2)
days = range(26,31)
days.remove(28)
#days = [29,30]

aerlat = -33.46
aerlon = -70.66

# calculate time difference
jetlag = timezone(aerlon)


ROIval = [[-40,-5,-105,-70], [-40,-5,-105,-70], [-40,-5,-105,-70], [-40,-5,-105,-70]]
#ROIval = [[-35,-25,-90,-70],[-35,-20,-100,-80],[-30,-15,-100,-85],[-20,-12,-105,-90]] 

crival = 1
res = 0.5

aodest = [] 
alhest = [] 
altest = [] 
factorest = []

LUTfile = 'coarse_LUT_%4i%02i%02i-%4i%02i%02i.h5' % (years[0], months[0], days[0], years[-1], months[-1], days[-1])     
LUT = tables.open_file(nobackupdir + LUTdir + LUTfile, 'r').root
rmse = LUT.Root_mean_suquared_error[:]
alh = LUT.Aerosol_layer_height[:]
alt = LUT.Aerosol_layer_thickness[:]
#aod = LUT.Aerosol_optical_thickness[:]
factor = LUT.Refractive_index_factor[:]

for il in range(len(rmse)): 
    ix, iy, iz = np.where(rmse[il,:,:,:] == np.nanmin(rmse[il,:,:,:]))
#    ix, iy = np.where(rmse[il,:,:] == np.nanmin(rmse[il,:,:]))
    factorest.append(factor[ix[0]])
#    aodest.append(aod[iy[0]])
    
    alhest.append(alh[iy[0]])
    altest.append(alt[iz[0]])
    print rmse[il,ix,iy,iz]
#    print rmse[il,ix,iy]

print 'coarse LUT: ', factorest, alhest, altest 

rms = np.zeros([len(days),5,5,3]) * np.nan 
alh = np.zeros([len(days),5])
alt = np.zeros([len(days),3])
factor = np.zeros([len(days),5])
dsmmean = np.zeros([len(days),5,5,3]) 
dsmmedian = np.zeros([len(days),5,5,3]) 
k = np.zeros([len(days),5,5,3]) 
R = np.zeros([len(days),5,5,3]) 

for iyear in years:
    for imon in months:
        for dd, iday in enumerate(days):
            factorval = np.array([factorest[dd] - 0.2, factorest[dd] - 0.1, factorest[dd], factorest[dd] + 0.1, factorest[dd] + 0.2])
            alhval = np.array([alhest[dd] - 0.5, alhest[dd] - 0.25, alhest[dd], alhest[dd] + 0.25, alhest[dd] + 0.5]) 
            altval = np.array([altest[dd] - 0.25,  altest[dd], altest[dd] + 0.25])
            
            factor[dd,:] = factorval
            alh[dd,:] = alhval 
            alt[dd,:] = altval 

            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROIval[dd], res, '3K')
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROIval[dd], res)
            latomi, lonomi, aaiomi0 = gridOMAERO(iyear, imon, iday, ROIval[dd], jetlag, 'AAI', res) 

            plumemsk = np.logical_or(aodmod0 < 0., aaigome0 < crival) 
            
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROIval[dd], jetlag, plumemsk, crival, res)

            aaiomi = np.ma.masked_array(aaiomi0,plumemskomi)  
            aaigome = np.ma.masked_array(aaigome0,plumemskomi)  
            aodmod = np.ma.masked_array(aodmod0,plumemskomi)  
            
            
            for ff, ifactor in enumerate(factorval):
                for hh, ialh in enumerate(alhval): 
                    for tt, ialt in enumerate(altval): 
                        if ialh > ialt: 
                            testname = '/fine_LUT_alh-%1.2f_alt-%1.2f_ni-%1.2f' % (ialh, ialt, ifactor)
                            print testname
                            
                            for imode in ['bimodal']:            
                                expname = 'AAI_sim_Mie_%s_%4i%02i%02i' %(imode, iyear, imon, iday)
                                
                                filelist = glob.glob(nobackupdir + expname + testname + '/*.h5') # OMI
                                
                                aaidsm = []
                                aaistd = [] 
                                aaiall = [] 
                                aaistdall = [] 
                                aaimean = [] 
                                aaimedian = [] 
                                aoddsm = []
                                for ip in range(len(lat)):
                                    for f in filelist:
                                        if f.find('%2.2f' % (lat[ip])) != -1 and f.find('%3.2f' % (lon[ip])) != -1: 
                                            outf = f
                                    disamarhdf = tables.open_file(outf, 'r')
                                    disamarout = disamarhdf.root
                                    aaidsm.append(disamarout.parameters.AAI[0])
                                    aaistd.append(disamarout.parameters.precisionAAI[0])
                                    disamarhdf.close() 
                                    
                                    aaiall.append(aaidsm)
                                    aaistdall.append(aaistd)
                                    aaimean.append(np.array(aaidsm).mean())
                                    aaimedian.append(np.median(np.array(aaidsm)))
                             
                                aaidsm = np.array(aaidsm)
                                aaiall = np.array(aaiall[0])
                                aaistdall = np.array(aaistdall)
                    
                                aaimean = np.array(aaimean)
                                aaimedian = np.array(aaimedian)
                                
                                # plot DSM ssimulated plume 
                                latnew = np.arange(ROIval[dd][0], ROIval[dd][1], res)
                                lonnew = np.arange(ROIval[dd][2], ROIval[dd][3], res)
                                 
                                x,y = np.meshgrid(lonnew,latnew)
                                aainew = griddata((lon, lat), np.array(aaiall), (x, y), method = 'linear')
    
                                aainew = np.ma.masked_array(aainew,plumemskomi)  
                        
                                # root mean squared error 
                                # AAI IQR 
                                diff = aaidsm - aai 
#                                
#                                 # AOD-AAI IQR 
#                                slope1, intercept1, r_value1, p_value, std_err = stats.linregress(aodmod[aodmod.mask==0],aaiomi[aaiomi.mask==0])
#                                aaiest = aodmod[aodmod.mask==0] * slope1 + intercept1 
#                                diff = aaiest - aai
                                

                                Q1 = np.percentile(diff,25) 
                                Q3 = np.percentile(diff,75)
#                                diff[np.isnan(diff)] = -32676. 
                                
                                valid = np.where((diff >= (Q1 - 1.5*(Q3-Q1))) & (diff <= (Q3 + 1.5*(Q3-Q1))))[0]
#                                rms[dd,ff,hh,tt] = np.sqrt(np.sum((((aaidsm[valid] - aai[valid]))**2)) / np.size(aaidsm[valid]))
                                rms[dd,ff,hh,tt] = np.sqrt(np.sum((((aaidsm - aai))**2)) / np.size(aaidsm))
                                dsmmean[dd,ff,hh,tt] = np.mean(aaidsm)
                                dsmmedian[dd,ff,hh,tt] = np.median(aaidsm)
                                slope, intercept, r_value, p_value, std_err = stats.linregress(aaidsm[valid],aai[valid])
                                k[dd,ff,hh,tt] = slope
                                R[dd,ff,hh,tt] = r_value
                                    
print 'finish the loop'        


# construct LUT
LUTfile = 'fine_LUT_%4i%02i%02i-%4i%02i%02i.h5' % (years[0], months[0], days[0], years[-1], months[-1], days[-1])     
f = tables.open_file('temp.h5', mode = 'w') 
f.create_array('/', 'Root_mean_suquared_error', rms)
f.create_array('/', 'Aerosol_layer_height', alh)
f.create_array('/', 'Aerosol_layer_thickness', alt)
f.create_array('/', 'Refractive_index_factor', factor)
f.create_array('/', 'DSM_AAI_mean', dsmmean)
f.create_array('/', 'DSM_AAI_median', dsmmedian)
f.close()

subprocess.call(["h5repack", "-f", "GZIP=2", "temp.h5", LUTfile])
subprocess.call(["rm", "temp.h5"])
 
# move file to nobackup  
if os.path.isfile(nobackupdir + LUTdir + LUTfile):
    os.remove(nobackupdir + LUTdir + LUTfile)    

shutil.move(LUTfile, nobackupdir + LUTdir + LUTfile)    
            


# read LUT and outputs 
factorest = []
alhest = [] 
altest = [] 
                    
                    
LUTfile = 'fine_LUT_1.5IQR_ROI_%4i%02i%02i-%4i%02i%02i.h5' % (years[0], months[0], days[0], years[-1], months[-1], days[-1])     
LUT = tables.open_file(nobackupdir + LUTdir + LUTfile, 'r').root
rmse = LUT.Root_mean_suquared_error[:]
alh = LUT.Aerosol_layer_height[:]
alt = LUT.Aerosol_layer_thickness[:]
factor = LUT.Refractive_index_factor[:]

for il in range(len(rmse)): 
    ix, iy, iz = np.where(rmse[il,:,:,:] == np.nanmin(rmse[il,:,:,:]))
    factorest.append(factor[il,ix[0]])
    alhest.append(alh[il,iy[0]])
    altest.append(alt[il,iz[0]])

    print 'RMSE: ', rmse[il,ix,iy,iz]          
print factorest, alhest, altest
                    
                    
