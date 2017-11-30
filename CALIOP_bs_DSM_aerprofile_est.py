# -*- coding: utf-8 -*-
"""
Plot CALIOP footprint


@author: sunj
"""


######### load python packages
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap
from scipy import stats
from OMI import rawOMAERO, gridOMAERO, OMAERO4cfg
from MODIS import gridMODIS, MODIS4cfg
from GOME2 import gridGOME2
from CALIPSO import trackCALIOP
from timezone import timezone



plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':12,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)

#casename = 'portugal201710/'
casename = 'chile201701/'
#casename = 'canada201708/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'


if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'+ casename


years = range(2017,2018)
months = range(1,2)
days = range(26,31)
days.remove(28)
#days=[26]

aerlat = -33.46
aerlon = -70.66
plat = -34.39
plon = -72.0
clat = -35.33 
clon = -72.42 


print 'Calculate time zone'
jetlag = timezone(aerlon)

#ROI = [45, 65, 0, 50]      # portugal 
ROI = [-40,-5,-105,-70] # chile
#ROI = [45.,80., -150.,-60.]      # canada

#lat0 = -22.0
#lon0 = -70.5
#nlay = 4 

crival = 2
omitot = []
gometot = [] 
omistat = []
gomestat = [] 
omitemp = np.zeros([3])
gometemp = np.zeros([3])

ROIval = [[-40,-30,-80,-70],[-35,-25,-90,-80],[-35,-25,-85,-75],[-25,-15,-90,-80]]

#trackval = [3] 
trackval = [2,0,2,2] 

# AAI IQR 
alhval = [5,4,5.25,2.5]
altval = [0.5,1,1,0.25]

## AOD-AAI IQR 
#alhval = [5.5,4.25,5.5,2.5]
#altval = [1.25,1.25,0.5,0.25]

res = 0.5
for iyear in years:
    for imon in months: 
        for i, iday in enumerate(days):
#            i = 3
#            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res,'3K')
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
            latomi, lonomi, aaiomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AAI',  res) 
            latomi, lonomi, aodomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AOD',  res) 
            latcal, loncal, bs532, z = trackCALIOP(iyear, imon, iday, ROI,'L1')


            plumemsk = np.logical_or(aaigome0 < crival, aaigome0 < crival)
            
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
                        
            aaiomi = np.ma.masked_array(aaiomi0, plumemskomi)
            aodomi = np.ma.masked_array(aodomi0, plumemskomi)
            aaigome = np.ma.masked_array(aaigome0, plumemskomi)
#            aodmod = np.ma.masked_array(aodmod0, plumemskomi)
            

            itrack = trackval[i]

            latm,zm = np.meshgrid(latcal[itrack],z[itrack]) 
            lonm,zm = np.meshgrid(loncal[itrack],z[itrack]) 
            
            msk = np.logical_or(bs532[itrack] == 0., bs532[itrack] < 0)
            msk = np.logical_or(msk, latm < ROIval[i][0])
            msk = np.logical_or(msk, latm > ROIval[i][1])
            msk = np.logical_or(msk, lonm < ROIval[i][2])
            msk = np.logical_or(msk, lonm > ROIval[i][3])
            bs532m = np.ma.masked_array(bs532[itrack],msk) 
            latm = np.ma.masked_array(latm,msk) 
            lonm = np.ma.masked_array(lonm,msk) 
            
            fig1 = plt.figure(figsize=(8.5,4.5))
            cbax = fig1.add_axes([0.825,0.175,0.03,0.675]) 
            ax1 = fig1.add_axes([0.1,0.175,0.7,0.675]) 
            ax2 = ax1.twiny()
            cb = ax1.pcolor(latm,zm,bs532m*1e3,cmap='CMRmap_r',norm=colors.LogNorm(vmin=1, vmax=100))
            ax2.pcolor(lonm,zm,bs532m*1e3,cmap='rainbow',norm=colors.LogNorm(vmin=1, vmax=100))
            cbar = plt.colorbar(cb,cax = cbax, ticks = np.arange(1,100.1))
            cbar.set_label(r'$\beta$ [ 10$^{-3}$km$^{-1}$sr$^{-1}$]')
            ax2.cla()
            ax1.set_xlabel('Latitude [$\circ$]')
            ax1.set_xticks(np.arange(int(latm.min()),int(latm.max())+0.1,2))
            ax2.set_xlabel('Longitude [$\circ$]')
            ax2.set_xticks(np.arange(int(lonm.min()),int(lonm.max())+0.1,1))
            ax1.set_xlim(ROIval[i][0],ROIval[i][1])
            plt.ylim(0,10)
            ax1.set_ylabel('Altitude [km]')
            ax1.tick_params(which='major', length=7, direction = 'in')
            ax2.tick_params(which='major', length=7, direction = 'in')
            plt.axhline(alhval[i],color = 'k',linewidth = 3)
            plt.axhline(y=alhval[i] - altval[i]/2., color = 'k', linestyle = '--',linewidth = 1.5)
            plt.axhline(y=alhval[i] + altval[i]/2., color = 'k', linestyle = '--',linewidth = 1.5)
            plt.title('%4i-%02i-%02i' %(iyear, imon, iday)) 
            
            plt.savefig(figdir + 'ALH_OMI_CALIOP_bs532_L1_%02i_%04i%02i%02i.png' % (itrack, iyear, imon, iday), dpi = 300, transparent = True) 
    
    
            plt.close('all')
           

