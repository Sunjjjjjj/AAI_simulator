# -*- coding: utf-8 -*-
"""
Intercomparison among datasets 

AOD: OMI vs MODIS vs AERONET 

@author: sunj
"""


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
from timezone import timezone
from AERONET import AERONET_dailymean

plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)




casename = 'chile201701/'
figdir = '/usr/people/sunj/Dropbox/EGU2017/Figure/'

if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir = '/usr/people/sunj/Dropbox/EGU2017/Figure/'+ casename


years = range(2017,2018)
months = range(10,11)
days = range(17,18)
#days.remove(28)
#days.remove(21)
#days = [27]
aerlat = -33.46
aerlon = -70.66

print 'Calculate time zone'
jetlag = timezone(aerlon)

#ROI = [40.,75., -150.,-45.]      # canada
ROI = [30, 65, -10, 80]      # portugal 
#ROI = [-40,-5,-105,-70]
data = 'L2'


crival = 0
omitot = []
modtot = [] 
aeronet = []
omistat = []
modstat = [] 
omitemp = np.zeros([3])
modtemp = np.zeros([3])

res = 0.5
for iyear in years: 
    for imon in months: 
        for i, iday in enumerate(days):
            
        
            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res, data)
            latomi, lonomi, aodomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AOD', res) 
#            aod550 = AERONET_dailymean(iyear,imon,iday, 'aod')
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
            
            plumemsk = np.logical_or(np.isnan(aodomi0), np.isnan(aodmod0))
            
    #        plumemsk = np.logical_or(plumemsk, aaiomi0 < crival)
#            plumemsk = np.logical_or( aodmod0 < 0, aodmod0 < 0)
    #        plumemsk = np.logical_or(plumemsk, aodomi0 < 0.)
    #        plumemsk = np.logical_or(plumemsk, np.isnan(aaiomi0))
            plumemsk = np.logical_or(plumemsk, aaigome0 < crival)
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
            
            x,y = np.meshgrid(lon,lat)
#            aodomi= griddata((lon, lat), np.array(aodomi0), (x, y), method = 'linear')
            aodomi = np.ma.masked_array(aodomi0,plumemskomi)  
            
    
            lat, lon, aodm, Ang, asy = MODIS4cfg(iyear, imon, iday, ROI, plumemskomi, res, data)
    
#            aodmod= griddata((lon, lat), np.array(aodmod0), (x, y), method = 'linear')
            aodmod = np.ma.masked_array(aodmod0,plumemskomi)  
            
        
            fig1 = plt.figure(figsize=(10.5,6.5), frameon = False)
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
            cbax = fig1.add_axes([0.9125,0.125,0.025,0.775])  
            plt.subplot(gs[0])
            map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], projection='cyl',resolution='c')
            map.drawcoastlines(color='gray',linewidth=1)
#            map.drawcountries(color = 'black',linewidth=0)
#            map.fillcontinents(color='lightgray')
            map.drawmeridians(np.arange(ROI[2]+5,ROI[3]+5,10),color='lightgray',fontsize = 10, labels = [True, False, False, True])
            map.drawparallels(np.arange(ROI[0],ROI[1],10),color='lightgray',fontsize = 10, labels = [True, False, False, True])        
            cb = plt.pcolor(lonomi,latomi,aodomi,cmap = 'rainbow', vmin=0, vmax=3)  
            plt.title('OMAERO %4i-%02i-%02i' %(iyear, imon, iday)) 
            cbar = plt.colorbar(cb,cax = cbax, ticks = np.arange(0,3.1,1))
            cbar.set_label('AOD')
        
            plt.subplot(gs[1])
            map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], projection='cyl',resolution='c')
            map.drawcoastlines(color='gray',linewidth=1)
#            map.drawcountries(color = 'black',linewidth=0)
#            map.fillcontinents(color='lightgray')
            map.drawmeridians(np.arange(ROI[2]+5,ROI[3]+5,10),color='lightgray',labels = [True, False, False, True])
            map.drawparallels(np.arange(ROI[0],ROI[1],10),color='lightgray',labels = [True, False, False, True])       
            cb = plt.pcolor(lonmod,latmod,aodmod,cmap = 'rainbow', vmin=0, vmax=3)  
            plt.title('MODIS %4i-%02i-%02i' %(iyear, imon, iday)) 
#            plt.savefig(figdir + 'AOD_OMI_MODIS_%04i%02i%02i.png' % (iyear, imon, iday), dpi = 300, transparent = False) 
#
#
#            # plot histogram 
#            slope, intercept, r_value, p_value, std_err = stats.linregress(aodomi[aodomi.mask==0],aodmod[aodmod.mask==0])
#            
#            fig2 = plt.figure(figsize=(10.5,4), frameon = False)
#            gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 1]) 
#            cbax = fig2.add_axes([0.915,0.125,0.025,0.775])  
#            plt.subplot(gs[0])
#            sns.kdeplot(aodomi[aodomi.mask==0], color = colorid[0], marker = 's', markevery = 20, label = 'OMAERO max = %1.2f' % (aodomi.max()))
#            sns.kdeplot(aodmod[aodmod.mask==0], color = colorid[1], marker = 's', markevery = 20, label = 'MODIS max = %1.2f' % (aodmod.max()))
#            plt.xlabel('AOD')
#            plt.ylabel('PDF')
#            plt.title('%4i-%02i-%02i' % (iyear, imon, iday))
#            plt.legend(frameon = False, loc = 1, ncol = 1)   
#
#            plt.subplot(gs[1])
#            hb = plt.hexbin(aodomi[aodomi.mask==0], aodomi[aodomi.mask==0], gridsize=50, norm=colors.LogNorm(), cmap='gray_r')
#            plt.plot([0,25],[0,25],'k')
#            plt.plot(np.arange(0,25),np.arange(0,25)*slope+intercept, 'k--')
#            plt.xlim(0,np.max([aodomi[aodomi.mask==0],aodomi[aodomi.mask==0]]))
#            plt.ylim(0,np.max([aodmod[aodmod.mask==0],aodmod[aodmod.mask==0]]))
#            plt.xticks(np.arange(0,np.max([aodomi[aodomi.mask==0],aodomi[aodomi.mask==0]]),2))
#            plt.yticks(np.arange(0,np.max([aodmod[aodmod.mask==0],aodmod[aodmod.mask==0]]),2))
#            cbar = plt.colorbar(hb,cax = cbax)
#            cbar.set_label('Counts')
#            plt.title('%4i-%02i-%02i R = %1.2f' %(iyear, imon, iday, np.corrcoef(aodomi[aodomi.mask==0], aodmod[aodmod.mask==0])[1][0])) 
#            plt.xlabel('OMAERO AOD [-]')
#            plt.ylabel('MODIS AOD [-]')
#            plt.savefig(figdir + 'AOD_OMI_MODIS_stats_%04i%02i%02i.png' % (iyear, imon, iday), dpi = 300, transparent = False) 
            
            
#            # time series 
#            omitemp = [aodomi[aodomi.mask==0].mean(),np.median(aodomi[aodomi.mask==0]),np.sum(aodomi[aodomi.mask==0])] 
#            modtemp = [aodmod[aodmod.mask==0].mean(),np.median(aodmod[aodmod.mask==0]),np.sum(aodmod[aodmod.mask==0])] 
#            omitot.append(aodomi[aodomi.mask==0].data)
#            modtot.append(aodmod[aodmod.mask==0].data)
#            omistat.append(omitemp)
#            modstat.append(modtemp)            
#            aeronet.append(aod550)
#            
#            plt.close('all')
#           
#
#days = np.array(days)
#omistat = np.array(omistat)
#modstat = np.array(modstat)
#aeronet = np.array(aeronet)

#fig3 = plt.figure(figsize=(10.5,4.5))    
#ax1 = fig3.add_axes([0.1,0.225,0.7,0.5])
#plt.plot(days, omistat[:,0], color = colorid[0], marker = 's', linewidth = 1, label = 'OMAERO mean')
#plt.plot(days, omistat[:,1], color = colorid[0], linestyle = '--', marker = 's', linewidth = 1, label = 'OMAERO median')
#plt.plot(days, modstat[:,0], color = colorid[1], marker = 's', linewidth = 1, label = 'MODIS mean')
#plt.plot(days, modstat[:,1], color = colorid[1], linestyle = '--', marker = 's', linewidth = 1, label = 'MYD04 median')
#plt.plot(days, aeronet, color = colorid[2], marker = 's', linewidth = 1, label = 'AERONET mean')
#plt.xlabel('Day ')
#plt.ylabel('AOD [-]')
#plt.xticks(days)
#plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
#plt.savefig(figdir + 'AOD_OMI_MODIS_ts.png', dpi = 300, transparent = False)    
#
#
#
#
#
#fig4 = plt.figure(figsize=(6,3))    
#ax1 = fig4.add_axes([0.1,0.2,0.85,0.6])
#bp1 = plt.boxplot(omitot, positions = days - 0.15, widths=0.3)
#plt.setp(bp1['boxes'], color=colorid[0], linewidth = 2)
#plt.setp(bp1['caps'], color=colorid[0], linewidth = 2)
#plt.setp(bp1['whiskers'], color=colorid[0])
#plt.setp(bp1['medians'], color=colorid[0])
#plt.setp(bp1['means'], color=colorid[0])
#bp2 = plt.boxplot(modtot, positions = days + 0.15, widths=0.3)
#plt.setp(bp2['boxes'], color=colorid[1], linewidth = 2)
#plt.setp(bp2['caps'], color=colorid[1], linewidth = 2)
#plt.setp(bp2['whiskers'], color=colorid[1])
#plt.setp(bp2['medians'], color=colorid[1])
#plt.setp(bp2['means'], color=colorid[1])
#plt.xlabel('Date')
#plt.ylabel(r'$\tau_{aer}$')
#hB, = plt.plot([1,1],color = colorid[0])
#hR, = plt.plot([1,1],color = colorid[1])
#plt.legend((hB, hR),('OMAERO', 'MYD04'))
#hB.set_visible(False)
#hR.set_visible(False)
#plt.ylim([0,2])
#ax1.set_xticks(days)
#ax1.set_xticklabels((days))
#plt.legend(frameon = False, loc = 6, ncol = 2)   
#plt.savefig(figdir + 'AOD_OMI_MODIS_ts_box.png', dpi = 300, transparent = False)      
#
#plt.close('all')