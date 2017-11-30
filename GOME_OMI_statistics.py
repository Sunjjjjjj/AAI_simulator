# -*- coding: utf-8 -*-
"""
Intercomparison among datasets 

- AAI: OMI vs GOME2

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
from AERONET import AERONET_day
from timezone import timezone
from CALIPSO import trackCALIOP 

plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)


#casename = 'chile201701/'
#casename = 'canada201708/'
casename = 'portugal201710/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'

if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir += casename




years = range(2017,2018)
months = range(10,11)
days = range(17,18)

#months = range(1,2)
#days = range(26,31)
#days.remove(28)

aerlat = -33.46
aerlon = -70.66


print 'Calculate time zone'
jetlag = timezone(aerlon)

#ROI = [-40,-5,-105,-70] # chile
ROI = [45, 65, 0, 50]      # portugal 
ROIaer = [50, 65, 5, 45]
#ROIaer = [50, 70, -150, -45] # canada

data ='L2'

crival = 2
omitot = []
gometot = [] 
omistat = []
gomestat = [] 
omitemp = np.zeros([3])
gometemp = np.zeros([3])

res = 0.5
for iyear in years:
    for imon in months: 
        for i, iday in enumerate(days):
#            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res, data)
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
            latomi, lonomi, aaiomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AAI',  res) 
            latomi, lonomi, aodomi0 = gridOMAERO(iyear, imon, iday, ROI, jetlag, 'AOD',  res) 
            latcal, loncal, bs532, z = trackCALIOP(iyear, imon, iday, ROI,'L2')
            
#            aerosite = AERONET_day(iyear, imon, iday, ROIaer)
##            plumemsk = np.zeros(aaiomi0.shape)
            plumemsk = np.logical_or(np.isnan(aaiomi0), np.isnan(aaigome0))
            plumemsk = np.logical_or(plumemsk, aaigome0 < crival)
            plumemsk = np.logical_or(plumemsk, aaiomi0 < crival)
#            
##            plumemsk = np.logical_or(aaigome0 < crival, aodmod0 < 0.)
##            plumemsk = np.logical_or(plumemsk, aodomi0 < 0.)
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
#            print aai.shape
            aaiomi = np.ma.masked_array(aaiomi0, plumemskomi)
            aaigome = np.ma.masked_array(aaigome0, plumemskomi)
#            aodmod = np.ma.masked_array(aodmod0, plumemskomi)
            
#            loc = np.column_stack([lon,lat])
#            kmeans = KMeans(n_clusters = 8)
#            kmeans.fit(loc)
#            loc_km = kmeans.predict(loc)
#            plt.figure()
#            plt.scatter(loc[:, 0], loc[:, 1], c=loc_km, s=50, cmap='viridis')
            
            
            for i in range(2,aerosite.shape[1]):
                plt.figure()
                plt.plot(aerosite[:,i],'s-')
                plt.plot(np.ones(aerosite.shape[0])*aerosite[:,i].mean())


               
            # plot spatical pattern 
            fig1 = plt.figure(figsize=(10.5,6.5), frameon = False)
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
            cbax = fig1.add_axes([0.9125,0.125,0.025,0.775])  
            plt.subplot(gs[0])
            map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], projection='cyl',resolution='c')
            map.drawcoastlines(color='gray',linewidth=1)
#            map.drawcountries(color = 'black',linewidth=1)
#            map.fillcontinents(color='lightgray')
            map.drawmeridians(np.arange(ROI[2]+5,ROI[3]+5,5), color='lightgray',fontsize = 12, labels = [True, False, False, True])
            map.drawparallels(np.arange(ROI[0],ROI[1],5),color='lightgray', fontsize = 12, labels = [True, False, False, True])        
            cb = plt.pcolor(lonomi, latomi, aaiomi, cmap = 'rainbow', vmin=0, vmax=8)          
            plt.title('OMAERO %4i-%02i-%02i' %(iyear, imon, iday)) 
            cbar = plt.colorbar(cb,cax = cbax, ticks = np.arange(0,8.1,2))
            cbar.set_label('AAI')
            for itrack in range(len(latcal)):
                plt.plot(loncal[itrack],latcal[itrack],color = colorid[itrack], linewidth =4)
#            plt.scatter(aerosite[:,1],aerosite[:,0],marker = 's', color = 'k')
        
            plt.subplot(gs[1])
            map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], lat_0 = 0, lon_0 = 0, projection='cyl',resolution='c')
            map.drawcoastlines(color='gray',linewidth=1)
#            map.drawcountries(color = 'black',linewidth=0)
#            map.fillcontinents(color='lightgray')    
            map.drawmeridians(np.arange(ROI[2]+5,ROI[3]+5,5),color='lightgray',fontsize = 12, labels = [True, False, False, True])
            map.drawparallels(np.arange(ROI[0],ROI[1],5),color='lightgray', fontsize = 12, labels = [True, False, False, True])  
            cb = plt.pcolor(longome,latgome,aaigome, cmap = 'rainbow', vmin=0, vmax=8)  
            plt.title('GOME-2 %4i-%02i-%02i' %(iyear, imon, iday))
            for itrack in range(len(latcal)):
                plt.plot(loncal[itrack],latcal[itrack],color = colorid[itrack], linewidth =4)
#            plt.scatter(aerosite[:,1],aerosite[:,0],marker = 's', color = 'k')

            plt.savefig(figdir + 'AAI_OMI_GOME2_%04i%02i%02i.png' % (iyear, imon, iday), dpi = 300, transparent = False) 
        
            
#            # plot histogram 
#            slope, intercept, r_value, p_value, std_err = stats.linregress(aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0])
#            
#            fig2 = plt.figure(figsize=(10.5,4), frameon = False)
#            gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 1]) 
#            cbax = fig2.add_axes([0.915,0.125,0.025,0.775])  
#            plt.subplot(gs[0])
#            sns.kdeplot(aaiomi[aaiomi.mask==0], color = colorid[0], marker = 's', markevery = 20,label = 'OMAERO max = %1.2f' % (aaiomi.max()))
#            sns.kdeplot(aaigome[aaigome.mask==0], color = colorid[1], marker = 's', markevery = 20,label = 'GOME2 max = %1.2f' % (aaigome.max()))
#            plt.xlabel('AAI [-]')
#            plt.ylabel('PDF')
#            plt.xticks(np.arange(0,12.1,2))
#            plt.title('%4i-%02i-%02i' % (iyear, imon, iday))
#            plt.legend(frameon = False, loc = 1, ncol = 1)   
#
#            plt.subplot(gs[1])
#            hb = plt.hexbin(aaiomi[aaiomi.mask==0], aaigome[aaigome.mask==0], gridsize=50, norm=colors.LogNorm(), cmap='gray_r')
#            plt.plot([0,25],[0,25],'k')
#            plt.plot(np.arange(0,25),np.arange(0,25)*slope+intercept, 'k--')
#            plt.xlim(0,np.max([aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0]]))
#            plt.ylim(0,np.max([aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0]]))
#            plt.xticks(np.arange(0,np.max([aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0]]),2))
#            plt.yticks(np.arange(0,np.max([aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0]]),2))
#            cbar = plt.colorbar(hb,cax = cbax)
#            cbar.set_label('Counts')
#            plt.title('%4i-%02i-%02i R = %1.2f' %(iyear, imon, iday, np.corrcoef(aaiomi[aaiomi.mask==0], aaigome[aaigome.mask==0])[1][0])) 
#            plt.xlabel('OMAERO AAI [-]')
#            plt.ylabel('GOME-2 AAI [-]')
#            plt.savefig(figdir + 'AAI_OMI_GOME2_stats_%04i%02i%02i.png' % (iyear, imon, iday), dpi = 300, transparent = False) 
#        
#            # time series 
#
#            omitemp = [aaiomi[aaiomi.mask==0].mean(),np.median(aaiomi[aaiomi.mask==0]),np.sum(aaiomi[aaiomi.mask==0])] 
#            gometemp = [aaigome[aaigome.mask==0].mean(),np.median(aaigome[aaigome.mask==0]),np.sum(aaigome[aaigome.mask==0])] 
#            omistat.append(omitemp)
#            gomestat.append(gometemp)
#            omitot.append(aaiomi[aaiomi.mask==0].data)
#            gometot.append(aaigome[aaigome.mask==0].data)      
            
#            plt.close('all')
#           
#
#days = np.array(days)
#omistat = np.array(omistat)
#gomestat = np.array(gomestat)
#
##fig3 = plt.figure(figsize=(10.5,4.5))    
##ax1 = fig3.add_axes([0.1,0.225,0.7,0.5])
##plt.plot(days, omistat[:,0], color = colorid[0], marker = 's', linewidth = 1, label = 'OMAERO mean')
##plt.plot(days, omistat[:,1], color = colorid[0], linestyle = '--', marker = 's', linewidth = 1, label = 'OMAERO median')
##plt.plot(days, gomestat[:,0], color = colorid[1], marker = 's', linewidth = 1, label = 'GOME mean')
##plt.plot(days, gomestat[:,1], color = colorid[1], linestyle = '--', marker = 's', linewidth = 1, label = 'GOME median')
##plt.xlabel('Day ')
##plt.ylabel('AAI [-]')
##plt.xticks(days)
##plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
##plt.savefig(figdir + 'AAI_OMI_GOME2_ts.png', dpi = 300, transparent = False)     
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
#bp2 = plt.boxplot(gometot, positions = days + 0.15, widths=0.3)
#plt.setp(bp2['boxes'], color=colorid[1], linewidth = 2)
#plt.setp(bp2['caps'], color=colorid[1], linewidth = 2)
#plt.setp(bp2['whiskers'], color=colorid[1])
#plt.setp(bp2['medians'], color=colorid[1])
#plt.setp(bp2['means'], color=colorid[1])
#plt.xlabel('Date')
#plt.ylabel('AAI')
#hB, = plt.plot([1,1],color = colorid[0])
#hR, = plt.plot([1,1],color = colorid[1])
#plt.legend((hB, hR),('OMAERO', 'GOME2'))
#hB.set_visible(False)
#hR.set_visible(False)
#plt.ylim([0,10])
#ax1.set_xticks(days)
#ax1.set_xticklabels((days))
#plt.legend(frameon = False, loc = 6, ncol = 2)   
##plt.savefig(figdir + 'AAI_OMI_GOME2_ts_box.png', dpi = 300, transparent = False)     
#
#plt.close('all')
            
            
