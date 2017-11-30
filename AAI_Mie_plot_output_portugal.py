"""
Plot AAI simulation:
- spatial pattern 
- statitics of the plume 

@author: sunj
"""


######### load python packages
import os
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
from CALIPSO import trackCALIOP 
from scipy import spatial 
from timezone import timezone
from scipy import signal



plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)


maindir = '/usr/people/sunj/Documents/DISAMAR'
outputdir = maindir + '/pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'
LUTdir = 'AAI_LUT/'

casename = 'portugal201710/'
#casename = 'chile201701/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'+ casename



# the following parameters using plume mean value

aaiomi = []
aaidsm = []


start=time.time() 

years = range(2017,2018)
months = range(10,11)
days = range(17,18)



# calculate time difference
jetlag = -5


ROI = [45, 65, 0, 50]      # portugal 
ROIaer = [50, 65, 5, 45]

crival = 2
res = 0.5


# AAI IQR 
factorval = [3]
aodval = [2.5] 
# AOD-AAI IQR 
#alhval = [5.5,4.25,5.5,2.5]
#altval = [1.25,1.25,0.5,0.25]
#factorval = [2.9,3.2,1.3,2.3]

# CALIOP track
trackval = [1,1,2,0] 


# initialization 
#rms = np.zeros([len(days),len(factorval),len(alhval), len(altval)]) 
#dsmmean = np.zeros([len(days),len(factorval),len(alhval), len(altval)]) 
#dsmmedian = np.zeros([len(days),len(factorval),len(alhval), len(altval)]) 
#k = np.zeros([len(days),len(factorval), len(alhval), len(altval)]) 

AAIomi = []
AAIdsm = [] 
AOD = [] 

AAIomiall = []
AAIdsmall = [] 
AODall = [] 
 

for iyear in years:
    for imon in months:
        for dd, iday in enumerate(days):

#            latmod, lonmod, aodmod0 = gridMODIS(iyear, imon, iday, ROI, res, '3K')
            latgome, longome, aaigome0 = gridGOME2(iyear, imon, iday, ROI, res)
            latomi, lonomi, aaiomi0 = gridOMAERO(iyear, imon, iday, ROI, -5, 'AAI', res) 
            latomi, lonomi, aodomi0 = gridOMAERO(iyear, imon, iday, ROI, -5, 'AOD', res) 
            latcal, loncal, bs532, z = trackCALIOP(iyear, imon, iday, ROI)
            
            plumemsk = np.logical_or(aaiomi0 < 0., aaigome0 < crival) 
#            plumemsk = np.logical_or(aodmod0 < 0., aodmod0 < 0.) 
            
            lat, lon, aai, aod, aodstd, sza, saa, vza, vaa, ssa, Ps, plumemskomi, alb = OMAERO4cfg(iyear, imon, iday, ROI, jetlag, plumemsk, crival, res)
            
            aaiomi = np.ma.masked_array(aaiomi0,plumemskomi)  
            aaigome = np.ma.masked_array(aaigome0,plumemskomi)  
            aodomi = np.ma.masked_array(aodomi0,plumemskomi)    
#            aodmod = np.ma.masked_array(aodmod0,plumemskomi)
            
            iaod = aodval[dd]
            ifactor = factorval[dd]

#            testname = '/fine_LUT_alh-%1.2f_alt-%1.2f_ni-%1.2f' % (ialh, ialt, ifactor)
            testname = '/aod-%1.2f_ni-%1.2f' % (iaod, ifactor)
            print testname
            
            for imode in ['bimodal']:            
                expname = 'AAI_sim_Mie_%s_%4i%02i%02i' %(imode, iyear, imon, iday)
    
                
                filelist = glob.glob(nobackupdir + expname + testname + '/*.h5') # OMI
                
                aaidsm = []
                aaistd = [] 
#                aaiall = [] 
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
                    
                    aaistdall.append(aaistd)
                    aaimean.append(np.array(aaidsm).mean())
                    aaimedian.append(np.median(np.array(aaidsm)))
             
                aaidsm = np.array(aaidsm)
                aaistdall = np.array(aaistdall)
    
                aaimean = np.array(aaimean)
                aaimedian = np.array(aaimedian)
                
                # plot DSM ssimulated plume 
                latnew = np.arange(ROI[0], ROI[1], res)
                lonnew = np.arange(ROI[2], ROI[3], res)
                 
                x,y = np.meshgrid(lonnew,latnew)
                aainew = griddata((lon, lat), np.array(aaidsm), (x, y), method = 'linear')
                aainew = np.ma.masked_array(aainew,plumemskomi)  
                
                mu = np.cos(sza / 180. * np.pi)
                mu0 = np.cos(vza / 180. * np.pi)
                raa = (saa + 180 - vaa) / 180 * np.pi
                theta1 = np.arccos(-mu*mu0 + np.sqrt(1-mu**2)*np.sqrt(1-mu0**2)*np.cos(raa)) / np.pi*180                           
   
                albnew = griddata((lon, lat), np.array(vaa), (x, y), method = 'linear')
                albnew = np.ma.masked_array(albnew,plumemskomi)  
#                #### quality control 
#                diff = aainew - aaiomi 
#                
#                slope1, intercept1, r_value1, p_value, std_err = stats.linregress(aodmod[aodmod.mask==0],aaiomi[aaiomi.mask==0])
#                aaiest = aodmod * slope1 + intercept1 
##                print slope1, intercept1, r_value1                
#                diff = aaiest - aaiomi
#                
#                Q1 = np.percentile(diff[diff.mask==0],25) 
#                Q3 = np.percentile(diff[diff.mask==0],75)
#                diff[np.isnan(diff)] = -32676. 
#                
#                invalid = np.logical_or(plumemsk, diff < (Q1 - 1.5 * (Q3-Q1))) 
#                invalid = np.logical_or(invalid, diff > (Q3 + 1.5 * (Q3-Q1))) 
#                
#
#                
#                aainewm = np.ma.masked_array(aainew,invalid)
#                aaiomim = np.ma.masked_array(aaiomi,invalid)
        
                # plot spatial pattern 
                fig1 = plt.figure(figsize=(10.5,6.5), frameon = False)
                cbax = fig1.add_axes([0.9125,0.125,0.025,0.775])  
                gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1]) 
                
                plt.subplot(gs[0])
                map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], projection='cyl',resolution='c')
                map.drawcoastlines(color='black',linewidth=1)
                map.drawmeridians(np.arange(ROI[2]+5,ROI[3]+5,10),color='lightgray', fontsize = 12, labels = [True, False, False, True])
                map.drawparallels(np.arange(ROI[0],ROI[1],10),color='lightgray', fontsize = 12, labels = [True, False, False, True])        
                cb = plt.pcolor(x,y,aaiomi, cmap = 'rainbow', vmin=0, vmax=8)    
                
#                for itrack in range(len(latcal)):
#                    plt.plot(loncal[itrack],latcal[itrack],color = 'lightgray', linewidth =4)
#                
#                itrack = trackval[dd]
#                idxROI = np.where((latcal[itrack]>ROIval[dd][0]) & (latcal[itrack]<ROIval[dd][1]) & (loncal[itrack]>ROIval[dd][2]) & (loncal[itrack]<ROIval[dd][3]))[0]
#    
#                plt.plot(loncal[itrack][idxROI],latcal[itrack][idxROI], marker = 'D', color = 'k', markersize = 4, linewidth = 0, markevery = 25)
#                plt.scatter(aerlon,aerlat,marker = '^',color = 'k',s = 80)
#                plt.scatter(plon,plat,marker = '^',color = 'r',s = 80)
#                plt.scatter(clon,clat,marker = '^',color = 'r',s = 80)
                plt.tick_params(which='major', length=7, direction = 'in')
                plt.title('OMAERO %4i-%02i-%02i' %(iyear,imon,iday))
                cbar = plt.colorbar(cb,cax = cbax, ticks = np.arange(0,8.1,2))
                cbar.set_label('AAI')


                plt.subplot(gs[1])
                map = Basemap(llcrnrlon=ROI[2],llcrnrlat=ROI[0],urcrnrlon=ROI[3],urcrnrlat=ROI[1], projection='cyl',resolution='c')
                map.drawcoastlines(color='black',linewidth=1)
                map.drawmeridians(np.arange(ROI[2]+5,ROI[3]+5,10),color='lightgray',labels = [True, False, False, True])
                map.drawparallels(np.arange(ROI[0],ROI[1],10),color='lightgray',labels = [True, False, False, True])        
#                for itrack in range(len(latcal)):
#                    plt.plot(loncal[itrack],latcal[itrack],color = 'lightgray', linewidth =4)
#                
#                for itrack in range(len(latcal)):
#                    plt.plot(loncal[itrack],latcal[itrack],color = 'lightgray', linewidth =4)
#                
#                itrack = trackval[dd]
#                idxROI = np.where((latcal[itrack]>ROI[0]) & (latcal[itrack]<ROI[1]) & (loncal[itrack]>ROI[2]) & (loncal[itrack]<ROI[3]))[0]
#    
#                idxROI = np.where((latcal[itrack]>ROIval[dd][0]) & (latcal[itrack]<ROIval[dd][1]) & (loncal[itrack]>ROIval[dd][2]) & (loncal[itrack]<ROIval[dd][3]))[0]
#    
#                plt.plot(loncal[itrack][idxROI],latcal[itrack][idxROI], marker = 'D', color = 'k', markersize = 4, linewidth = 0, markevery = 25)

                cb = plt.pcolor(x,y,aainew, cmap = 'rainbow', vmin=0, vmax=8)  
#                plt.scatter(aerlon,aerlat,marker = '^',color = 'k',s = 80)
#                plt.scatter(plon,plat,marker = '^',color = 'r',s = 80)
#                plt.scatter(clon,clat,marker = '^',color = 'r',s = 80)
                plt.title('DISAMAR %4i-%02i-%02i' %(iyear,imon,iday))
                plt.savefig(figdir + 'AAI_sim_Mie_%s_%04i%02i%02i_alh-%1.2f_alt-%1.2f_ni-%1.2f.png' % (imode,iyear, imon, iday, ialh, ialt, ifactor), dpi = 300, transparent = False) 

                
                
#                 plot statistics
                diff = aaidsm - aai
                
#                slope1, intercept1, r_value1, p_value, std_err = stats.linregress(aodmod[aodmod.mask==0].data,aai)
#                aaiest = aodmod[aodmod.mask==0].data * slope1 + intercept1 
#                diff = aaiest - aai
                
                Q1 = np.percentile(diff,25) 
                Q3 = np.percentile(diff,75)

                valid = np.where((diff >= (Q1 - 1.5*(Q3-Q1))) & (diff <= (Q3 + 1.5*(Q3-Q1))) & \
                        (lat>ROI[0]) & (lat<ROI[1]) & (lon>ROI[2]) & (lon<ROI[3]))[0]
                
                aaidsmv = aaidsm[valid]
                aaiomiv = aai[valid]   
#                aodmodv = aodmod[aodmod.mask==0].data[valid]
#                aodomiv = aodomi[aodmod.mask==0].data[valid]
#                aaigomev = aaigome[aodmod.mask==0].data[valid] 
                
                print 'AAI median difference:       DSM = %1.4f, OMI = %1.4f, diff = %1.2f' %(np.median(aaidsm), np.median(aai), (np.median(aaidsm) - np.median(aai)) / np.median(aai) *1e2)
#                print 'AAI median difference valid: DSM = %1.4f, OMI = %1.4f, diff = %1.2f' %(np.median(aaidsmv),np.median(aaiomiv), (np.median(aaidsmv) - np.median(aaiomiv))/np.median(aaiomiv)*1e2)
#                print 'AOD median difference:', np.median(aodmod[aodmod.mask==0]), np.median(aodomi[aodomi.mask==0]), (np.median(aodmod[aodmod.mask==0]) - np.median(aodomi[aodomi.mask==0])) 
#                print 'AOD median difference valid: diff = %1.4f' % (np.median(aodmodv) -  np.median(aodomiv))

#                print 'AAI difference median:', np.median(aaidsm - aai)
#                print 'AAI difference median valid:', np.median(aaidsmv - aaiomiv) 
#                print 'AOD difference median:', np.median(aodmod[aodmod.mask==0] - aodomi[aodomi.mask==0])
#                print 'AOD difference median valid: diff = %1.4f' % (np.median(aodmodv - aodomiv))


                slope1, intercept1, r_value1, p_value, std_err = stats.linregress(aaiomi[aaiomi.mask==0],aaidsm)
#                print 'R:', r_value1

                
#                fig2 = plt.figure(figsize=(13.5,4.5), frameon = False)
#                gs = gridspec.GridSpec(1,3, width_ratios=[0.85,0.001,1.25])
#                plt.subplot(gs[2])
#                plt.hist(aodmodv - aodomiv, bins = 10, linewidth=1 ,color = 'k', histtype='step', normed = True, label = 'MODIS - OMAERO')
#                plt.hist(aaidsmv - aaiomiv, bins = 10, linewidth=2 ,color = colorid[0], histtype='step', normed = True, label = 'DISAMAR - OMAERO')
#                plt.hist(aaigomev - aaiomiv, bins = 10, linewidth=2, histtype='step',  normed = True, color = colorid[5], label = 'GOME2 - OMAERO')
#                                
#                plt.title('%4i-%02i-%02i' % (iyear, imon, iday))
#                plt.xlabel(r'AAI / $\tau_{aer}$ difference')
#                plt.ylabel('Probability')
#                plt.xlim(-3,3)
#                plt.ylim(0.6)
#                plt.yticks(np.arange(0,6.1,1))
#                plt.legend(frameon=False)
#                
#                plt.subplot(gs[0])
#                slope1, intercept1, r_value1, p_value, std_err = stats.linregress(aaiomiv,aaidsmv)
#                slope2, intercept2, r_value2, p_value, std_err = stats.linregress(aaiomi[aaiomi.mask==0],aaidsm)
#                plt.scatter(aai,aaidsm, s = 12, marker = 's',facecolor='none',linewidth=1, color = colorid[5])
#                plt.scatter(aaiomiv,aaidsmv, s = 12, facecolor='none',linewidth=1, color = colorid[0])
#                plt.plot([0,25],[0,25],'k')
#                plt.plot(np.arange(0,25),np.arange(0,25)*slope2+intercept2, linestyle = '--', linewidth =2, color = colorid[5], label = 'k = %1.2f b = %1.2f R = %1.2f' % (slope2, intercept2, r_value2))                
#                plt.plot(np.arange(0,25),np.arange(0,25)*slope1+intercept1, linestyle = '--', linewidth =2, color = colorid[0], label = 'Data in IQR \n k = %1.2f b = %1.2f R = %1.2f' % (slope1, intercept1, r_value1))
#                plt.xlim(0, np.max([aaidsm,aaiomi[aaiomi.mask==0]]))
#                plt.ylim(0, np.max([aaidsm,aaiomi[aaiomi.mask==0]]))
#                plt.ylabel('DISAMAR AAI')
#                plt.xlabel('OMAERO AAI')
#                plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(0,0.85))   
#                plt.savefig(figdir + 'AAI_sim_Mie_stats_MOD_%s_%04i%02i%02i_alh-%1.2f_alt-%1.2f_ni-%1.2f.png' % (imode, iyear, imon, iday, ialh, ialt, ifactor), dpi = 300, transparent = False) 

                
#            plt.subplot(gs[2])
#            slope1, intercept1, r_value1, p_value, std_err = stats.linregress(aaiomi[aaiomi.mask==0][valid],aaigome[aaigome.mask==0][valid])
#            slope2, intercept2, r_value2, p_value, std_err = stats.linregress(aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0])
#            plt.scatter(aai, aaigome[aaigome.mask==0], facecolor='none',linewidth=1, marker='s',s = 14, color = colorid[5])
#            plt.scatter(aai[valid], aaigome[aaigome.mask==0][valid], facecolor='none',linewidth=1, marker='s',s = 14, color = colorid[0])
#            plt.plot([0,25],[0,25],'k')
#            plt.plot(np.arange(0,25),np.arange(0,25)*slope2+intercept2, linestyle = '--',linewidth=2, color = colorid[5], label = 'k = %1.2f b = %1.2f R = %1.2f' % (slope2, intercept2, r_value2))
#            plt.plot(np.arange(0,25),np.arange(0,25)*slope1+intercept1, linestyle = '--',linewidth=2, color = colorid[0], label = 'Data in IQR \n k = %1.2f b = %1.2f R = %1.2f' % (slope1, intercept1, r_value1))
#            plt.xlim(0, np.max([aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0]]))
#            plt.ylim(0, np.max([aaiomi[aaiomi.mask==0],aaigome[aaigome.mask==0]]))
#            plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(0,0.85))   
#            plt.xlabel('OMAERO AAI')
#            plt.ylabel('GOME2 AAI')
#                cbar = plt.colorbar(hb,cax = cbax)
#                cbar.set_label('Counts')


                
                AAIomi = AAIomi + list(aai[valid])
                AAIdsm = AAIdsm + list(aaidsm[valid])
                AOD = AOD + list(aodmod[aodmod.mask==0][valid].data[:])
                
                AAIomiall = AAIomiall + list(aai)
                AAIdsmall = AAIdsmall + list(aaidsm)
                AODall = AODall + list(aodmod[aodmod.mask==0].data[:])
                
                
                
plt.figure(figsize = (10.5,4.5)) 
gs = gridspec.GridSpec(1,3, width_ratios=[1.25,0.001,1.25]) 

plt.subplot(gs[0])
slope4, intercept4, r_value4, p_value, std_err = stats.linregress(AODall,AAIomiall)    
plt.plot(np.arange(-25,25),np.arange(-25,25)*slope4+intercept4, linestyle = '-', linewidth =1, color = colorid[5], label = 'k = %1.2f b = %1.2f R = %1.2f' % (slope4, intercept4, r_value4))            
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(AOD,AAIomi)    
plt.plot(np.arange(-25,25),np.arange(-25,25)*slope3+intercept3, linestyle = '-', linewidth =1, color = colorid[0], label = 'Data in IQR \n k = %1.2f b = %1.2f R = %1.2f' % (slope3, intercept3, r_value3))      
plt.scatter(AODall,AAIomiall,s = 12, color = colorid[5], marker='s',linewidth =1, facecolor ='none')
plt.scatter(AOD,AAIomi,s = 12, color = colorid[0], linewidth =1, facecolor ='none')
plt.xlim(0,4)
plt.ylim(0,10)
plt.xticks(np.arange(0,4.1,1))
plt.yticks(np.arange(0,10.1,2))
plt.xlabel(r'$\tau_{aer}$')
plt.ylabel('OMAERO AAI')
plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(0.25,0.15))   


plt.subplot(gs[2])
slope3, intercept3, r_value3, p_value, std_err = stats.linregress(AODall,AAIdsmall)    
plt.plot(np.arange(-25,25),np.arange(-25,25)*slope3+intercept3, linestyle = '-', linewidth =1, color = colorid[5], label = 'k = %1.2f b = %1.2f R = %1.2f' % (slope3, intercept3, r_value3))      
slope4, intercept4, r_value4, p_value, std_err = stats.linregress(AOD,AAIdsm)    
plt.plot(np.arange(-25,25),np.arange(-25,25)*slope4+intercept4, linestyle = '-', linewidth =1, color = colorid[0], label = 'Data in IQR \n k = %1.2f b = %1.2f R = %1.2f' % (slope4, intercept4, r_value4))    
plt.scatter(AODall,AAIdsmall,s = 12, marker='s', linewidth = 1, color = colorid[5], facecolor ='none')
plt.scatter(AOD,AAIdsm,s = 12, linewidth = 1, color = colorid[0], facecolor ='none')
plt.xlim(0,4)
plt.ylim(0,10)
plt.xticks(np.arange(0,4.1,1))
plt.yticks(np.arange(0,10.1,2))
plt.xlabel(r'$\tau_{aer}$')
plt.ylabel('DISAMAR AAI')
plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(0.25,0.15))   
#plt.savefig(figdir + 'AAI_sim_Mie_stats_AOD_AAI_MOD.png', dpi = 300, transparent = False) 

                    
#plt.close('all')
                    

                    
                    
