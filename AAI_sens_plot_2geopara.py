"""
Plot AAI sensitivity for combined geometry

SZA + VZA for different VAA (HG + Mie)

@author: sunj
"""


######### load python packages
import os
import numpy as np
import matplotlib.pyplot as plt
import tables
from math import sqrt
import seaborn as sns
from DSM_Mie_func import mphyMie

plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':14})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)

# directory 
maindir = '/usr/people/sunj/Documents/DISAMAR'
outputdir = '/nobackup/users/sunj/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'

casename = 'AAI_sensitivity/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir += casename


# initialization 
aaiHG = [] 
aaiMie = []

ssaval = np.arange(0.7,0.96,0.05)             # single scattering albedo
ssaval = [0.75] 

Angval = np.arange(0.5,2.1,0.5)               # angstrom exponent
Angval = [1.]

asyval = np.array([0.]+list(np.arange(0.5,1.0,0.1)))                       # asymmetry factor
asyval = [0.]

rgval = np.arange(0.1,0.21,0.02)
rgval = [0.14] 

reval = np.arange(1.300, 1.501, 0.05)
reval = [1.5]

imval = np.arange(0.01, 0.06, 0.01)
imval = [0.06]      # geometry  

szaval = np.arange(0.,76.,15)
#szaval = [30.]

saaval = np.arange(0.,181.,45)
saaval = [0.]

vzaval = np.arange(0.,76.,15)
#vzaval = [0.]

vaaval = np.arange(0.,181.,45)
vaaval = [0]

aodval = np.arange(0.5,2.1,0.5)
aodval = [1.]

alhval = np.arange(2.5,8.6,2)
alhval = [4.5]

altval = np.arange(0.5,2.1,0.5)
altval = [1.]

#dpval = np.arange(8,25,4)    

"""

NOTICE: there are 2 approaches to calculate the scattering angle (Theta); 
Also depends how to define the relative azimuth angle.
  


""" 


for scheme in ['HG','Mie']: 
    expname = 'AAI_sens_'+scheme+'_geom_4' 
    
    sca1 = []       # MODIS & DISAMAR scattering angle calculation 
    sca2 = []       # others scattering angle calculation  
    
    for iAng in Angval:
        for iasy in asyval:
            for issa in ssaval: 
                for irg in rgval: 
                    for iim in imval: 
                        for ire in reval:
                            for isza in szaval:
                                for isaa in saaval:
                                    for ivza in vzaval: 
                                        for ivaa in vaaval: 
                                            for iaod in aodval:
                                                for ialt in altval: 
                                                    for ialh in alhval: 
    
                                                        # calculate cosine of zenith angle 
                                                        mu = np.cos(isza / 180. * np.pi)
                                                        mu0 = np.cos(ivza / 180. * np.pi)
                                                        
                                                        # MODIS & DISAMAR scattering angle calculation 
                                                        iraa = (isaa + 180 - ivaa) / 180 * np.pi
                                                        theta1 = np.arccos(-mu*mu0 + np.sqrt(1-mu**2)*np.sqrt(1-mu0**2)*np.cos(iraa)) / np.pi*180                           
                                                        sca1.append(theta1)
                                                        
                                                        # others scattering angle calculation 
                                                        raa = (isaa - ivaa) / 180 * np.pi
                                                        theta2 = np.arccos(mu*mu0 + np.sqrt(1-mu**2)*np.sqrt(1-mu0**2)*np.cos(raa)) / np.pi*180
                                                        sca2.append(theta2)
                                                        
                                                        if scheme == 'HG':
                                                            src = '/AAI_sens_saa-%1.2f_Ang-%1.2f_asy-%1.2f_sza-%03i_saa-%03i_vza-%03i_vaa-%03i_aod-%1.2f_alh-%1.2f_alt-%1.2f_dp-%02i.h5' \
                                                                              %(issa,iAng,iasy,isza,isaa,ivza,ivaa,iaod,ialh,ialt,8)
                                                            disamarhdf = tables.open_file(outputdir + expname + src,'r')
                                                            
                                                            disamarout = disamarhdf.root                    
                                                            aairetr = disamarout.parameters.AAI[0]
                                                            aaiHG.append(aairetr)
                                                            disamarhdf.close()          
                                                        
                                                        if scheme == 'Mie':
                                                            src = '/AAI_sens_rg-%1.2f_re-%1.4f_im-%1.4f_sza-%03i_saa-%03i_vza-%03i_vaa-%03i_aod-%1.2f_alh-%1.2f_alt-%1.2f_dp-%02i.h5' \
                                                                            %(irg,ire,iim,isza,isaa,ivza,ivaa,iaod,ialh,ialt,20)
                                                            disamarhdf = tables.open_file(outputdir + expname + src,'r')
                                                            
                                                            disamarout = disamarhdf.root                    
                                                            aairetr = disamarout.parameters.AAI[0]
                                                            aaiMie.append(aairetr)
                                                            disamarhdf.close()

                                        
# choose 2 parameters 
para1 = szaval
para2 = vzaval

aaiHGre = np.array(aaiHG).reshape(vzaval.shape[0],szaval.shape[0])
aaiMiere = np.array(aaiMie).reshape(vzaval.shape[0],szaval.shape[0])
sca1 = np.array(sca1).reshape(vzaval.shape[0],szaval.shape[0])       # MODIS & DISAMAR scattering angle calculation  
sca2 = np.array(sca2).reshape(vzaval.shape[0],szaval.shape[0])       # others scattering angle calculation 

# plot sensitivity    
fig1 = plt.figure(figsize=(10.5,4.5))   

v1 = np.arange(-0, 1.01, 0.025) 
v2 = np.arange(-0, 1.01, 0.02)
cbax = fig1.add_axes([0.925,0.1,0.025,0.8])  

plt.subplot(1,2,1)
P2, P1 = np.meshgrid(para2, para1)
plt.contourf(P2,P1,aaiHGre / abs(aaiHGre).max(), v1, cmap = 'rainbow')
#cbar = plt.colorbar(ticks = np.arange(0, np.max((aaiHG,aaiMie))+0.2,1), cax = cbax)
cbar = plt.colorbar(ticks = np.arange(0,1.1,1), cax = cbax)
cbar.set_label('Normalized AAI')
cs = plt.contour(P2, P1, sca1, colors = 'k', linestyles = 'dashed') 
plt.clabel(cs, colors = 'k', inline = 1, fontsize = 10, fmt = '%3i')
plt.plot([0,75],[0,75],'w--')
plt.title(r'HG $\phi$ = '+str(vaaval[0])+r' [$\circ$]')
plt.xlabel(r'$\theta$ [$\circ$]')
plt.ylabel(r'$\theta_0$ [$\circ$]')


plt.subplot(1,2,2)
P2, P1 = np.meshgrid(para2, para1)
plt.contourf(P2,P1,aaiMiere / aaiMiere.max(), v2, cmap = 'rainbow')
cbar.set_label('Normalized AAI')
cs = plt.contour(P2, P1, sca1, colors = 'k', linestyles = 'dashed') 
plt.clabel(cs, colors = 'k', inline = 1, fontsize = 10, fmt = '%3i')
plt.plot([0,75],[0,75],'w--')
plt.title(r'Mie $\phi$ = '+str(vaaval[0])+r' [$\circ$]')
plt.xlabel(r'$\theta$ [$\circ$]')
plt.ylabel(r'$\theta_0$ [$\circ$]')
#plt.savefig(figdir + 'AAI_sens_geom_g-%1.2f_vaa-%03i.png' % (iasy, ivaa), dpi = 300, transparent = True)  


expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (irg, ire,iim)
wvl, ssamie, ang, scamat = mphyMie(expfilename)

fig2 = plt.figure(figsize=(6,3))  
ax = fig2.add_axes([0.15,0.225,0.5,0.7])
axes = [ax, ax.twinx()]
fig2.subplots_adjust(right=0.75)
axes[-1].set_frame_on(True)
axes[-1].patch.set_visible(False)
axes[0].plot(np.diag(sca1[:-2]), np.diag(aaiHGre[:-2] / aaiHGre.max()), color = colorid[0], marker='s',label = 'HG')
axes[0].plot(np.diag(sca1), np.diag(aaiMiere / aaiMiere.max()), color = colorid[1], marker='s', label = 'Mie')
axes[1].plot(ang, scamat[np.where(wvl==354)[0][0],1:] / scamat[np.where(wvl==388)[0][0],1:], color = colorid[2], marker='s', markevery =20)
#axes[0].set_ylim(0.8,1.)
axes[0].set_xlabel('$\Theta [\circ]$')
axes[0].set_xticks(np.arange(0,181,30))
axes[0].set_ylabel('Normalized AAI', color = colorid[0])
axes[1].set_ylabel('p($\Theta$) contrast', color = colorid[2])
axes[0].set_yticks(np.arange(0.8,1.01,0.05))
axes[1].set_yticks(np.arange(0.9,1.21,0.1))
axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[1].tick_params(axis = 'y',colors=colorid[2])
ax.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1.2,0.5))   
#plt.savefig(figdir + 'AAI_sens_geom_mphy.png', dpi = 300, transparent = True)  


#plt.close('all')
