"""
Plot AAI sensitivity for combined aerosol relevant parameters

SSA + Ang   (HG)
SSA + g     (HG)
Ang + g     (HG)
AOD + ALH   (HG + Mie)
ALH + ALT   (HG + Mie)
nr + ni     (Mie)

@author: sunj
"""


######### load python packages
import os
import numpy as np
import matplotlib.pyplot as plt
import tables
import seaborn as sns
from DSM_Mie_func import mphyMie 

plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':14})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)


# directory 
maindir = '/usr/people/sunj/Documents/DISAMAR/disamar'
outputdir = '/nobackup/users/sunj/'
expcoefdir = maindir + '/expCoefFiles/'

figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'

casename = 'AAI_sensitivity/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir += casename

# initialization 


ssaval = np.arange(0.7,0.96,0.05)             # single scattering albedo
ssaval = [0.75] 

Angval = np.arange(0,2.1,0.5)               # angstrom exponent
Angval = [1.]

asyval = np.array([0.]+list(np.arange(0.5,1.0,0.1)))                       # asymmetry factor
asyval = [0.8]

#rgval = np.arange(0.1,0.21,0.05)
rgval = np.arange(0.1,0.41,0.05)
#rgval = [0.05,0.1,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.30,0.35,0.40] 

reval = np.arange(1.300, 1.501, 0.05)
#reval = [1.40]

imval = np.arange(0.04, 0.11, 0.02)
#imval = [0.06] 

szaval = np.arange(0.,76.,15)
szaval = [30.]

saaval = np.arange(0.,181.,45)
saaval = [0.]

vzaval = np.arange(0.,76.,15)
vzaval = [0.]

vaaval = np.arange(0.,181.,45)
vaaval = [180.]

aodval = np.arange(0.5,2.1,0.5)
aodval = [1.]

alhval = np.arange(2.5,8.6,2)
alhval = [4.5]

altval = np.arange(0.5,2.1,0.5)
altval = [1.]

#dpval = np.arange(8,25,4)    

wvl = np.arange(340,677,2)

"""

NOTICE: put para2 (x-axis) loop should be after para1 (legend) loop. The aai list is reshapred as [para1, para2], m-by-n
para1(1): para2(1), para2(2), ... , para2(n)
para1(2): para2(1), para2(2), ... , para2(n)
  :     :    :    ,    :    ,  :  ,      : 
para1(m): para2(1), para2(2), ... , para2(n)

""" 
aaiMie = np.zeros([len(rgval),len(imval),len(reval)]) 
refMie = np.zeros([len(rgval),len(imval),len(reval),2]) 
ssauv = np.zeros([len(rgval),len(imval),len(reval),2]) 
guv = np.zeros([len(rgval),len(imval),len(reval),2]) 

Robs = np.zeros([len(rgval),len(imval),len(reval)]) 
Rray = np.zeros([len(rgval),len(imval),len(reval)]) 
Rref = np.zeros([len(rgval),len(imval),len(reval)]) 



for scheme in ['Mie']: 
    
    expname = 'AAI_sens_'+scheme+'_mphy_6' 
    
    for iasy in asyval: 
        for iAng in Angval:
            for issa in ssaval: 
                for x, irg in enumerate(rgval): 
                     for y, iim in enumerate(imval):
                         for z, ire in enumerate(reval):
                            for isza in szaval:
                                for isaa in saaval:
                                    for ivza in vzaval: 
                                        for ivaa in vaaval: 
                                            for iaod in aodval:
                                                for ialt in altval: 
                                                    for ialh in alhval: 
    
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
                                                            rad = disamarout.additional_output_sim.earth_radiance_band_1[:] 
                                                            irr = disamarout.additional_output_sim.solar_irradiance_band_1[:]
                                                            refMie[x,y,z,0] = ( rad[0] / irr[0] )
                                                            refMie[x,y,z,1] = ( rad[1] / irr[1] )
                                                            Robs[x,y,z] = (disamarout.sun_normalized_radiance.measured_spectrum[0])
                                                            Rray[x,y,z] = (disamarout.sun_normalized_radiance.fitted_spectrum[0])
                                                            Rref[x,y,z] = (disamarout.sun_normalized_radiance.fitted_spectrum[1])
                                                            aaiMie[x,y,z] = aairetr
                                                            disamarhdf.close()
#                                                            print isza, ivza, ire, aairetr
                                                            
                                                            
                                                        expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (irg,ire,iim)
                                                        expfile = expcoefdir + expfilename + '.dat'
                                                        with open(expfile,'r') as f: 
                                                            content = f.readlines()
                                                        
                                                        gmie = [] 
                                                        for ig, il in enumerate(content): 
                                                            if il.find('alpha1') != -1:
                                                                temp1 = content[ig+2]
                                                                temp2 = float(temp1.split('  ')[2]) / 3.
                                                                gmie.append(temp2)
                                                                
                                                        guv[x,y,z,0] = gmie[np.where(wvl == 354)[0][0]]
                                                        guv[x,y,z,1] = gmie[np.where(wvl == 388)[0][0]]
                                                        
                                                        ssamie = [] 
                                                        for ig, il in enumerate(content): 
                                                            if il.find('ssa') != -1:
                                                                ssamie.append(float(il[5:11]))
                                                        ssauv[x,y,z,0] = ssamie[np.where(wvl==354)[0][0]]
                                                        ssauv[x,y,z,1] = ssamie[np.where(wvl==388)[0][0]]
                                                        
#print 'ssa', ssauv[:,:,:,0]
#print 'g', guv[:,:,:,0]
        
      
para1 = rgval        
para2 = imval
para3 = reval


temp1 = []
for irg in rgval: 
    expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (irg,1.5,0.06)
    wvl, ssamie, ang, scamat = mphyMie(expfilename)
    temp1.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]] / scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 

fig1 = plt.figure(figsize=(6,3))  
ax = fig1.add_axes([0.125,0.225,0.5,0.7])
axes = [ax, ax.twinx()]
fig1.subplots_adjust(right=0.75)
#axes[-1].spines['right'].set_position(('axes', 1.5))
#axes[-2].spines['right'].set_position(('axes', 1.25))
#axes[-2].set_frame_on(True)
#axes[-2].patch.set_visible(False)
axes[0].plot(para1, ssauv[:,1,-1,0], color = colorid[5], marker='s')
axes[0].plot(para1, ssauv[:,1,-1,1], color = colorid[5], linestyle= '--',marker='s')
axes[1].plot(para1, guv[:,1,-1,0], color = colorid[0], marker='s')
axes[1].plot(para1, guv[:,1,-1,1], color = colorid[0], linestyle= '--', marker='s')
#axes[3].plot(para1, temp1, color = 'k', linestyle= '-', marker='s')
#axes[0].plot(para1, aaiMie[:,1,-1] / aaiMie[:,1,-1].max(), color = colorid[0], marker='s')
axes[0].set_xlabel('r$_g$ [$\mu m$]')
axes[0].set_xticks(para1)
axes[0].set_xlim(0.1,0.40000001)
#axes[0].set_ylabel('Normalized AAI', color = colorid[0])
axes[0].set_ylabel('$\omega_0$', color = colorid[5])
axes[1].set_ylabel('g', color = colorid[0])
#axes[3].set_ylabel(r'P(150$^{\circ}$) 354 / 388', color = 'k')
#axes[0].set_yticks(np.arange(0.8,1.1,0.1))
#axes[0].set_ylim(0.8,1)
axes[0].set_yticks(np.arange(0.5,0.91,0.1))
axes[0].set_ylim(0.5,0.9)
axes[1].set_yticks(np.arange(0.6,1.1,0.1))
axes[1].set_ylim(0.6,1)
#axes[3].set_ylim(0.9,1)
#axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[0].tick_params(axis = 'y',colors=colorid[5])
axes[1].tick_params(axis = 'y',colors=colorid[0])
plt.savefig(figdir + 'AAI_sens_mie_rg.png', dpi = 300, transparent = True)


fig2 = plt.figure(figsize=(6,3))
ax = fig2.add_axes([0.15,0.225,0.5,0.7])  
axes = [ax, ax.twinx()]
fig2.subplots_adjust(right=0.75)
#axes[-1].spines['right'].set_position(('axes', 1.25))
#axes[-1].set_frame_on(True)
#axes[-1].patch.set_visible(False)
#axes[0].plot(para2, aaiMie[1,:,-1] / aaiMie[1,:,-1].max(), color = colorid[0], marker='s')
axes[0].plot(para2, ssauv[1,:,-1,0], color = colorid[5], marker='s')
axes[0].plot(para2, ssauv[1,:,-1,1], color = colorid[5], linestyle = '--', marker='s')
axes[1].plot(para2, guv[1,:,-1,0], color = colorid[0], marker='s')
axes[1].plot(para2, guv[1,:,-1,1], color = colorid[0], linestyle = '--', marker='s')
axes[0].set_xlabel('n$_i$ [-]')
axes[0].set_xticks(para2)
#axes[0].set_ylabel('Normalized AAI', color = colorid[0])
axes[0].set_ylabel('$\omega_0$', color = colorid[5])
axes[1].set_ylabel('g', color = colorid[0])
#axes[0].set_yticks(np.arange(0.6,1.1,0.1))
#axes[0].set_ylim(0.6,1)
axes[0].set_yticks(np.arange(0.5,0.91,0.1))
axes[0].set_ylim(0.5,0.9)
axes[1].set_yticks(np.arange(0.6,1.1,0.1))
axes[1].set_ylim(0.6,1)
#axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[0].tick_params(axis = 'y',colors=colorid[5])
axes[1].tick_params(axis = 'y',colors=colorid[0])
plt.savefig(figdir + 'AAI_sens_mie_ni.png', dpi = 300, transparent = True)


fig3 = plt.figure(figsize=(6,3))  
ax = fig3.add_axes([0.15,0.225,0.5,0.7])
axes = [ax, ax.twinx()]
fig3.subplots_adjust(right=0.75)
#axes[-1].spines['right'].set_position(('axes', 1.25))
#axes[-1].set_frame_on(True)
#axes[-1].patch.set_visible(False)
#axes[0].plot(para3, aaiMie[1,1,:] / aaiMie[1,1,:].max(), color = colorid[0], marker='s')
axes[0].plot(para3, ssauv[1,1,:,0], color = colorid[5], marker='s')
axes[0].plot(para3, ssauv[1,1,:,1], color = colorid[5], linestyle = '--', marker='s')
axes[1].plot(para3, guv[1,1,:,0], color = colorid[0], marker='s')
axes[1].plot(para3, guv[1,1,:,1], color = colorid[0], linestyle = '--', marker='s')
axes[0].set_xlabel('n$_r$ [-]')
axes[0].set_xticks(para3)
axes[0].set_xlim(1.3,1.500000001)
#axes[0].set_ylabel('Normalized AAI', color = colorid[0])
axes[0].set_ylabel('$\omega_0$', color = colorid[5])
axes[1].set_ylabel('g', color = colorid[0])
#axes[0].set_yticks(np.arange(0.7,1.1,0.1))
#axes[0].set_ylim(0.7,1)
axes[0].set_yticks(np.arange(0.5,0.91,0.1))
axes[0].set_ylim(0.5,0.9)
axes[1].set_yticks(np.arange(0.6,1.1,0.1))
axes[1].set_ylim(0.6,1)
#axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[0].tick_params(axis = 'y',colors=colorid[5])
axes[1].tick_params(axis = 'y',colors=colorid[0])
plt.savefig(figdir + 'AAI_sens_mie_nr.png', dpi = 300, transparent = True)




temp1 = []
temp2 = []
temp3 = [] 
for irg in rgval:
    print irg
    expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (irg,1.5,0.06)
    wvl, ssamie, ang, scamat = mphyMie(expfilename)
    temp1.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]] / scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
    temp2.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]])
    temp3.append(scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 


fig4 = plt.figure(figsize=(6,3))  
ax = fig4.add_axes([0.15,0.225,0.5,0.7])
axes = [ax, ax.twinx(), ax.twinx()]
fig4.subplots_adjust(right=0.75)
axes[-1].spines['right'].set_position(('axes', 1.3))
#axes[-1].spines['right'].set_position(('axes', 1.5))
#axes[-1].set_frame_on(True)
#axes[-1].patch.set_visible(False)
#axes[3].plot(para1, temp1, color = 'k', marker='s')
axes[0].plot(para1,aaiMie[:,1,-1] / aaiMie[:,1,-1].max() , color = 'k', marker ='s')
axes[1].plot(para1,(Rray[:,1,-1] - Robs[:,1,-1]) / (Rray[:,1,-1] - Robs[:,1,-1]).max(), color = colorid[5], marker ='s')
axes[2].plot(para1,(Robs[:,1,-1])/(Robs[:,1,-1]).max(), color = colorid[0], marker ='s')


#axes[1].plot(para1, temp2, color = colorid[6], marker='s')
#axes[1].plot(para1, temp3, color = colorid[6], linestyle= '--', marker='s')
axes[0].set_xlabel('r$_g$ [$\mu m$]')
axes[0].set_xlim(0.1,0.40000001)
axes[0].set_ylabel('Normalized AAI', color = 'k')
axes[1].set_ylabel(r'Normalized $\Delta I_{\lambda_1}$', color = colorid[5])
axes[2].set_ylabel(r'Normalized $I^{obs}_{\lambda_1}$', color = colorid[0])
#axes[3].set_ylabel(r'P(150$^{\circ}$) 354/388', color = 'k')
axes[0].set_yticks(np.arange(0.8,1.01,0.1))
axes[0].set_ylim(0.8,1)
axes[1].set_yticks(np.arange(0.7,1.01,0.1))
axes[1].set_ylim(0.7,1)
axes[2].set_yticks(np.arange(0.8,1.01,0.1))
axes[2].set_ylim(0.8,1)
#axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[1].tick_params(axis = 'y',colors=colorid[5])
axes[2].tick_params(axis = 'y',colors=colorid[0])

plt.savefig(figdir + 'AAI_sens_mie_rg_analysis.png', dpi = 300, transparent = True)







fig5 = plt.figure(figsize=(6,3))  
ax = fig5.add_axes([0.15,0.225,0.5,0.7])
axes = [ax, ax.twinx(), ax.twinx()]
fig5.subplots_adjust(right=0.75)
axes[-1].spines['right'].set_position(('axes', 1.3))
#axes[-1].spines['right'].set_position(('axes', 1.5))
axes[-1].set_frame_on(True)
axes[-1].patch.set_visible(False)
#axes[3].plot(para1, temp1, color = 'k', marker='s')
axes[0].plot(para2,aaiMie[1,:,-1] / aaiMie[1,:,-1].max() , color = 'k', marker ='s')
axes[1].plot(para2,(Rray[1,:,-1] - Robs[1,:,-1]) / (Rray[1,:,-1] - Robs[1,:,-1]).max(), color = colorid[5], marker ='s')
axes[2].plot(para2,(Robs[1,:,-1]) / Robs[1,:,-1].max(), color = colorid[0], marker ='s')

#axes[1].plot(para1, temp2, color = colorid[6], marker='s')
#axes[1].plot(para1, temp3, color = colorid[6], linestyle= '--', marker='s')
axes[0].set_xlabel('n$_i$')
#axes[0].set_xlim(0.1,0.40000001)
axes[0].set_ylabel('Normalized AAI', color = 'k')
axes[1].set_ylabel(r'Normalized $\Delta I_{\lambda1}$', color = colorid[5])
axes[2].set_ylabel(r'Normalized $I^{obs}_{\lambda1}$', color = colorid[0])
#axes[3].set_ylabel(r'P(150$^{\circ}$) 354/388', color = 'k')
axes[0].set_yticks(np.arange(0.6,1.01,0.1))
axes[0].set_ylim(0.6,1)
axes[1].set_yticks(np.arange(0.7,1.01,0.1))
axes[1].set_ylim(0.7,1)
axes[2].set_yticks(np.arange(0.8,1.01,0.1))
axes[2].set_ylim(0.8,1)
#axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[1].tick_params(axis = 'y',colors=colorid[5])
axes[2].tick_params(axis = 'y',colors=colorid[0])

plt.savefig(figdir + 'AAI_sens_mie_ni_analysis.png', dpi = 300, transparent = True)




fig6 = plt.figure(figsize=(6,3))  
ax = fig6.add_axes([0.15,0.225,0.5,0.7])
axes = [ax, ax.twinx(), ax.twinx()]
fig6.subplots_adjust(right=0.75)
axes[-1].spines['right'].set_position(('axes', 1.3))
#axes[-1].spines['right'].set_position(('axes', 1.5))
axes[-1].set_frame_on(True)
axes[-1].patch.set_visible(False)
#axes[3].plot(para1, temp1, color = 'k', marker='s')
axes[0].plot(para3,aaiMie[1,1,:] / aaiMie[1,1,:].max() , color = 'k', marker ='s')
axes[1].plot(para3,(Rray[1,1,:] - Robs[1,1,:]) / (Rray[1,1,:] - Robs[1,1,:]).max(), color = colorid[5], marker ='s')
axes[2].plot(para3,(Robs[1,1,:])/(Robs[1,1,:]).max() , color = colorid[0], marker ='s')
#axes[1].plot(para1, temp2, color = colorid[6], marker='s')
#axes[1].plot(para1, temp3, color = colorid[6], linestyle= '--', marker='s')
axes[0].set_xlabel('n$_r$')
axes[0].set_xlim(1.3,1.5)
axes[0].set_ylabel('Normalized AAI', color = 'k')
axes[1].set_ylabel(r'Normalized $\Delta I_{\lambda1}$', color = colorid[5])
axes[2].set_ylabel(r'Normalized $I^{obs}_{\lambda1}$', color = colorid[0])
#axes[3].set_ylabel(r'P(150$^{\circ}$) 354/388', color = 'k')
axes[0].set_yticks(np.arange(0.6,1.01,0.1))
axes[0].set_ylim(0.6,1)
axes[1].set_yticks(np.arange(0.7,1.01,0.1))
axes[1].set_ylim(0.7,1)
axes[2].set_yticks(np.arange(0.8,1.01,0.1))
axes[2].set_ylim(0.8,1)
#axes[0].tick_params(axis = 'y',colors=colorid[0])
axes[1].tick_params(axis = 'y',colors=colorid[5])
axes[2].tick_params(axis = 'y',colors=colorid[0])

plt.savefig(figdir + 'AAI_sens_mie_nr_analysis.png', dpi = 300, transparent = True)


#temp1 = []
#temp2 = []
#temp3 = [] 
#for iim in imval: 
#    expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (0.14,1.5,iim)
#    wvl, ssamie, ang, scamat = mphyMie(expfilename)
#    temp1.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]] / scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
#    temp2.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]])
#    temp3.append(scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
#
#
#fig4 = plt.figure(figsize=(6,3))  
#ax = fig4.add_axes([0.15,0.225,0.5,0.7])
#axes = [ax, ax.twinx()]
#fig4.subplots_adjust(right=0.75)
##axes[-1].spines['right'].set_position(('axes', 1.25))
#axes[-1].set_frame_on(True)
#axes[-1].patch.set_visible(False)
#axes[0].plot(para2, temp1, color = 'k', marker='s')
#axes[1].plot(para2, temp2, color = colorid[6], marker='s')
#axes[1].plot(para2, temp3, color = colorid[6], linestyle= '--', marker='s')
#axes[0].set_xlabel('r$_g$ [$\mu m$]')
#axes[0].set_xlim(0.1,0.40000001)
#axes[0].set_ylabel('P(150$^{\circ}$) 354 / 388', color = 'k')
#axes[1].set_ylabel('P(150$^{\circ}$)', color = colorid[6])
#axes[0].set_yticks(np.arange(0.8,1.01,0.1))
#axes[0].set_ylim(0.8,1)
#axes[1].set_yticks(np.arange(0.055,0.076,0.005))
#axes[1].set_ylim(0.055,0.075)
#axes[0].tick_params(axis = 'y',colors='k')
#axes[1].tick_params(axis = 'y',colors=colorid[6])

#temp1 = []
#temp2 = []
#temp3 = [] 
#for ire in reval: 
#    expfilename = 'rg-%1.2f_nr-%1.4f_ni-%1.4f' % (0.14,ire,0.06)
#    wvl, ssamie, ang, scamat = mphyMie(expfilename)
#    temp1.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]] / scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
#    temp2.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]])
#    temp3.append(scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
#
#
#fig4 = plt.figure(figsize=(6,3))  
#ax = fig4.add_axes([0.15,0.225,0.5,0.7])
#axes = [ax, ax.twinx()]
#fig4.subplots_adjust(right=0.75)
##axes[-1].spines['right'].set_position(('axes', 1.25))
#axes[-1].set_frame_on(True)
#axes[-1].patch.set_visible(False)
#axes[0].plot(para3, temp1, color = 'k', marker='s')
#axes[1].plot(para3, temp2, color = colorid[6], marker='s')
#axes[1].plot(para3, temp3, color = colorid[6], linestyle= '--', marker='s')


#plt.close('all')  