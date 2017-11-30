"""
Calculate AAI uncertainty

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
import shutil
import glob
import numpy as np
import matplotlib.pyplot as plt
import tables
import seaborn as sns
from scipy import stats
from scipy.optimize import fsolve

plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':14})
sns.set_style('white')
colorid = sns.color_palette("hls", 10)
#colorid = [colorid[0],colorid[1],colorid[3],colorid[5]]


# directory 
maindir = '/usr/people/sunj/Documents/DISAMAR/'
outputdir = maindir + 'pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'

exedir = maindir + 'disamar/createExpCoefFiles/'
expindir = exedir + 'inputFiles/'
expoutdir = exedir + 'expCoefFiles/'
os.chdir(exedir)

expcoefdir = maindir + 'disamar/expCoefFiles/'
template = expindir + 'smoke_template.in'


casename = 'AAI_sensitivity/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir += casename

# initialization 
aaiHG = [] 
aaiMie = [] 

sigma = 1.5

ssaval = np.arange(0.7,0.96,0.05)             # single scattering albedo
ssaval = [0.75] 

Angval = np.arange(0.5,2.1,0.5)               # angstrom exponent
Angval = [1.]

asyval = np.arange(0.5,1.0,0.1)                # asymmetry factor
asyval = [0.5]

rgval = np.arange(0.1,0.21,0.02)
rgval = [0.14] 

reval = np.arange(1.400, 1.601, 0.05)
reval = [1.5]

imval = np.arange(0.01, 0.09, 0.01)
#imval = [0.03]

szaval = np.arange(0.,76.,15)
szaval = [30.]

saaval = np.arange(0.,181.,45)
saaval = [0.]

vzaval = np.arange(0.,76.,15)
vzaval = [0.]

vaaval = np.arange(0.,181.,45)
vaaval = [180.]

#aodval = np.arange(0.5,2.1,0.5)
aodval = np.arange(0.2,0.9,0.1)
#aodval = [0.5]

alhval = np.arange(3,5.1,0.25)
#alhval = [4]

altval = np.arange(0.25,2.1,0.25)
#altval = [1]

#dpval = np.arange(8,25,4)    


# the following parameters using plume mean value
wvl = np.arange(340,677,2)
ssaval = np.arange(0.7,1.01,0.05)

aaitemp = [] 
aaiMie = np.ones([len(reval), len(imval), len(aodval), len(alhval), len(altval)]) * np.nan 
#ni_est = np.ones([len(reval), len(imval), len(aodval), len(alhval), len(altval)]) * np.nan 
ssa = np.ones([len(imval)]) * np.nan 
#aai_ssa = np.ones([len(ssaval), len(aodval), len(alhval), len(altval)]) * np.nan 


for scheme in ['Mie']: 
    
    expname = 'AAI_sens_'+scheme+'_uncertainty_4' 
    filelist = glob.glob(nobackupdir + expname + '/*')
    for iAng in Angval:
        for issa in ssaval: 
            for iasy in asyval: 
                for r,irg in enumerate(rgval): 
                    for m,iim in enumerate(imval):
                        for ire in reval:
                            for isza in szaval:
                                for isaa in saaval:
                                    for ivza in vzaval: 
                                        for ivaa in vaaval: 
                                            for i,iaod in enumerate(aodval):
                                                for j,ialh in enumerate(alhval): 
                                                    for k,ialt in enumerate(altval):
                                                        src = '/AAI_sens_rg-%1.2f_re-%1.4f_im-%1.4f_sza-%03i_saa-%03i_vza-%03i_vaa-%03i_aod-%1.2f_alh-%1.2f_alt-%1.2f_dp-%02i.h5' \
                                                                        %(irg, ire,iim,isza,isaa,ivza,ivaa,iaod,ialh,ialt,20)
                                                        if (nobackupdir + expname + src) in filelist: 
                                                            disamarhdf = tables.open_file(nobackupdir + expname + src,'r')
                                                            
                                                            disamarout = disamarhdf.root                    
                                                            aairetr = disamarout.parameters.AAI[0]
                                                            aaiMie[r,m,i,j,k] = aairetr
                                                            disamarhdf.close() 
                                            

    # calculate uncertainty 
    
    p_aai = np.polyfit(imval,aaiMie[0,:,3,6,2],2)
    aai_im = p_aai[0] * imval**2 + p_aai[1] * imval + p_aai[2]


    for m, iim in enumerate(imval): 
        expfiledat = expcoefdir + 'rg-%1.2f_nr-%1.4f_ni-%1.4f.dat' % (irg, ire, iim)
        with open(expfiledat,'r') as f: 
            content = f.readlines()
    
        ssamie = [] 
        for ig, il in enumerate(content): 
            if il.find('ssa') != -1:
                ssamie.append(float(il[5:11]))
        ssamie = np.array(ssamie)
        
        ssa[m] = ssamie[np.where(wvl==550)[0][0]]  
   
    p_ssa = np.polyfit(ssa, imval,2)                    # ni as function of ssa
    ni_ssa = p_ssa[0] * ssaval**2 + p_ssa[1] * ssaval + p_ssa[2]     
    plt.figure()
    plt.plot(ssaval,ni_ssa)
    plt.xlabel('ssa')
    plt.ylabel('ni')
 

    
    
   

    aod_ssa = np.zeros([len(ssaval),len(aodval)]) 
#    plt.figure() 
    for i in range(len(aodval)): 
        p_aai = np.polyfit(imval,aaiMie[0,:,i,6,2],2)     # ni as function of aai 
        aod_ssa[:,i] = p_aai[0] * ni_ssa**2 + p_aai[1] * ni_ssa + p_aai[2]      # aai as function of ssa 
#        plt.plot(ssaval, aod_ssa[:,i], color = colorid[i],label=r'$\tau_{aer}$ = %1.1f' % (aodval[i]))
#    plt.xlabel(r'$\omega_0$ @ 550 nm')
#    plt.ylabel('AAI')
#    plt.legend() 
    
    fig1 = plt.figure(figsize = (8,4)) 
    ax1 = fig1.add_axes([0.15,0.225,0.5,0.7])
    for i in range(len(aodval)): 
        p_aai = np.polyfit(imval,aaiMie[0,:,i,6,2],2)     # ni as function of aai 
        if i != 3:
            plt.plot(ssaval - 0.85, aod_ssa[:,i] - aod_ssa[3,3], 's-', color = colorid[i],label=r'$\tau_{aer}$ = %1.1f' % (aodval[i]))
        if i == 3:
            plt.plot(ssaval - 0.85, aod_ssa[:,i] - aod_ssa[3,3], 's-', color = 'k',label=r'$\tau_{aer}$ = %1.1f' % (aodval[i]))
#    plt.axhline(0,color = 'k',linestyle = '--')
    plt.xlabel(r'$\omega_0$ @ 550 nm')
    plt.ylabel('AAI')
    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    plt.grid(True)
    plt.xlim(-0.15,0.15)
    plt.xticks(np.arange(-0.15,0.16,0.05))
    plt.savefig(figdir + 'AAI_ssa_uncertainty_aod.png', dpi = 300, transparent = True)
    
    
#    plt.figure()     
    alh_ssa = np.zeros([len(ssaval),len(alhval)]) 
    for i in range(len(alhval)): 
        p_aai = np.polyfit(imval,aaiMie[0,:,3,i,2],2)     # ni as function of aai 
        alh_ssa[:,i] = p_aai[0] * ni_ssa**2 + p_aai[1] * ni_ssa + p_aai[2]      # aai as function of alh 
#        plt.plot(ssaval, alh_ssa[:,i], color = colorid[i],label=r'$z_{aer}$ = %1.2f' % (alhval[i]))
#    plt.xlabel(r'$\omega_0$ @ 550 nm')
#    plt.ylabel('AAI')
#    plt.legend() 
        
    fig2 = plt.figure(figsize = (8,4)) 
    ax2 = fig2.add_axes([0.15,0.225,0.5,0.7])
    for i in range(len(alhval)): 
        p_aai = np.polyfit(imval,aaiMie[0,:,3,i,2],2)     # ni as function of aai 
        if i != 6:
            plt.plot(ssaval-0.85, alh_ssa[:,i] - alh_ssa[3,6], 's-', color = colorid[i],label=r'$z_{aer}$ = %1.2f [km]' % (alhval[i]))
        if i == 6:
            plt.plot(ssaval-0.85, alh_ssa[:,i] - alh_ssa[3,6], 's-', color = 'k',label=r'$z_{aer}$ = %1.2f [km]' % (alhval[i]))
#    plt.axhline(0,color = 'k',linestyle = '--')
    plt.xlabel(r'$\omega_0$ @ 550 nm')
    plt.ylabel('AAI')
    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    plt.grid(True)
    plt.xlim(-0.15,0.15)
    plt.xticks(np.arange(-0.15,0.16,0.05))
    plt.savefig(figdir + 'AAI_ssa_uncertainty_alh.png', dpi = 300, transparent = True)
    
    
    

#    plt.figure()  
    alt_ssa = np.zeros([len(ssaval),len(altval)]) 
    for i in range(len(altval)): 
        p_aai = np.polyfit(imval,aaiMie[0,:,3,6,i],2)     # ni as function of aai 
        alt_ssa[:,i] = p_aai[0] * ni_ssa**2 + p_aai[1] * ni_ssa + p_aai[2]      # aai as function of alt 
#        plt.plot(ssaval, alt_ssa[:,i], color = colorid[i],label=r'$\Delta z$ = %1.2f' % (altval[i]))
#    plt.xlabel(r'$\omega_0$ @ 550 nm')
#    plt.ylabel('AAI')
#    plt.legend() 
        
    fig3 = plt.figure(figsize = (8,4)) 
    ax3 = fig3.add_axes([0.15,0.225,0.5,0.7])
    for i in range(len(altval)): 
        p_aai = np.polyfit(imval,aaiMie[0,:,3,6,i],2)     # ni as function of aai 
        if i != 2:
            plt.plot(ssaval-0.85, alt_ssa[:,i]-alt_ssa[3,2], 's-', color = colorid[i],label=r'$\Delta z$ = %1.2f [km]' % (altval[i]))
        if i == 2: 
            plt.plot(ssaval-0.85, alt_ssa[:,i]-alt_ssa[3,2], 's-', color = 'k',label=r'$\Delta z$ = %1.2f [km]' % (altval[i]))
#    plt.axhline(0,color = 'k',linestyle = '--')
    plt.xlabel(r'$\omega_0$ @ 550 nm')
    plt.ylabel('AAI')
    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    plt.grid(True)
    plt.xlim(-0.15,0.15)
    plt.xticks(np.arange(-0.15,0.16,0.05))
    plt.savefig(figdir + 'AAI_ssa_uncertainty_alt.png', dpi = 300, transparent = True)
        


#    for i in range(len(aodval)): 
#        def aai_im(x):
#            return p_im[0] * x**2 + p_im[1] * x + p_im[2] - aaiMie[2,i,4,3]
#        ni_est[i,4,3] = fsolve(aai_im,0)
#    print ni_est[:,4,3]
#    
#    for j in range(len(alhval)): 
#        def aai_im(x):
#            return p_im[0] * x**2 + p_im[1] * x + p_im[2] - aaiMie[2,4,j,3]
#        ni_est[3,j,3] = fsolve(aai_im,0)
#    print ni_est[3,:,3]
#    
#    for k in range(len(altval)): 
#        def aai_im(x):
#            return p_im[0] * x**2 + p_im[1] * x + p_im[2] - aaiMie[2,4,4,k]
#        ni_est[3,4,k] = fsolve(aai_im,0)
#    print ni_est[3,4,:]
    
#    for m in range(len(imval)): 
#        for i in range(len(aodval)): 
#            for j in range(len(alhval)): 
#                for k in range(len(altval)): 
#                    if np.isnan(aaiMie[m,i,j,k]) == False:
#                        def aai_im(x):
#                            return p_aai[0] * x**2 + p_aai[1] * x + p_aai[2] - aaiMie[m,i,j,k]
#                        ni_est[m,i,j,k] = fsolve(aai_im, 0)
                        



#    p_ssa = np.polyfit(imval,ssa,2)
#    
#    ni = ni_est[:,3,4,3]
#    ssa_ni = p_ssa[0] * ni**2 + p_ssa[1] * ni + p_ssa[2]            
#    plt.figure(figsize=(4.5,3))
##    plt.subplot(3,1,1)    
##    plt.plot(imval,aaiMie[:,3,4,3],'s-')
##    plt.subplot(3,1,2)    
##    plt.plot(imval,ni_est[:,3,4,3],'s-')
##    plt.subplot(3,1,3)    
##    plt.plot(imval,ssa_ni,'s-')
#    plt.plot(imval[0:5], ssa_ni[0:5] - ssa_ni[2],'s-')
#    plt.axhline(0,color = 'k', linestyle = '--')
#    plt.xlabel('ni')
#    plt.ylabel(r'$\omega_{0}$')
#    plt.xlim(imval[0],imval[4])
#    plt.xticks(imval[0:5])
#    plt.ylim(-0.15,0.15)
#    plt.yticks(np.arange(-0.15,0.16,0.05))
#    
#    ni = ni_est[2,:,4,3]
#    ssa_aod = p_ssa[0] * ni**2 + p_ssa[1] * ni + p_ssa[2]                
#    plt.figure(figsize=(4.5,3))
##    plt.subplot(3,1,1)    
##    plt.plot(aodval,aaiMie[2,:,4,3],'s-')
##    plt.subplot(3,1,2)    
##    plt.plot(aodval,ni_est[2,:,4,3],'s-')
##    plt.subplot(3,1,3)    
##    plt.plot(aodval,ssa_ni,'s-')
#    plt.plot(aodval, ssa_aod - ssa_aod[3],'s-')
#    plt.axhline(0,color = 'k', linestyle = '--')
#    plt.xlabel(r'$\tau_{aer}$')
#    plt.ylabel(r'$\omega_{0}$')
#    plt.xticks(aodval)
#    plt.ylim(-0.15,0.15)
#    
#    
#    ni = ni_est[2,3,:,3]
#    ssa_alh = p_ssa[0] * ni**2 + p_ssa[1] * ni + p_ssa[2]     
#    plt.figure(figsize=(4.5,3))
##    plt.subplot(3,1,1)    
##    plt.plot(alhval,aaiMie[2,3,:,3],'s-')
##    plt.subplot(3,1,2)    
##    plt.plot(alhval,ni_est[2,3,:,3],'s-')
##    plt.subplot(3,1,3)    
##    plt.plot(alhval,ssa_ni,'s-')
#    plt.plot(alhval, ssa_alh - ssa_alh[4],'s-')
#    plt.axhline(0,color = 'k', linestyle = '--')
#    plt.xlabel(r'$z_{aer}$ [km]')
#    plt.ylabel(r'$\omega_{0}$')
#    plt.xticks(alhval)
#    plt.ylim(-0.1,0.1)
#
#
#
#    ni = ni_est[2,3,4,:]
#    ssa_alt = p_ssa[0] * ni**2 + p_ssa[1] * ni + p_ssa[2]     
#    plt.figure(figsize=(4.5,3))
##    plt.subplot(3,1,1)    
##    plt.plot(altval,aaiMie[2,3,4,:],'s-')
##    plt.subplot(3,1,2)    
##    plt.plot(altval,ni_est[2,3,4,:],'s-')
##    plt.subplot(3,1,3)    
##    plt.plot(altval,ssa_ni,'s-')
#    plt.plot(altval, ssa_alt - ssa_alt[3],'s-')
#    plt.axhline(0,color = 'k', linestyle = '--')
#    plt.xlabel(r'$\Delta z_{aer} [km]$')
#    plt.ylabel(r'$\omega_{0}$')
#    plt.xticks(altval)
#    plt.ylim(-0.3,0.3)

#fig1 = plt.figure(figsize=(6,5))  
#ax = fig1.add_axes([0.15,0.25,0.75,0.5])
#axes = [ax, ax.twiny(), ax.twiny(), ax.twiny()]
#fig1.subplots_adjust(bottom=0.30)
#
#axes[-1].spines['top'].set_position(('axes', 1.25))
#axes[-1].set_frame_on(True)
#axes[-1].patch.set_visible(False)
#axes[2].spines['bottom'].set_position(('axes', -0.25))
#axes[2].set_frame_on(True)
#axes[2].patch.set_visible(False)
#axes[0].plot(imval[0:5], ssa_ni[0:5] - ssa_ni[2], color = colorid[0], marker='s')
#axes[2].plot(aodval, ssa_aod - ssa_aod[3], color = colorid[1], marker='s')
#axes[1].plot(alhval, ssa_alh - ssa_alh[4], color = colorid[3], marker='s')
#axes[3].plot(altval, ssa_alt - ssa_alt[3], color = colorid[5], marker='s')
#axes[0].axhline(0,color = 'k',linestyle = '--')
#axes[0].set_xlabel(r'$n_{i} @ 354 nm$',color = colorid[0])
#axes[0].set_xticks(imval[0:5])
#axes[0].set_xlim(imval[0],imval[4])
#axes[0].set_ylabel('$\omega_0$ @ 550 nm')
#axes[0].set_yticks(np.arange(-0.3,0.4,0.1))
#axes[0].set_ylim(-0.3,0.31)
#axes[2].set_xlabel(r'$\tau_{aer} @ 550 nm$',color = colorid[1])
#axes[2].set_xticks(aodval)
#axes[2].set_xlim(aodval[0],aodval[-1])
#axes[1].set_xlabel(r'$z_{aer}$ [km]',color = colorid[3])
#axes[1].set_xticks(alhval)
#axes[1].set_xlim(alhval[0],alhval[-1])
#axes[-1].set_xlabel(r'$\Delta z_{aer} [km]$',color = colorid[5])
#axes[-1].set_xticks(altval)
#axes[-1].set_xlim(altval[0],altval[-1])
#axes[2].xaxis.set_ticks_position('bottom')
#axes[2].xaxis.set_label_position('bottom')
#
#axes[0].tick_params(axis = 'x',colors=colorid[0])
#axes[2].tick_params(axis = 'x',colors=colorid[1])
#axes[1].tick_params(axis = 'x',colors=colorid[3])
#axes[-1].tick_params(axis = 'x',colors=colorid[5])
#plt.savefig(figdir + 'AAI_sens_mie_uncertainty.png', dpi = 300, transparent = True)
