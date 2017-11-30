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

plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':14})
sns.set_style('white')
colorid = sns.color_palette("hls", 8)
#colorid = [colorid[0],colorid[1],colorid[3],colorid[5]]


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

asyval = np.arange(0.5,1.0,0.1)                # asymmetry factor
asyval = [0.5]

rgval = np.arange(0.1,0.21,0.02)
rgval = [0.14] 

reval = np.arange(1.300, 1.501, 0.05)
reval = [1.4]

imval = np.arange(0.01, 0.06, 0.01)
imval = [0.06] 

szaval = np.arange(0.,76.,15)
szaval = [30.]

saaval = np.arange(0.,181.,45)
saaval = [0.]

vzaval = np.arange(0.,76.,15)
vzaval = [0.]

vaaval = np.arange(0.,181.,45)
vaaval = [180.]

aodval = np.arange(0.5,2.1,0.5)
#aodval = [1.]

alhval = np.arange(2.5,8.6,2)
#alhval = [4.5]

altval = np.arange(0.5,2.1,0.5)
altval = [1.]

#dpval = np.arange(8,25,4)    

Robs = [] 
Rray = []
Rref = []

"""

NOTICE: put para2 (x-axis) loop should be after para1 (legend) loop. The aai list is reshapred as [para1, para2], m-by-n
para1(1): para2(1), para2(2), ... , para2(n)
para1(2): para2(1), para2(2), ... , para2(n)
  :     :    :    ,    :    ,  :  ,      : 
para1(m): para2(1), para2(2), ... , para2(n)

""" 

fig1 = plt.figure(figsize=(6,3))    
ax1 = fig1.add_axes([0.15,0.225,0.5,0.7])

fig2 = plt.figure(figsize=(6,3))    
fig2.subplots_adjust(right=0.75)
ax2 = fig2.add_axes([0.15,0.225,0.5,0.7])
axes = [ax2, ax2.twinx(), ax2.twinx()]
axes[-1].spines['right'].set_position(('axes', 1.3))



for scheme in ['Mie']: 
    
    expname = 'AAI_sens_'+scheme+'_mphy_4' 
    for iAng in Angval:
        for issa in ssaval: 
            for iasy in asyval: 
                for irg in rgval: 
                    for iim in imval:
                        for ire in reval:
                            for isza in szaval:
                                for isaa in saaval:
                                    for ivza in vzaval: 
                                        for ivaa in vaaval: 
                                            for ialt in altval: 
                                                for ialh in alhval: 
                                                    for iaod in aodval:
                                                    
    
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
                                                                            %(irg, ire,iim,isza,isaa,ivza,ivaa,iaod,ialh,ialt,20)
                                                            disamarhdf = tables.open_file(outputdir + expname + src,'r')
                                                            
                                                            disamarout = disamarhdf.root                    
                                                            aairetr = disamarout.parameters.AAI[0]
                                                            Robs.append(disamarout.sun_normalized_radiance.measured_spectrum[0])
                                                            Rray.append(disamarout.sun_normalized_radiance.fitted_spectrum[0])
                                                            Rref.append(disamarout.sun_normalized_radiance.fitted_spectrum[1])
                                                            aaiMie.append(aairetr)
                                                            disamarhdf.close()
                                            
        

    
    # choose 2 parameters 
    para1 = alhval        
    para2 = aodval
    
    # plot sensitivity 
    for i in range(len(para1)):
        
        if scheme == 'HG':
            aaiHGre = np.array(aaiHG).reshape(len(para1),len(para2))  
            ax1.plot(para2, aaiHGre[i,:] / aaiHGre.max(), '-', color = colorid[i], marker = 's', linewidth = 1, label = r'$\alpha$ = '+str(para1[i]))
            
        if scheme == 'Mie':
            aaiMiere = np.array(aaiMie).reshape(len(para1),len(para2))    
            Robs = np.array(Robs).reshape(len(para1),len(para2))    
            Rray = np.array(Rray).reshape(len(para1),len(para2))    
            ax1.plot(para2, aaiMiere[i,:] / aaiMiere.max(), '-', color = colorid[i], marker = 's', linewidth = 1, label = r'z$_{aer}$ = '+str(para1[i])+' [km]')
#            axes[0].plot(para2, (Rray[i,:] - Robs[i,:]) / (Rray - Robs).max(), '--', color = colorid[i], marker = 's', linewidth = 1, label = r'z$_{aer}$ = '+str(para1[i])+' [km]')
#            axes[1].plot(para2, Robs[i,:] / Robs.max(), '-', color = colorid[i], marker = 's', linewidth = 1, label = r'z$_{aer}$ = '+str(para1[i])+' [km]')

    ax1.set_xlim([(para2).min(),(para2).max()])
    ax1.set_xticks(para2)
    ax1.set_yticks(np.arange(0.,1.1,0.2))
    ax1.set_xlabel(r'$\tau_{aer}$')
    ax1.set_ylabel('Normalized AAI')
    ax1.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    
    i = 1 
    axes[0].plot(para2, aaiMiere[i,:] / aaiMiere[i,:].max(), '-', color = 'k', marker = 's', linewidth = 1)
    axes[1].plot(para2, (Rray[i,:] - Robs[i,:]) / (Rray[i,:] - Robs[i,:]).max(), '-', color = colorid[5], marker = 's', linewidth = 1)
    axes[2].plot(para2, Robs[i,:] / Robs[i,:].max(), '-', color = colorid[0], marker = 's', linewidth = 1)    
    ax2.set_xlim([(para2).min(),(para2).max()])
    ax2.set_xticks(para2)
    axes[1].set_yticks(np.arange(0.4,1.01,0.2))
    axes[2].set_yticks(np.arange(0.7,1.01,0.1))
    ax2.set_xlabel(r'$\tau_{aer}$')
    axes[0].set_ylabel(r'Normalized AAI')
    axes[1].set_ylabel(r'Normalized $\Delta I_{\lambda1}$', color = colorid[5])
    axes[2].set_ylabel(r'Normalized $I^{obs}_{\lambda1}$', color = colorid[0])
    axes[1].tick_params(axis = 'y',colors=colorid[5])
    axes[2].tick_params(axis = 'y',colors= colorid[0])
    plt.title(r'z$_{aer}$ = %1.1f [km] $\Delta$z = %1.1f [km]' % (4.5,1))
#    ax2.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    plt.savefig(figdir + 'AAI_sens_aod_analysis.png', dpi = 300, transparent = True)
#plt.close('all')  