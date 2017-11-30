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
colorid = [colorid[0],colorid[1],colorid[3],colorid[5]]


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

"""

NOTICE: put para2 (x-axis) loop should be after para1 (legend) loop. The aai list is reshapred as [para1, para2], m-by-n
para1(1): para2(1), para2(2), ... , para2(n)
para1(2): para2(1), para2(2), ... , para2(n)
  :     :    :    ,    :    ,  :  ,      : 
para1(m): para2(1), para2(2), ... , para2(n)

""" 

fig1 = plt.figure(figsize=(6,3))    
ax1 = fig1.add_axes([0.15,0.225,0.5,0.7])

for scheme in ['Mie']: 
    
    expname = 'AAI_sens_'+scheme+'_uncertainty' 
    for iAng in Angval:
        for issa in ssaval: 
            for iasy in asyval: 
                for irg in rgval: 
                    for ire in reval:
                        for isza in szaval:
                            for isaa in saaval:
                                for ivza in vzaval: 
                                    for ivaa in vaaval: 
                                        for ialt in altval: 
                                            for ialh in alhval: 
                                                for iaod in aodval:
                                                    for iim in imval:
    
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
                                                            aaiMie.append(aairetr)
                                                            disamarhdf.close()
                                            
        

    
    # choose 2 parameters 
    para1 = aodval        
    para2 = imval
    
    # plot sensitivity 
    for i in range(len(para1)):
        
        if scheme == 'Mie':
            aaiMiere = np.array(aaiMie).reshape(len(para1),len(para2))    
            ax1.plot(para2, aaiMiere[i,:] , '-', color = colorid[i], marker = 's', linewidth = 1, label = r'$\tau_{aer}$ = '+str(para1[i]))

    plt.xlim([(para2).min(),(para2).max()])
    plt.xticks(para2)
    plt.yticks(np.arange(0.,1.1,0.2))
    plt.xlabel('ni')
    plt.ylabel('AAI')
    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
#    plt.savefig(figdir + 'AAI_sens_uncertainty.png', dpi = 300, transparent = True)
#plt.close('all')  