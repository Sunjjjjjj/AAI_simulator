"""
Display DISAMAR micro-physical parameter in form of optical parameters

size distribution + refractive index ==> ssa, g 

@author: sunj
"""


######### load python packages
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from AERONET import AERONET_day, AERONET_dailymean
from DSM_Mie_func import create_Expcoef, mphyMie


# figure initialization =======================================================
sns.set_style('white')
plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,\
                    'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11})
colorid = sns.color_palette("hls", 8)
markers = ['s']*3

# directory ===================================================================
casename = 'chile201701/'

maindir = '/usr/people/sunj/Documents/DISAMAR/disamar'
exedir = maindir + '/createExpCoefFiles'
expindir = exedir + '/inputFiles/'
expoutdir = exedir + '/expCoefFiles/'
scaoutdir = exedir +'/scaMatFiles/'
expcoefdir = maindir + '/expCoefFiles/'

figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir += casename

# =============================================================================
years = range(2017,2018)
months = range(1,2)
days = range(26,31)
days.remove(28)
#days = [26]



ROIaer = [50, 65, 5, 45]
# AAI IQR
#factors = [2.8,2,2.2,3]
factors = [1.9,3,2.3,2.3]
# AOD-AAI IQR 
#factors = [2.9,3.2,1.3,2.3]

#factors = [3]

wvl = np.arange(340,677,2)
ssamea = [] 
ssaest = []
for iyear in years:
    for imon in months:
        for dd, iday in enumerate(days):
            
            ssa_AERONET = AERONET_dailymean(iyear, imon, iday, 'ssa')
#            aerosite = AERONET_day(iyear, imon, iday, ROIaer)
#            ssa_AERONET = aerosite[:,8].mean()
            print 'AERONET SSA: ', ssa_AERONET
            ssamea.append(ssa_AERONET)
            
            fine, coarse, bimodal = create_Expcoef(iyear, imon, iday, factors[dd])   
            
            expfilename = expcoefdir + bimodal+'.dat'
            with open(expfilename,'r') as f: 
                content = f.readlines()

            ssamie = [] 
            for ig, il in enumerate(content): 
                if il.find('ssa') != -1:
                    ssamie.append(float(il[5:11]))
            ssamie = np.array(ssamie)
            
            print 'DISAMAR SSA: ', ssamie[np.where(wvl==550)[0][0]]
            ssaest.append(ssamie[np.where(wvl==550)[0][0]])
            
            print 'Relative difference: ', (ssamie[np.where(wvl==550)[0][0]] - ssa_AERONET) / ssa_AERONET
            
fig = plt.figure(figsize = (4.5,3))
ax1 = fig.add_axes([0.15,0.2,0.7,0.7])
plt.plot(days, ssamea,'s-', color = colorid[0], label = 'AERONET')
plt.axhline(np.mean(ssamea), linestyle = '--', color = colorid[0])
plt.plot(days, ssaest,'s-', color = colorid[5], label = 'Retrieved')
plt.axhline(np.mean(ssaest), linestyle ='--', color = colorid[5])
plt.xlabel('Date')
plt.ylabel(r'$\omega_0$')
plt.xticks(days)
plt.yticks(np.arange(0.75,0.96,0.05))
plt.ylim(0.75,0.96)
plt.legend(frameon=False,loc=4)
plt.savefig(figdir + 'SSA_AERONET_DSM.png', dpi = 300, transparent = True)

######## DISAMAR micro-physics
#plt.figure()
#temp1 = []
#temp2 = []
#temp3 = [] 
#rgval = np.arange(0.1,0.41,0.05)
#for irg in rgval: 
#    expfilename = 'AAI_sim_Mie_bimodal' % (ialh, ialt,ifactor)
#    wvl, ssamie, ang, scamat = mphyMie(expfilename)
#    temp1.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]] / scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
#    temp2.append(scamat[np.where(scamat[:,0] == 354)[0],np.where(ang==150)[0][0]])
#    temp3.append(scamat[np.where(scamat[:,0] == 388)[0],np.where(ang==150)[0][0]]) 
#
#plt.plot(temp1,'s')
#plt.ylim(0.5,1.5)
#plt.plot(temp2)
#plt.plot(temp3)

plt.close('all')

