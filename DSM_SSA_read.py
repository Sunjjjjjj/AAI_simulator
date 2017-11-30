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
from AERONET import AERONET_dailymean
from DSM_Mie_func import create_Expcoef, mphyMie
from scipy import stats


sns.set_style('white')
plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':1,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11})
colorid = sns.color_palette("hls", 8)
markers = ['s','s','s']

# directory 
maindir = '/usr/people/sunj/Documents/DISAMAR/disamar'
exedir = maindir + '/createExpCoefFiles'
expindir = exedir + '/inputFiles/'
expoutdir = exedir + '/expCoefFiles/'
scaoutdir = exedir +'/scaMatFiles/'
expcoefdir = maindir + '/expCoefFiles/'
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'

casename = 'AAI_sensitivity/'
if not os.path.exists(figdir + casename):
    os.makedirs(figdir + casename)   
figdir = '/usr/people/sunj/Dropbox/Paper_Figure/'+ casename




reval = np.arange(1.300, 1.501, 0.05)
reval = [1.4]

imval = np.arange(0.02, 0.0401, 0.005) 
imval = [0.0546, 0.0303, 0.0171]
#imval = [0.04424488,  0.03,  0.01720842]

ssa = [] 

wvl = np.arange(340,677,2)



for iim in imval:
    expfilename = expcoefdir + 'rg-0.14_nr-1.4000_ni-%1.4f.dat' % (iim)
    with open(expfilename,'r') as f: 
        content = f.readlines()

    ssamie = [] 
    for ig, il in enumerate(content): 
        if il.find('ssa') != -1:
            ssamie.append(float(il[5:11]))
    ssamie = np.array(ssamie)
    
    print 'ni = ', iim, '  DISAMAR SSA: ', ssamie[np.where(wvl==550)[0][0]]
    ssa.append(ssamie[np.where(wvl==550)[0][0]])
        
plt.figure()
plt.plot(imval, ssa,'s-')            

slope1, intercept1, r_value1, p_value, std_err = stats.linregress(imval,ssa)
print 'k = ', slope1
print (ssa - ssa[1])
print (ssa - ssa[1]) / ssa[1] * 1e2, '%'









#reval = np.arange(1.300, 1.501, 0.05)
##reval = [1.4]
#
#imval = np.arange(0.02, 0.0401, 0.005)
#imval = [0.03]
#
#ssa = [] 
#
#for ire in reval:
#    expfilename = expcoefdir + 'rg-0.14_nr-%1.4f_ni-0.0300.dat' % (ire)
#    with open(expfilename,'r') as f: 
#        content = f.readlines()
#
#    ssamie = [] 
#    for ig, il in enumerate(content): 
#        if il.find('ssa') != -1:
#            ssamie.append(float(il[5:11]))
#    ssamie = np.array(ssamie)
#    
#    print 'ni = ', iim, '  DISAMAR SSA: ', ssamie[np.where(wvl==550)[0][0]]
#    ssa.append(ssamie[np.where(wvl==550)[0][0]])
#        
#plt.figure()
#plt.plot(reval, ssa,'s-')            
#
#slope1, intercept1, r_value1, p_value, std_err = stats.linregress(reval,ssa)
#print 'k = ', slope1
#print (ssa - ssa[2]) / ssa[2] *1e2, '%'



