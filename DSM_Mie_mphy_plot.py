# -*- coding: utf-8 -*-
"""
Intercomparison among satellite data

@author: sunj
"""


######### load python packages
import sys, os
import shutil
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 
import glob
from scipy import ndimage
from mpl_toolkits.basemap import Basemap
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn.metrics import mean_squared_error
from scipy.optimize import leastsq
from shutil import copyfile
import subprocess 
from scipy import optimize
from math import sqrt
from OMI import gridOMAERO, OMAERO4cfg
from MODIS import gridMODIS, MODIS4cfg
from GOME2 import gridGOME2




# directory 
exedir = '/usr/people/sunj/Documents/DISAMAR/disamar/createExpCoefFiles/'
expindir = exedir + 'inputFiles/'
expoutdir = exedir + 'expCoefFiles/'
scaoutdir = exedir +'scaMatFiles/'
expcoefdir = '/usr/people/sunj/Documents/DISAMAR/disamar/expCoefFiles/'


maindir = '/usr/people/sunj/Documents/DISAMAR'
outputdir = maindir + '/pydisamar/output/'
figdir = '/usr/people/sunj/Dropbox/EGU2017/Figure/'


plt.rcParams.update({'xtick.color': 'k','ytick.color':'k','text.color':'k','axes.labelcolor':'k', 'lines.linewidth':2,'font.size':14,'axes.labelsize':14,'axes.titlesize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':14})

######## sensitivity study result 
rgval = np.array([0.14])
sigmaval = np.array([1.28,1.38,1.48,1.58,1.68])

#imval = [ 0.008, 0.01, 0.012, 0.014, 0.016]

#rgval = [0.06]
c0 = 3 
rg0 = 0.16

sns.set_style('white')
fig1 = plt.figure(1,figsize=(9.2,2.5))    


for irg in rgval: 

    expfilename = 'AAI_sens_Mie_geom'
    expfile = expoutdir + expfilename + '.dat'
    with open(expfile,'r') as f: 
        content = f.readlines()
    
    ssamie = [] 
    for il in content: 
        if il.find('ssa') != -1:
            ssamie.append(float(il[5:11]))
    wvl = np.arange(340,441,2)
    ssamie = np.array(ssamie)

    cmap = sns.cubehelix_palette(9, start=.5, rot=.3)
    cmap = cmap[1:]
    sns.set_palette(cmap)

    ax1 = fig1.add_axes([0.1,0.2,0.4,0.7])
    ax1.plot(wvl,ssamie)
    plt.ylabel('Single scattering albedo [-]')
    plt.xlabel('Wavelength [nm]')    
    plt.xlim([340,440])
    plt.ylim([0.8,1])
    plt.legend(loc =1, frameon = False)


    
    scafile = scaoutdir + expfilename + '.dat'
#    scafile = 'strong_abs_320-400.dat'
    with open(scafile,'r') as f: 
        content = f.readlines()
        
    scamat = []
    for il in content: 
        if il.find('F11') == 0:
            scamat.append(il[5:-1])
            
    scamat = np.loadtxt(scamat)
    
    ang = np.arange(0,181)
    i354 = np.where(scamat[:,0] == 354)[0]
    i388 = np.where(scamat[:,0] == 388)[0]
#    i550 = np.where(scamat[:,0] == 550)[0]  
    ax2 = fig1.add_axes([0.575,0.2,0.4,0.7])
    ax2.plot(ang,scamat[i354,1:][0], 'k-', label = 'Mie' )
#    ax2.plot(ang,scamat[i388,1:][0], 'k--',label = '388 nm' )
    plt.xlabel('Scattering angle [deg]')
    plt.ylabel('Phase function')    
    plt.xticks([0,30,60,90,120,150,180])
    plt.yticks([0,10,20,30,40,50])
    ax2.set_yscale('log')
#    plt.title('Phase function')

    
    g = 0.7
    P_Ray = (1-g**2) / (1+g**2-2*g*np.cos(ang/180.*np.pi))**(3/2)
    plt.plot(ang, P_Ray, 'r',label = 'HG')    
    
    plt.legend(loc =1, frameon = False)    
#    plt.savefig(figdir + 'DSM_Mie_rg_ssa.png', dpi = 300)    
#    plt.savefig('/usr/people/sunj/Dropbox/EGU2017/Figures/DSM_Mie_rg_ssa.png', dpi = 300, transparent = True)    
    

            

fig2 = plt.figure(figsize=(5,5))

ax = plt.subplot(111, projection = 'polar')    

for irg in rgval: 

    expfilename = 'AAI_sens_Mie_geom'
    expfile = expoutdir + expfilename + '.dat'
    with open(expfile,'r') as f: 
        content = f.readlines()
    
    
    scafile = scaoutdir + expfilename + '.dat'
#    scafile = 'strong_abs_320-400.dat'
    with open(scafile,'r') as f: 
        content = f.readlines()
        
    scamat = []
    for il in content: 
        if il.find('F11') == 0:
            scamat.append(il[5:-1])
            
    scamat = np.loadtxt(scamat)
    
    ang = np.arange(0,181)
    i354 = np.where(scamat[:,0] == 354)[0]
    i388 = np.where(scamat[:,0] == 388)[0]
#    i550 = np.where(scamat[:,0] == 550)[0]  
#    ir354 = i354[]

    ax.plot(ang/180. * np.pi, scamat[i354,1:][0], 'k-', linewidth = 2, label = 'Mie' )
    ax.plot((ang+180.)/180. * np.pi, scamat[i354,1:][0][::-1], 'k-')
##    ax2.plot(ang,scamat[i388,1:][0], 'k--',label = '388 nm' )
#    plt.xlabel('Scattering angle [deg]')
#    plt.ylabel('Phase function')    
##    plt.xticks([0,30,60,90,120,150,180])
##    plt.yticks([0,10,20,30,40,50])
#    ax.set_yscale('log')
#    plt.title('Phase function')

    g = 0.7
    ang = np.arange(0,361) 
    P_Ray = (1-g**2) / (1+g**2-2*g*np.cos(ang/180.*np.pi))**(3/2)
    plt.plot( ang / 180.*np.pi,  P_Ray, 'r', linewidth = 2, label = 'HG')  
    ax.set_rticks([0,1,10])
    ax.set_rscale('log')
    ax.grid(True)
    

    
    plt.legend(loc = 6, frameon = False)    
    plt.savefig(figdir + 'phase_function.png', dpi = 300)    
#    plt.savefig('/usr/people/sunj/Dropbox/EGU2017/Figures/DSM_Mie_rg_ssa.png', dpi = 300, transparent = True)    

#r = np.arange(0, 2, 0.01)
#theta = 2 * np.pi * r
#
#ax = plt.subplot(111, projection='polar')
#ax.plot(theta, r)



#fig2 = plt.figure(2,figsize=(9.2,2.5))    
#for isigma in sigmaval: 
#
##    expfilename = 'Mongu_sens_sigma_%1.4f' % (isigma)
#    expfile = expoutdir + expfilename + '.dat'
#    with open(expfile,'r') as f: 
#        content = f.readlines()
#    
#    ssamie = [] 
#    for il in content: 
#        if il.find('ssa') != -1:
#            ssamie.append(float(il[5:11]))
#    wvl = np.arange(340,441,2)
#    ssamie = np.array(ssamie)
#
#    cmap = sns.cubehelix_palette(9, start=.5, rot=.3)
#    cmap = cmap[1:]
#    sns.set_palette(cmap)
#
#    ax1 = fig2.add_axes([0.1,0.2,0.4,0.7])
#    ax1.plot(wvl,ssamie)
#    plt.ylabel('Single scattering albedo [-]')
#    plt.xlabel('Wavelength [nm]')    
#    plt.xlim([340,550])
#    plt.ylim([0.8,1])
#    plt.legend(loc =1, frameon = False)
#
#
#    
#    scafile = scaoutdir + expfilename + '.dat'
#    with open(scafile,'r') as f: 
#        content = f.readlines()
#        
#    scamat = []
#    for il in content: 
#        if il.find('F11') == 0:
#            scamat.append(il[5:-1])
#            
#    scamat = np.loadtxt(scamat)
#    
#    ang = np.arange(0,181)
#    i340 = np.where(scamat[:,0] == 340)[0]
#    i380 = np.where(scamat[:,0] == 380)[0]
##    i550 = np.where(scamat[:,0] == 550)[0]
#    ax2 = fig2.add_axes([0.575,0.2,0.4,0.7])
#    ax2.plot(ang,scamat[i380,1:][0], label = 'vg = %1.2f' %isigma )
#    plt.xlabel('Scattering angle [deg]')
#    plt.ylabel('Phase function')    
#    plt.xticks([0,30,60,90,120,150,180])
#    plt.yticks([0,10,20,30,40,50])
##    plt.title('Phase function')
#    plt.legend(loc =1, frameon = False)
##    plt.savefig(figdir + 'DSM_Mie_sigma_ssa.png', dpi = 300)    
##    plt.savefig('/usr/people/sunj/Dropbox/EGU2017/Figures/DSM_Mie_sigma_ssa.png', dpi = 300, transparent = True)    