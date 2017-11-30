# -*- coding: utf-8 -*-
"""
Run AAI simulation Mie scheme

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
#import tables
#from scipy import ndimage
from math import sqrt
from DISAMAR import rundisamar
from AERONET import AERONET_day, AERONET_dailymean
from scipy import spatial 
from timezone import timezone



maindir = '/usr/people/sunj/Documents/DISAMAR/'
outputdir = maindir + 'pydisamar/output/'
nobackupdir = '/nobackup/users/sunj/'

exedir = maindir + 'disamar/createExpCoefFiles/'
expindir = exedir + 'inputFiles/'
expoutdir = exedir + 'expCoefFiles/'
scaoutdir = exedir +'/scaMatFiles/'
os.chdir(exedir)

expcoefdir = maindir + 'disamar/expCoefFiles/'
template = expindir + 'smoke_template.in'



#def create_Expcoef(year, month, day, aerosite, factor):
def create_Expcoef(year, month, day, factor):

    sigma_f = 1.5
    sigma_c = 2.0

#    refr354, refi354, refr388, refi388, refr550, refi550, ssa550, r_f, r_c, wfc_n, wfc_v = aerosite[:,2:].mean(axis=0)
    
    r_f, r_c, wfc_n, wfc_v = AERONET_dailymean(year,month,day,'sizfunc')
    refr354, refi354, refr388, refi388, refr550, refi550 = AERONET_dailymean(year,month,day,'refidx')       
    
#    refr354 -= 0.0
#    refr388 -= 0.0
#    refr550 -= 0.0
#    print refi354, refi388,refi550
    refi354 *= factor
    refi388 *= factor
    refi550 *= factor
#    print refi354, refi388,refi550
    
    t0 = time.time() 
    print ('**** creating fine mode') 
    imode = 'AAI_sim_Mie_fine_%4i%02i%02i' %(year, month, day)
    finemode = imode
    with open(template,'r') as f: 
        content = f.read()
#                        
    content = content.replace('XnameX',imode)
    content = content.replace('XrgX','%1.2f' %(r_f))
    content = content.replace('XsigmaX','%1.2f' %(sigma_f))   
    content = content.replace('Xre354X','%1.4f' %(refr354))
    content = content.replace('Xim354X','%1.4f' %(refi354))
    content = content.replace('Xre388X','%1.4f' %(refr388))
    content = content.replace('Xim388X','%1.4f' %(refi388))        
    content = content.replace('Xre550X','%1.4f' %(refr550))
    content = content.replace('Xim550X','%1.4f' %(refi550))  
    
    with open(expindir + imode+ '.in','w') as f:
        f.write(content)                        
    
    if os.path.isfile(exedir + 'Config.in'):
        os.remove(exedir + 'Config.in')
    os.symlink(expindir + imode + '.in', exedir + 'Config.in')
    os.system(exedir + "createExpCoefFile.exe") 
    
    shutil.copy2(expoutdir + imode+'.dat', expcoefdir + imode + '.dat')
    
    t1 = time.time() 
    print 'time for fine mode:', t1 - t0 
    
    print ('**** creating coarse mode') 
    imode = 'AAI_sim_Mie_coarse_%4i%02i%02i' %(year, month, day) 
    coarsemode = imode
    with open(template,'r') as f: 
        content = f.read()
#    
    content = content.replace('XnameX',imode)
    content = content.replace('XrgX','%1.2f' %(r_c))
    content = content.replace('XsigmaX','%1.2f' %(sigma_c))   
    content = content.replace('Xre354X','%1.4f' %(refr354))
    content = content.replace('Xim354X','%1.4f' %(refi354))
    content = content.replace('Xre388X','%1.4f' %(refr388))
    content = content.replace('Xim388X','%1.4f' %(refi388))        
    content = content.replace('Xre550X','%1.4f' %(refr550))
    content = content.replace('Xim550X','%1.4f' %(refi550))  
    
    with open(expindir + imode+ '.in','w') as f:
        f.write(content)                        
    
    if os.path.isfile(exedir + 'Config.in'):
        os.remove(exedir + 'Config.in')
    os.symlink(expindir + imode + '.in', exedir + 'Config.in')
    os.system(exedir + "createExpCoefFile.exe") 
    
    shutil.copy2(expoutdir + imode+'.dat', expcoefdir + imode + '.dat')
    t2 = time.time() 
    print 'time for coarse mode:', t2 - t1 
    
    print('**** creating bimode') 
    bimodal = 'AAI_sim_Mie_bimodal_%4i%02i%02i' %(year, month, day) 
    filename_out = combine_Expcoef(expcoefdir + finemode + '.dat', expcoefdir + coarsemode + '.dat', wfc_n)
    os.rename(filename_out, expcoefdir + bimodal + '.dat')
    
    t3 = time.time()
    print 'time for bimodal:', t3 - t2 
    
    return finemode, coarsemode, bimodal
    
    
    
    
    
def combine_Expcoef(filename_1,filename_2, wfc):  
    
    
    w1 = wfc
    w2 = 1. - wfc
    
    
    filename_out =  expcoefdir + 'AAI_sim_Mie_bimodal.dat'
    
    fo = '' # output string
    
    
    #-----------------------------------------------------------------------------------------
    # Read the two files into lists f1 and f2
     
    f = open(filename_1, 'r')
    f1 = (f.read()).split('\n')
    f.close()
    
    f = open(filename_2, 'r')
    f2 = (f.read()).split('\n')
    f.close()
    
    
    #-----------------------------------------------------------------------------------------
    # combine the two headers
    
    nheader_1 = int( (f1[0]).split('=')[0] )
    nheader_2 = int( (f2[0]).split('=')[0] )
    
    n_header = nheader_1 + nheader_2 + 10
    
    fo = fo + '%4i  = number of header lines after this line\n' % n_header
    
    fo = fo + '\nCombination of %s and %s\n' % (filename_1, filename_2)
    fo = fo + 'weight file_1: %f, weight file_2: %f\n' % (w1, w2)
    fo = fo + '---------------------------------------------------------------------------\n'
    
    fo = fo + '\nHeader file %s\n' % filename_1
    fo = fo + '---------------------------------------------------------------------------\n'
    
    for i in range(1, nheader_1 + 1):
        fo = fo + f1[i] + '\n'
    
    fo = fo + '\nHeader file %s\n' % filename_2
    fo = fo + '---------------------------------------------------------------------------\n'
    
    for i in range(1, nheader_2 + 1):
        fo = fo + f2[i] + '\n'
    
    #-----------------------------------------------------------------------------------------
    # read and write the general info
    
    gnrl_1 = list()
    gnrl_2 = list()
    gnrl_labels = list()
    
    for i in range(6):
        gnrl_1.append( float( (f1[nheader_1+1+i]).split()[0] ) )
    
    for i in range(6):
        gnrl_2.append( float( (f2[nheader_2+1+i]).split()[0] ) )
        
    
    for i in range(6):
        gnrl_labels.append( ' '.join( (f1[nheader_1+1+i]).split()[1:] ) )
    
    # check if the number of expansion coefficients and wavelengths match
    if ( int(gnrl_1[0]) != int(gnrl_2[0]) ):
        raise ValueError('%s do not match: %i, %i' % ( gnrl_labels[0], int(gnrl_1[0]), int(gnrl_2[0]) ) )
        
    if ( int(gnrl_1[1]) != int(gnrl_2[1]) ):
        raise ValueError('%s do not match: %i, %i' % ( gnrl_labels[1], int(gnrl_1[1]), int(gnrl_2[1]) ) )
    
    ncoef  = int(gnrl_1[0])
    nwavel = int(gnrl_1[1])
    
    # write the general info
    fo = fo + '%i %s\n' % ( int(gnrl_1[0]), gnrl_labels[0])
    fo = fo + '%i %s\n' % ( int(gnrl_1[1]), gnrl_labels[1])
    
    for i in range(2,6):
        fo = fo + '%16.10e %s\n' % ( w1*gnrl_1[i] + w2*gnrl_2[i], gnrl_labels[i])
        
    #-----------------------------------------------------------------------------------------
    # wavelength loop
    for iwavel in range(nwavel):
        o1 = 1 + nheader_1 + len(gnrl_labels) + iwavel * ( 6 + ncoef + 1)
        o2 = 1 + nheader_2 + len(gnrl_labels) + iwavel * ( 6 + ncoef + 1)
    
        wavel_1 = float( (f1[o1+1]).split()[0] )
        ncext_1 = float( (f1[o1+2]).split()[0] )
        ssa_1   = float( (f1[o1+3]).split()[0] )
    
        wavel_2 = float( (f2[o2+1]).split()[0] )
        ncext_2 = float( (f2[o2+2]).split()[0] )
        ssa_2   = float( (f2[o2+3]).split()[0] )
        
        # check if the wavelengths match
        if ( abs(wavel_1-wavel_2) > 1e-3 ):
            raise ValueError('Wavelengths do not match: %f, %f' % ( wavel_1, wavel_2 ) )
    
        ew1 = (w1 * ncext_1) / (w1 * ncext_1 + w2 * ncext_2 ) # extinction weigths
        ew2 = 1. - ew1
            
        sw1 = (w1 * ssa_1 * ncext_1) / (w1 * ssa_1 * ncext_1 + w2 * ssa_2 * ncext_2 ) # scattering weights
        sw2 = 1. - sw1
    
        # write the wavelength header
        fo = fo + '\n'
        fo = fo + '%14.8f   wavelength\n' % (wavel_1)
        fo = fo + '%14.8f   Cext/Cext550\n' % (w1*ncext_1 + w2*ncext_2)
        fo = fo + '%14.8f   ssa\n' % (ew1*ssa_1 + ew2*ssa_2)
        fo = fo + '\n'
        fo = fo + f1[o1+5]+ '\n'
        
        # compute and write the coefs
        for icoef in range(ncoef+1):
        
            fo = fo + '%4i' % icoef
        
            scoefs_1 = (f1[o1+6+icoef]).split()
            scoefs_2 = (f2[o2+6+icoef]).split()
            
            for i in range(1, len(scoefs_1) ):
                val = sw1 * float(scoefs_1[i]) + sw2 * float(scoefs_2[i])
                fo = fo + '%15.6e' % val     # sunj ' '*3 
    
            fo = fo + '\n'
    
    #-----------------------------------------------------------------------------------------
    # write the output file
    f = open(filename_out, 'w')
    f.write(fo)
    f.close()
    
    print('DSM_Mie_combine_Expcoef Done!')
    return filename_out     
    






######## DISAMAR micro-physics

def mphyMie(expfilename): 
    
    # HG SSA
    wvl = np.arange(340,677,2)
    ssaHG = 0.75 * np.ones(len(wvl))
    
    # HG phase function
    ang = np.arange(0,181)
    g = 0.8
    phafuncHG = (1-g**2) / (1+g**2-2*g*np.cos(ang/180.*np.pi))**(3/2.) / (4*np.pi)
    
        
     
    g550 = [] 
    ssa = [] 
    
    
    expfile = expcoefdir + expfilename + '.dat'
    print expfile
    with open(expfile,'r') as f: 
        content = f.readlines()
    
    gmie = [] 
    for ig, il in enumerate(content): 
        if il.find('alpha1') != -1:
            temp1 = content[ig+2]
            temp2 = float(temp1.split('  ')[2]) / 3.
            gmie.append(temp2)
            
    g550.append(gmie[np.where(wvl == 354)[0][0]])
    
    ssamie = [] 
    for ig, il in enumerate(content): 
        if il.find('ssa') != -1:
            ssamie.append(float(il[5:11]))
    ssamie = np.array(ssamie)
    
#    plt.plot(wvl,ssaHG,color = colorid[0], marker = 's', markevery = 20, linewidth = 1, label = 'HG g@354=%1.2f' % (g))       
#    plt.plot(wvl,ssamie,color = colorid[1], marker = 's', markevery = 40, linewidth = 1, label = 'Mie g@354=%1.2f' % (g550[0]))        
#    plt.ylabel('SSA [-]')
#    plt.xlabel('Wavelength [nm]')
#    plt.xlim([340,675])
#    plt.ylim([0.8,1.0])
#    plt.xticks(np.arange(340,677,50))
#    plt.yticks(np.arange(0.5,0.95,0.1))
#    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
#    plt.savefig(figdir + 'SSA.png' , dpi = 300)    
    
    # plot phase function 
    expfile = expcoefdir + expfilename + '.dat'
#    plt.yscale('log')
    scafile = scaoutdir + expfilename + '.dat'
    
    with open(scafile,'r') as f: 
        content = f.readlines()
        
    scamat = []
    for il in content: 
        if il.find('F11') == 0:
            scamat.append(il[5:-1])            
    scamat = np.loadtxt(scamat)
    
#    i354 = np.where(scamat[:,0] == 354)[0]
#    i388 = np.where(scamat[:,0] == 388)[0]
#    i550 = np.where(scamat[:,0] == 550)[0]
#    
#    fig2 = plt.figure(figsize=(6,3))  
#    ax1 = fig2.add_axes([0.15,0.225,0.5,0.7])            
#    plt.plot(ang,phafuncHG, color = colorid[0], marker = 's', markevery = 20, linewidth = 1, label = 'HG')
#    plt.plot(ang,scamat[i354,1:][0], color = colorid[1], marker = 's', markevery = 20, linewidth = 1, label = 'Mie @354')
#    plt.plot(ang,scamat[i388,1:][0], color = colorid[2], marker = 's', markevery = 20, linewidth = 1, label = 'Mie @388')
#    plt.xlabel('Scattering angle [deg]')
#    plt.ylabel('Phase function [-]')    
#    plt.yscale('log')
#    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
#    plt.savefig(figdir + 'phafunc.png' , dpi = 300)     
#    
#    fig3 = plt.figure(figsize=(6,3))  
#    ax1 = fig3.add_axes([0.15,0.225,0.5,0.7])            
#    plt.plot(ang,scamat[i354,1:][0] / scamat[i388,1:][0] , color = colorid[1], marker = 's', markevery = 20, linewidth = 1, label = 'Mie 354/388')
#    plt.ylim(0.5,1.5)    
#    plt.xlabel('Scattering angle [deg]')
#    plt.ylabel('Phase function [-]')    
#    plt.legend(frameon = False, loc = 6, ncol = 1, bbox_to_anchor=(1,0.5))   
    
    
    
    return wvl, ssamie, ang, scamat   
    
#year = 2017
#month = 1
#day = 27
#    
#create_Expcoef(year, month, day,1)    
