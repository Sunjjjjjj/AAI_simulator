# -*- coding: utf-8 -*-
"""
reading AERONET data

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
from mpl_toolkits.basemap import Basemap
import pandas as pd
from pandas import Series, DataFrame, Panel
from scipy import spatial
from timezone import timezone
from scipy.signal import argrelextrema
from numpy import trapz
from scipy import interpolate 





def AERONET_day(year, month, day, ROI): 
    count = 0 
    aerdir = '/nobackup/users/sunj/AERONET/INV/LEV15/DUBOV/DAILY/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340
    filelist = glob.glob( aerdir + '*.dubovikday')

    vldsite = []
    
    for ff in filelist[:]:
        count += 1 
        
        f = open(ff, 'r')
        lines = f.readlines()
        headerlines = lines[0]
        paraname = lines[3]
        
        aerdata = pd.read_csv(ff, sep=",", header = 3)
        
        
        
        # get refrective index 
        idx1 = paraname.find('REFR')
        wvl1 = int(paraname[idx1+5:idx1+8])
        wvl2 = int(paraname[idx1+15:idx1+18])
        
        refr1 = aerdata['REFR(%3i)' % (wvl1)]
        refr2 = aerdata['REFR(%3i)' % (wvl2)]
        refi1 = aerdata['REFI(%3i)' % (wvl1)]
        refi2 = aerdata['REFI(%3i)' % (wvl2)]
        
        ar =  (refr2 - refr1) / (wvl2 - wvl1)
        ai =  (refi2 - refi1) / (wvl2 - wvl1)
        
        refr354 = refr1 - ar * (wvl1 - 354.)      
        refi354 = refi1 - ai * (wvl1 - 354.)   
        refr388 = refr1 - ar * (wvl1 - 388.)      
        refi388 = refi1 - ai * (wvl1 - 388.)   
        refr550 = refr1 + ar * (550. - wvl1)
        refi550 = refi1 + ai * (550. - wvl1)
    
        refidx = np.column_stack([refr354, refi354, refr388, refi388, refr550, refi550])
        
        
        # SSA @ 550 nm 
        idx1 = paraname.find('SSA')
        wvl1 = int(paraname[idx1+3:idx1+6])
        wvl2 = int(paraname[idx1+12:idx1+15])    
    
        ssa1 = aerdata['SSA%3i-T' % (wvl1)]
        ssa2 = aerdata['SSA%3i-T' % (wvl2)]
        
        ssa550 = (ssa2 - ssa1) / (wvl2 - wvl1) * (550. - wvl1) + ssa1  
        
        
        # size distribution function 
        idx1 = paraname.find('0.050000')
        idx2 = paraname.find(',Inflection_Point')
        
        size = paraname[idx1:idx2].split(',')
        
        
        for i,s in enumerate(size):
            size[i] = float(s)
    
        size = np.array(size)
        sizeint = np.arange(np.min(size),np.max(size),0.01) 
    
        paraname = paraname.split(',')
        prob = np.array(aerdata)[:,62:84] 
        probint = np.zeros([len(prob),len(sizeint)])
    
        r_f = np.ones(len(ssa550)) * np.nan
        r_c = np.ones(len(ssa550)) * np.nan
        wfc_n = np.ones(len(ssa550)) * np.nan
        wfc_v = np.ones(len(ssa550)) * np.nan
    
        # covert volume distribution to number distribution
        # r_n = r_v / exp(3*sigma_v**2)
        # sigma_v = sigma_n = ln(sigma) 
        # sigma is the geometric standard deviation (also input in MMP and DISAMAR), typical range for aerosols is 1.5 to 2. 
        sigma_f = 1.5
        sigma_c = 2.0 
                
        for i in range(len(prob)):
     
            f = interpolate.interp1d( size, np.array(list(prob[i])), kind ='linear')
            probint[i,:] = f (sizeint)
            peaks = argrelextrema(probint[i,:],np.greater) 
    
            
            # only one model 
            if peaks[0].size > 1: 
                
            
                fine = sizeint[peaks[0][0]]
                coarse = sizeint[peaks[0][-1]]
                
                
                fcidx = argrelextrema(probint[i,:],np.less)[0]
                
                if fcidx.size > 0: 
                    fcidx = argrelextrema(probint[i,:],np.less)[0][0]
                    
                    Cv_f = trapz(probint[i,0:fcidx+1],dx = 0.01)
                    Cv_c = trapz(probint[i,fcidx:-1], dx = 0.01)
        
                    
                    r_f[i] = fine / np.exp(3*(np.log(sigma_f)**2))
                    r_c[i] = coarse / np.exp(3*(np.log(sigma_c)**2))
                    
                    
                    rfc_Cn = (Cv_f / Cv_c) * (r_c[i] / r_f[i])**3 * np.exp(-4.5 * (sigma_f**2-sigma_c**2))  
                    
                    rfc_Cv = Cv_f / Cv_c
                    
                    wfc_n[i] = rfc_Cn / (rfc_Cn + 1.) 
                    wfc_v[i] = rfc_Cv / (rfc_Cv + 1.)        
    #            else:
    #                print i 
    
            if peaks[0].size == 1: 
                if sizeint[peaks[0][0]] < 1:
                    fine = sizeint[peaks[0][0]]
                    coarse = 0
    
                    r_f[i] = fine / np.exp(3*(np.log(sigma_f)**2))
                    r_c[i] = coarse / np.exp(3*(np.log(sigma_c)**2))
    
                    wfc_n[i] = 1 
                    wfc_v[i] = 1 
                    
                else:
                    fine = 0 
                    coarse = sizeint[peaks[0][0]]
                    
                    r_f[i] = fine / np.exp(3*(np.log(sigma_f)**2))
                    r_c[i] = coarse / np.exp(3*(np.log(sigma_c)**2))
    
                    wfc_n[i] = 0 
                    wfc_v[i] = 0
                    
                
        sizfunc = np.column_stack([r_f, r_c, wfc_n, wfc_v])    
        
        paradata = np.column_stack([refidx, ssa550, sizfunc])
        
        
        # location filtering 
        idx1 = headerlines.find('long')
        idx2 = headerlines.find(',elev')
        
        headerlines = headerlines[idx1:idx2]
        
        idx1 = headerlines.find('long=')
        idx2 = headerlines.find('lat=')
        
        aerlon = float(headerlines[idx1+5:idx2-1])
        aerlat = float(headerlines[idx2+4:])
        
    
      
        # time filtering 
        if (aerlat > ROI[0]) and (aerlat < ROI[1]) and (aerlon > ROI[2]) and (aerlon < ROI[3]):
            aerdata = pd.read_csv(ff, sep=",", header = 3)
            date = aerdata['Date(dd-mm-yyyy)']
    
    
            if year % 4 == 0:
                dayspermon = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
            else:
                dayspermon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    
            
            for idate in date:
                if (float(idate[-4:]) == year) and (float(idate[3:5]) == month) and (float(idate[0:2]) == day):
                    
                    print count, 'reading %s' % (ff)
                    
                    ts = aerdata['Julian_Day']
                    
                    if month > 1:
                        doy = dayspermon[0:month-1].sum() + day 
                    else:
                        doy = day 
                    
                    idx = np.where(ts == doy)[0][0]
                    print idate
                    
                    sitedata = np.column_stack([aerlat, aerlon, paradata[idx,:].reshape(1,len(paradata[idx,:]))])
                    vldsite.append(sitedata.reshape(sitedata.size))
                    
    #                if para == 'aod' or para == 'AE':
    #                    idx = np.where(ts == doy)[0][0]
    #                else:
    #            #        print abs(ts-doy-13./24)
    #                    idx = np.argmin(abs(ts-doy-13./24))
                
    return np.array(vldsite) 
          

#year = 2017
#month = 8 
#day = 18
#ROI = [60.,70., -65.,-45.] 
#
#count = 0 
#count1 = 0
#
#data = AERONET_day(year, month, day, ROI)
#
#plt.figure() 
#i=8
#plt.plot(data[:,i],'s-')
#plt.plot(range(19),np.ones(19)*data[:,i].mean())
#plt.figure()
#
#plt.plot(data[:,1], data[:,0], 's')
#plt.xlim(-160,-45)
#plt.ylim(40,75)











"""
read AERONET site for Chile201701
temporary script 
""" 

def AERONET_dailymean(year, month, day, para):    
    
    if year % 4 == 0:
        dayspermon = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    else:
        dayspermon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        
    ######## AERONET 
    print '**** Get %s from AERONET %02i-%02i-%4i' % (para, day, month, year)
    
    aerdir = '/nobackup/users/sunj/AERONET/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340
    
    if para == 'aod':
        aerfile = glob.glob( aerdir + '*.lev15')[0]
        aerdata = pd.read_csv(aerfile, sep=",", header = 4)
        
        ts = aerdata['Julian_Day']
        aod340 = aerdata['AOT_340']
        aod440 = aerdata['AOT_440']
        aod500 = aerdata['AOT_500']
        alpha1 = aerdata['500-870Angstrom']
        
    
        print '**** Applying Angstrom exponent to calculate AOD @ 550nm'
        paradata = aod500 * (550. / 500.) ** (-alpha1)
    
#        print '**** Plotting AERONET'
#        plt.figure()
#        plt.plot(ts, paradata, 'k', linewidth = 1)
#        plt.title('AOD')
#        plt.xlabel('Day of year')
        
        
    if para == 'AE':
        aerfile = glob.glob( aerdir + '*.lev15')[0]
        aerdata = pd.read_csv(aerfile, sep=",", header = 4)
        
        ts = aerdata['Julian_Day']
        alpha1 = aerdata['500-870Angstrom']
        
        paradata = alpha1
 
    if para == 'refidx':
        aerfile = glob.glob( aerdir + '*.rin')[0]
        aerdata = pd.read_csv(aerfile, sep=",", header = 3)    
        
        ts = aerdata['Julian_Day']
        refr440 = aerdata['REFR(440)']
        refi440 = aerdata['REFI(440)']
        refr675 = aerdata['REFR(675)']
        refi675 = aerdata['REFI(675)']
        
        
        ar =  (refr675 - refr440) / (675 - 440.)
        ai =  (refi675 - refi440) / (675 - 440.)
        
#        print ar, ai
        refr354 = refr440 - ar * (440 - 354.)      
        refi354 = refi440 - ai * (440 - 354.)   
        refr388 = refr440 - ar * (440 - 388.)      
        refi388 = refi440 - ai * (440 - 388.)   
        refr550 = refr440 + ar * (550 - 440.)
        refi550 = refi440 + ai * (550 - 440.)
    
        paradata = np.column_stack([refr354, refi354, refr388, refi388, refr550, refi550])
        
#        print '**** Plotting AERONET'
#        plt.figure(figsize = (12,4))
#        plt.subplot(1,2,1)
#        plt.plot(ts,refr354, 'b', linewidth = 1, label = '@ 354 nm')
#        plt.plot(ts,refr388, 'g', linewidth = 1, label = '@ 388 nm')
#        plt.plot(ts,refr550, 'r', linewidth = 1, label = '@ 550 nm')
#        plt.title('Real')
#        plt.xlabel('Day of year')
#
#        plt.subplot(1,2,2)
#        plt.plot(ts,refi354, 'b', linewidth = 1, label = '@ 354 nm')
#        plt.plot(ts,refi388, 'g', linewidth = 1, label = '@ 388 nm')
#        plt.plot(ts,refi550, 'r', linewidth = 1, label = '@ 550 nm')
#        plt.title('Imaginary')
#        plt.xlabel('Day of year')
#        plt.legend()
        
    if para == 'asy': 
        aerfile = glob.glob( aerdir + '*.asy')[0]
        aerdata = pd.read_csv(aerfile, sep=",", header = 3)    
        
        ts = aerdata['Julian_Day']     
        asy440 = aerdata['ASYM440-T']
        asy675 = aerdata['ASYM675-T']
        
        asy550 = (asy675 - asy440) / (675 - 440.) * (550 - 440.) + asy440        
        
        paradata = asy550 
        
#        print '**** Plotting AERONET'
#        plt.figure()
#        plt.plot(ts, asy440, 'k', linewidth = 1)
#        plt.plot(ts, asy675, 'b', linewidth = 1)
#        plt.plot(ts, asy550, 'r', linewidth = 1)
#        plt.title('Asmmetry factor')
#        plt.xlabel('Day of year')


    if para == 'ssa': 
        aerfile = glob.glob( aerdir + '*.ssa')[0]
        aerdata = pd.read_csv(aerfile, sep=",", header = 3)    
        
        ts = aerdata['Julian_Day']     
        ssa440 = aerdata['SSA440-T']
        ssa675 = aerdata['SSA675-T']
        
        ssa550 = (ssa675 - ssa440) / (675 - 440.) * (550 - 440.) + ssa440         
        
        paradata = ssa550
#        paradata = ssa440
        
#        print '**** Plotting AERONET'
#        plt.figure()
#        plt.plot(ts, ssa440, 'k', linewidth = 1)
#        plt.plot(ts, ssa675, 'b', linewidth = 1)
#        plt.plot(ts, ssa550, 'r', linewidth = 1)
#        plt.title('SSA')
#        plt.xlabel('Day of year')



    if para == 'sizfunc': 
        aerfile = glob.glob( aerdir + '*.siz')[0]
        aerdata = pd.read_csv(aerfile, sep=",", header = 3)    
        
        ts = aerdata['Julian_Day']  
        prob = np.array(aerdata)[:, 3:25]
        f = open(aerfile)
        lines = f.readlines()
        size = lines[3][43:242].split(',')
 
        
        for i,s in enumerate(size):
            size[i] = float(s)
        
#        print type(size)
        size = np.array(size)
        
        sizeint = np.arange(np.min(size),np.max(size),0.01) 
        
        fine = np.zeros(len(ts))
        coarse = np.zeros(len(ts))
        Cv_f = np.zeros(len(ts))
        Cv_c = np.zeros(len(ts))
        
        probint = np.zeros([len(prob),len(sizeint)])
        
        for i in range(len(prob)):
#            print i 
#        for i in range(337,338):
            f = interpolate.interp1d( size, np.array(list(prob[i])), kind ='linear')
            probint[i,:] = f (sizeint)
            peaks = argrelextrema(probint[i,:],np.greater) 
            
            
            fine[i] = sizeint[peaks[0][0]]
            coarse[i] = sizeint[peaks[0][-1]]
            
            fcidx = argrelextrema(probint[i,:],np.less)[0]
            
            if fcidx.size == 0:
                f = interpolate.interp1d( size, np.array(list(prob[i])), kind ='cubic')
                probint[i,:] = f (sizeint)
                
            fcidx = argrelextrema(probint[i,:],np.less)[0][0]
#            print fcidx
                
            Cv_f[i] = trapz(probint[i,0:fcidx+1],dx = 0.01)
            Cv_c[i] = trapz(probint[i,fcidx:-1], dx = 0.01) 
                        
            
            
        # covert volume distribution to number distribution
        # r_n = r_v / exp(3*sigma_v**2)
        # sigma_v = sigma_n = ln(sigma) 
        # sigma is the geometric standard deviation (also input in MMP and DISAMAR), typical range for aerosols is 1.5 to 2. 
        sigma_f = 1.5
        sigma_c = 2.0 
        
        r_f = fine / np.exp(3*(np.log(sigma_f)**2))
        r_c = coarse / np.exp(3*(np.log(sigma_c)**2))
        rfc_Cn = (Cv_f / Cv_c) * (r_c / r_f)**3 * np.exp(-4.5 * (sigma_f**2-sigma_c**2))  
        
        rfc_Cv = Cv_f / Cv_c
        
        wfc_n = rfc_Cn / (rfc_Cn + 1.) 
        wfc_v = rfc_Cv / (rfc_Cv + 1.)
         
        paradata = np.column_stack([r_f, r_c, wfc_n, wfc_v]) 
#            paradata.append(size[np.argmax(prob[i,0:12])])  
#        print paradata
#        print np.mean(paradata)
        
#        print '**** Plotting AERONET'
#        fig1 = plt.figure()
#        ax1 = fig1.add_axes([0.1,0.1,0.8,0.8])
##        for i in range(len(prob)): 
#        for i in range(332,333): 
#            plt.plot(size, prob[i],'x', linewidth = 1, label = 'AERONET')
#            plt.plot(sizeint, probint[i], linewidth = 1, label = 'interpolated')
#            plt.title('sizfunc')
#            plt.xlabel('r [mum]')
#            ax1.set_xscale('log')
#        plt.legend(frameon = False)
#    
    if month > 1:
        doy = dayspermon[0:month-1].sum() + day 
    else:
        doy = day 
        
    
    if para == 'aod' or para == 'AE':
        idx = np.where(ts == doy)[0][0]
    else:
#        print abs(ts-doy-13./24)
        idx = np.argmin(abs(ts-doy-13./24))
    
#    print idx
    
    
#    print np.argmax(paradata[idx])
#    print size[np.argmax(paradata[idx])]
        
    return paradata[idx]

#iday = 26
#r_f,r_c,w_n,w_v = AERONET_dailymean(2017,1,iday, 'sizfunc')
#print w_n,w_v

#aod550 = AERONET_dailymean(2017,1,iday, 'refidx')
#print aod550