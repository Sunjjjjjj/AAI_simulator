
"""
Download OMI by wget 

@author: sunj
"""

import sys, os
import shutil
import glob
import time
import numpy as np

dataset = 'OMAERO.003'
dayspermon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

for iyear in range(2017,2018):
    
    if iyear%4 == 0:
        dayspermon[1] = 29 
    else:
        dayspermon[1] = 28
    
    for doy in range(230,234):
        
        print 'start downloading day-of-year: ', doy
        
        url = "https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/"+dataset+"/%4i/%03i/" % (iyear, doy)
        os.system('wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -r -c -nH -nd -np -A he5 %s' % url)
        
        for imon in range(0,12):
            if (doy > dayspermon[0:imon].sum()) & (doy <= dayspermon[0:imon+1].sum()):
                imonth = imon+1
                iday = doy - dayspermon[0:imon].sum()
                print imonth, iday
        
        directory = '/nobackup/users/sunj/OMI/OMAERO/%4i/%02i/%02i/' % (iyear, imonth, iday)
        
        if not os.path.exists(directory):
            os.makedirs(directory)   
            
        filelist = glob.glob('*.he5')
        
        start=time.time() 
        for f in filelist:  
            if os.path.isfile(directory+f):
                os.remove(directory+f)
            shutil.move(f, directory)
        end = time.time()    
        print 'Time period of range downloaded files:',end-start,'s'
