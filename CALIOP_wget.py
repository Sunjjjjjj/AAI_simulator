 
"""
Download CALIOP by wget 

@author: sunj
"""

import sys, os
import shutil
import glob
import time
import numpy as np

ordernr = "201700624323910139674"
#dataset = 'CAL_LID_L2_VFM-ValStage1-V3-40'
dataset = 'L2-05kmAPro_Prov-V3-40'
#dataset = 'CAL_LID_L1-Standard-V4-10'

dayspermon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

#url = 'ftp://xfr140.larc.nasa.gov/'+ordernr+'/*.hdf'
url = 'ftp://xfr139.larc.nasa.gov/'+ordernr+'/*.hdf'
os.system('wget -r --retr-symlinks %s' % url)

filelist = glob.glob("/nobackup/users/sunj/xfr139.larc.nasa.gov/"+ordernr+"/*.hdf")

for f in filelist:  
    idx = f.find('T')
    iyear = float(f[idx-10:idx-6])
    imonth = float(f[idx-5:idx-3])
    iday = float(f[idx-2:idx])
    
#    if iyear%4 == 0:
#        dayspermon[1] = 29 
#    else:
#        dayspermon[1] = 28
#    
#    for imon in range(0,12):
#        if (doy > dayspermon[0:imon].sum()) & (doy <= dayspermon[0:imon+1].sum()):
#            imonth = imon+1
#            iday = doy - dayspermon[0:imon].sum()
    
    directory = '/nobackup/users/sunj/CALIPSO/'+dataset+'/%4i/%02i/%02i/' % (iyear, imonth, iday)
        
    if not os.path.exists(directory):
        os.makedirs(directory)   
        
    if os.path.isfile(directory+f):
        os.remove(directory+f)
    
    shutil.move(f, directory)
