
"""
Download MODIS by wget 

@author: sunj
"""

import sys, os
import shutil
import glob
import time
import numpy as np

ordernr = "501175154"
platform = 'AQUA'
dataset = 'MYD04'
dayspermon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])


##url = "ftp://ladsweb.nascom.nasa.gov/orders/"+ordernr+"/"  
#url = "ftp://ladsweb.modaps.eosdis.nasa.gov/orders/"+ordernr+"/"   
#os.system('wget -r --retr-symlinks %s' % url)

#filelist = glob.glob("/nobackup/users/sunj/ladsweb.nascom.nasa.gov/orders/"+ordernr+"/*.hdf")   
filelist = glob.glob("/nobackup/users/sunj/ladsweb.modaps.eosdis.nasa.gov/orders/"+ordernr+"/*.hdf")
print filelist

for f in filelist:  
#    iyear = float(f[72:76])
#    doy = float(f[76:79])

    iyear = float(f[79:83])
    doy = float(f[83:86])
    
    if iyear%4 == 0:
        dayspermon[1] = 29 
    else:
        dayspermon[1] = 28
    
    for imon in range(0,12):
        if (doy > dayspermon[0:imon].sum()) & (doy <= dayspermon[0:imon+1].sum()):
            imonth = imon+1
            iday = doy - dayspermon[0:imon].sum()
    
    directory = '/nobackup/users/sunj/MODIS/'+platform+'/'+dataset+'/%4i/%02i/%02i/' % (iyear, imonth, iday)
        
    if not os.path.exists(directory):
        os.makedirs(directory)   
        
    if os.path.isfile(directory+f):
        os.remove(directory+f)
    
    shutil.move(f, directory)
