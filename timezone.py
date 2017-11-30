# -*- coding: utf-8 -*-
"""
Time zone calcualtion

@author: sunj
"""

def timezone(lon): 
    if abs(lon)%15 < 7.5:
        jetlag = abs(lon)//15
    else:
        jetlag = abs(lon)//15+1
        
    if lon<0:
        jetlag = - jetlag 
        
    print 'Time zone:', jetlag 
    
    return jetlag 

