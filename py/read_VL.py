#Read the Via Lactea subhalo data

import re
from numpy import *

def read_VL(path='/home/users/jb2777/sommboost/VL/vltwosubs.txt'):
    """Read the Via Lactea II subhalo data

    Input:
       path    - path to the data

    Output:
       array with all of the data in it
    """
    VLfile= open(path,'r')

    expr= re.compile(r"-?[0-9]+(\.[0-9]*)?(e\+?-?[0-9]+)?")

    rawdata= []
    nline= 0
    nvalues= 15#number of columns in the VL data
    for line in VLfile:
        if line[0] == '#':
            continue
        nline+= 1
        values= expr.finditer(line)
        nvalue= 0
        for i in values:
            rawdata.append(float(i.group()))
            nvalue+= 1
        if nvalue != nvalues:
            print "Warning, nvalues not equal to 15"

    VLdata= zeros((nline,nvalues))
    for ii in range(nline):
        for jj in range(15):
            VLdata[ii,jj]= rawdata[ii*nvalues+jj]

    return VLdata

                           
                                 
        
