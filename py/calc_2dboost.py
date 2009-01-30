#Calculate the boost for a halo, by manually parallalizing the calculation
#Makes matrix of total boost for a halo of mass Mh
#function of m_V/m and alpha
#argument1: log_10 Mh
#argument2: grid (default: 64)
#argument3: bin number (bins go vertical)

from math import *
from numpy import *
#import matplotlib
#matplotlib.use('PS')
#from pylab import *
#from matplotlib.pyplot import *
from enhance import avg_enhance
from boost import boost
import os, sys
import pickle
#import matplotlib.cm as cm
#from matplotlib import rc


#parameters
if len(sys.argv) > 1:
    logMh= int(sys.argv[1])
else:
    print "No Mh given, returning..."
    sys.exit(-1)
Mh=10**logMh
vdisp= 2.5*10**-2*pow(Mh,1./3.)/(3*10**5)
m=1.0
if len(sys.argv) > 2:
    nms= int(sys.argv[2])
    nas= int(sys.argv[2])
else:
    nms= 64
    nas= 64
nbins= 16#Hard-coded!!
if len(sys.argv) > 3:
    bin= int(sys.argv[3])
else:
    bin= 1
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a
    
#construct filenames
savefilename= '2dboost'+str(logMh)+'_'+str(nas)+'_'+str(bin)+'.sav'



if os.path.exists(savefilename):
    savefile=open(savefilename,'r+')
    data=pickle.load(savefile)
    ii= data['ii']
    jj= data['jj']
    if data['nms'] != nms:
        print "Warning: suppplied nms does not equal savefile nms"
        sys.exit(-1)
    if data['nas'] != nas:
        print "Warning: suppplied nas does not equal savefile nas"
        sys.exit(-1)
    mv= data['mv']
    a= data['a']
    S= data['S']
    print 'restarting at ii= '+str(ii)+', jj= '+str(jj)+', at iteration '+str(ii*nms+jj+1)+'/'+str(nms*nas/nbins)
else:
    S= zeros((nas,nms))
    ii= 0
    jj= (bin-1)*nas/nbins
    savefile=open(savefilename,'w')

while jj < bin*nas/nbins:
    while ii < nms:
        boo= boost(Mh,mv[ii],a[jj])
        so= avg_enhance(mv[ii],m,a[jj],vdisp)
#        print mv[ii]
#        print a[jj]
#        print boo
#        print so
        S[ii,jj]= log10(so+boo)
#        print S[ii,jj]
        sys.stdout.write('\r'+str(ii+(jj-(bin-1)*nas/nbins)*nas/nbins+1)+'/'+str(nms*nas/nbins))
        sys.stdout.flush()
        ii= ii+1
        data={'mv':mv,'a':a,'S':S,'nms':nms,'nas':nas,'ii':ii,'jj':jj}
        savefile.seek(0,0)
        pickle.dump(data,savefile)
        savefile.flush()
    jj= jj+1
    ii= 0
sys.stdout.write('\n')

savefile.close()


