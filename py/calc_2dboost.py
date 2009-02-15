#Calculate the boost for a halo, by manually parallalizing the calculation
#Makes matrix of total boost for a halo of mass Mh
#function of m_V/m and alpha
#argument1: log_10 Mh
#argument2: grid (default: 64), should be divisible by 16!
#argument3: bin number (bins go vertical)
#argument4: Einasto alpha, or gamma NFW (based on the value)
#argument5: slope of the subhalo abundance power-law (default -0.9)

from math import *
from numpy import *
from enhance import avg_enhance
from boost import boost
import os, sys
import pickle

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
    nms= int(sys.argv[2])/2*3+1
    nas= int(sys.argv[2])+1
else:
    nms= 97
    nas= 65
nbins= 16#Hard-coded!!, actually +1
if len(sys.argv) > 3:
    bin= int(sys.argv[3])
else:
    bin= 1
if len(sys.argv) > 4:
    gamma= double(sys.argv[4])
    if gamma > 0.75:
        einasto= 0
    else:
        einasto= 1
else:
    gamma= 1.
    einasto= 0
if len(sys.argv) > 5:
    n= double(sys.argv[5])
    #Also recalculate the normalization of the subhalo abundance function
    q= 0.1
    g= 10**-5
    if n != -1:
        A= 0.1*(n+1.)/(pow(q,n+1.)-pow(g,n+1.))
    else:
        A= 0.1/log(q/g)
else:
    n= -0.9
    q= 0.1
    g= 10**-5
    A= 0.1*(n+1.)/(pow(q,n+1.)-pow(g,n+1.))
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a
    
#construct filenames
savefilename= '2dboost'+str(logMh)+'_'+str(nas)
if len(sys.argv) > 4:
    savefilename+= '_'+str(gamma)
if len(sys.argv) > 5:
    savefilename+= '_'+str(n)
savefilename+='_'+str(bin)+'.sav'

if os.path.exists(savefilename):
    savefile=open(savefilename,'r+')
    data=pickle.load(savefile)
    ii= data['ii']
    jj= data['jj']
    if data['nms'] != nms:
        print "Warning: suplied nms does not equal savefile nms"
        sys.exit(-1)
    if data['nas'] != nas:
        print "Warning: suplied nas does not equal savefile nas"
        sys.exit(-1)
    mv= data['mv']
    a= data['a']
    S= data['S']
    print 'restarting at ii= '+str(ii)+', jj= '+str(jj)+', at iteration '+str(ii+(jj-(bin-1)*(nas-1)/nbins)*nms+1)+'/'+str(nms*nas/nbins)
else:
    S= zeros((nms,nas))
    ii= 0
    if bin != 17:
        jj= (bin-1)*(nas-1)/nbins
    else:
        jj= nas-1
    savefile=open(savefilename,'w')

if bin != 17:
    while jj < bin*(nas-1)/nbins:
        while ii < nms:
            boo= boost(Mh,mv[ii],a[jj],0,gamma,einasto,10**-6,0.,n,0.1,A)
            so= avg_enhance(mv[ii],m,a[jj],vdisp)
            S[ii,jj]= log10(so+boo)
            sys.stdout.write('\r'+str(ii+(jj-(bin-1)*(nas-1)/nbins)*nms+1)+'/'+str(nms*(nas-1)/nbins))
            sys.stdout.flush()
            ii= ii+1
            data={'mv':mv,'a':a,'S':S,'nms':nms,'nas':nas,'ii':ii,'jj':jj}
            savefile.seek(0,0)
            pickle.dump(data,savefile)
            savefile.flush()
        jj= jj+1
        ii= 0
    sys.stdout.write('\n')
else:
    while ii < nms:
        boo= boost(Mh,mv[ii],a[jj],0,gamma,einasto,10**-6,0.,n)
        so= avg_enhance(mv[ii],m,a[jj],vdisp)
        S[ii,jj]= log10(so+boo)
        sys.stdout.write('\r'+str(ii+(jj-(bin-1)*(nas-1)/nbins)*nms+1)+'/'+str(nms*(nas-1)/nbins))
        sys.stdout.flush()
        ii= ii+1
        data={'mv':mv,'a':a,'S':S,'nms':nms,'nas':nas,'ii':ii,'jj':jj}
        savefile.seek(0,0)
        pickle.dump(data,savefile)
        savefile.flush()
    sys.stdout.write('\n')
    
savefile.close()


