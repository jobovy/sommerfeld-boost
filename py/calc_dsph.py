#Calculate the flux for a dSph, by manually parallalizing the calculation
#Makes matrix of total flux as a
#function of m_V/m and alpha
#argument1: dsph
#argument2: grid (default: 64), should be divisible by 16!
#argument3: bin number (bins go vertical)
#argument4: Einasto alpha, or gamma NFW (based on the value)

from math import *
from numpy import *
from enhance import avg_enhance
from boost import boost
import os, sys
import pickle

#parameters
if len(sys.argv) > 1:
    if sys.argv[1] == "fornax":
        rs= 1.85
        rhos= 10**7
        M= 1.54*10**8
        D= 138.
        vdisp= 10.5
    elif sys.argv[1] == "leoI":
        rs= 2.77
        rhos= 10**7
        M= 5.16*10**8
        D= 250.
        vdisp= 8.8
    elif sys.argv[1] == "sculptor":
        rs= 2.77
        rhos= 10**7
        M= 5.16*10**8
        D= 79.
        vdisp= 6.6
    elif sys.argv[1] == "leoII":
        rs= 1.85
        rhos= 10**7
        M= 1.54*10**8
        D= 205.
        vdisp= 6.7
    elif sys.argv[1] == "sextans":
        rs= 0.92
        rhos= 10**7
        M= 1.9*10**7
        D= 86.
        vdisp= 6.6
    elif sys.argv[1] == "carina":
        rs= 1.39
        rhos= 1.13*10**7
        M= 7.37*10**7
        D= 101.
        vdisp= 6.8
    elif sys.argv[1] == "ursaminor":
        rs= 2.77
        rhos= 10**7
        M= 5.15*10**8
        D= 66.
        vdisp= 9.3
    elif sys.argv[1] == "draco":
        rs= 3.7
        rhos= 10**7
        M= 1.23*10**8
        D= 82.
        vdisp= 9.5
    elif sys.argv[1] == "coma":
        rs= 0.3
        rhos= 3.*10**8
        M= 1.97*10**8
        D= 44.
        vdisp= 4.6
    elif sys.argv[1] == "umaii":
        rs= 0.3
        rhos= 3*10**8
        M= 1.97*10**8
        D= 32.
        vdisp= 6.7
else:
    print "No dsph given, returning..."
    sys.exit(-1)
vdisp= vdisp/(3*10**5)
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
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a
    
#construct filenames
savefilename= 'flux_'+sys.argv[1]+'_'+str(nas)
if len(sys.argv) > 4:
    savefilename+= '_'+str(gamma)
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
            boo= boost(M,mv[ii],a[jj],0,gamma,einasto)
            so= avg_enhance(mv[ii],m,a[jj],vdisp)
            S[ii,jj]= log10((so+boo)*7.*pi/6.*pow(rhos,2.)*pow(rs,3.)/pow(D,2)*4.63*10**-24)
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
        boo= boost(M,mv[ii],a[jj],0,gamma,einasto)
        so= avg_enhance(mv[ii],m,a[jj],vdisp)
        S[ii,jj]= log10((so+boo)*7.*pi/6.*pow(rhos,2.)*pow(rs,3.)/pow(D,2)*4.63*10**-24)
        sys.stdout.write('\r'+str(ii+(jj-(bin-1)*(nas-1)/nbins)*nms+1)+'/'+str(nms*(nas-1)/nbins))
        sys.stdout.flush()
        ii= ii+1
        data={'mv':mv,'a':a,'S':S,'nms':nms,'nas':nas,'ii':ii,'jj':jj}
        savefile.seek(0,0)
        pickle.dump(data,savefile)
        savefile.flush()
    sys.stdout.write('\n')
    
savefile.close()


