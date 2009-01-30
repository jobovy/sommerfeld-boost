#! /usr/bin/env python
from enhance import enhance
from math import exp, cos, sin, sqrt
from scipy import integrate
from numpy import *
from matplotlib.pyplot import *
import os
import pickle


nb=5
nsamples=20

savefilename='1d.sav'
if os.path.exists(savefilename):
    savefile=open(savefilename,'r+')
    data=pickle.load(savefile)
else:
    mv=.09
    a=1/30.0
    
    m= linspace(0,2,nsamples)
    m= 10**m
    b= array([.1,.01,.001,.0001,.00001])
    S= zeros((nsamples,nb))

    for ii in range(nsamples):
        for jj in range(nb):
            S[ii,jj]= enhance(mv,m[ii],a,b[jj])
            print str(ii*nb+jj)+'/'+str(nb*nsamples)

    savefile=open(savefilename,'w')
    data={'m':m,'S':S}
    pickle.dump(data,savefile)


            
for jj in range(nb):
    loglog(data['m'],data['S'][:,jj])

axis([1,100,1,10**6])
xlabel('m (TeV)')
ylabel('S')
savefig("enhancement.eps",format='eps')
