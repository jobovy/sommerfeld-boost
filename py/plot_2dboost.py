#Makes density plot of total boost for a halo of mass Mh
#function of m_V/m and alpha
#argument1: log_10 Mh
#argument2: grid (default: 64)
#First calculate this density matrix using calc_2dboost.py!!!

from math import *
from numpy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
from enhance import avg_enhance
from boost import boost
import os, sys
import pickle
import matplotlib.cm as cm
from matplotlib import rc


#parameters
if len(sys.argv) > 1:
    logMh= int(sys.argv[1])
else:
    print "No Mh given, returning..."
    sys.exit(-1)
Mh=10**logMh
#vdisp= 2.5*10**-2*pow(Mh,1./3.)/(3*10**5)
#m=1.0
if len(sys.argv) > 2:
    nms= int(sys.argv[2])
    nas= int(sys.argv[2])
else:
    nms= 64
    nas= 64
nbins= 16#Hard-coded!!
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a

boostS= zeros((nas,nms))

#Loop over bins, restore and put in the right matrix
for bb in range(16):
    #construct filenames
    bin= bb+1
    restorefilename= '2dboost'+str(logMh)+'_'+str(nas)+'_'+str(bin)+'.sav'
    restorefile=open(restorefilename,'r')
    data=pickle.load(restorefile)
    S= data['S']
    #Load the boosts from S
    ii= 0
    jj= (bin-1)*nas/nbins
    while jj < bin*nas/nbins:
        while ii < nms:
            boostS[ii,jj]= S[ii,jj]
            ii= ii+1
        jj= jj+1
        ii= 0
restorefile.close()


plotfilename= '2dboost'+str(logMh)+'_'+str(nas)+'.eps'


#Plotting parameters
fig_width = 3.25  # width in inches
fig_height = 2.9      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize':10,
          'ytick.labelsize':10,
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)
#rc('text',usetex=True)
fig= figure()
im= imshow(boostS,origin='lower',#cmap=cm.gray,
       extent=[-3,-1,-5,-2],aspect=2./3.)
axis([-3,-1,-5,-2])
ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
xlabel(r'$\log_{10}\alpha$')#,fontsize=16)
text(-2,-4.75,r'$M = 10^{'+str(logMh)+r'} M_{\odot}$')#,fontsize=16)
CB1= colorbar(im,shrink=0.87)
levels=linspace(0,5,11)
#cont= contour(data['S'],levels,origin='lower',linewidths=1,colors='k',
#               extent=[-3,-1,-5,-2],aspect=2./3.)
#clabel(cont, fontsize=9, inline=1)
#CB2= colorbar(cont,shrink=0.8)
#Move the 2 colorbars on top of each other
#ll,bb,ww,hh = CB2.ax.get_position()
#CB1.ax.set_position([ll, bb, ww, hh])
CB1.set_label(r'$\log_{10}\mathcal{B}$')#,fontsize=16)

savefig(plotfilename,format='eps')



