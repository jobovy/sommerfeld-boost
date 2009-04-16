#Makes density plot of total boost for a halo of mass Mh
#function of m_V/m and alpha
#argument1: log_10 Mh
#argument2: grid (default: 64)
#argument3: DM density profile inner slope (GNFW or Einasto based on value)
#argument4: substructure abundance power-law exponent (default -0.9)
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
    nms= int(sys.argv[2])/2*3+1
    nas= int(sys.argv[2])+1
else:
    nms= 97
    nas= 65
if len(sys.argv) > 3:
    gamma= double(sys.argv[3])
    if gamma > 0.75:
        einasto= 0
    else:
        einasto= 1
else:
    gamma= 1.
    einasto= 0
if len(sys.argv) > 4:
    n= double(sys.argv[4])
else:
    n= -0.9
nbins= 16#Hard-coded!!
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a

boostS= zeros((nms,nas))

#Loop over bins, restore and put in the right matrix
for bb in range(nbins):
#    if (bb > 0 and bb < 8):
#        continue
    #construct filenames
    bin= bb+1
    restorefilename= '2dboost'+str(logMh)+'_'+str(nas)
    if len(sys.argv) > 3:
        restorefilename+= '_'+str(gamma)
    if len(sys.argv) > 4:
        restorefilename+= '_'+str(n)
    restorefilename+='_'+str(bin)+'.sav'
    restorefile=open(restorefilename,'r')
    data=pickle.load(restorefile)
    S= data['S']
    #Load the boosts from S
    ii= 0
    jj= (bin-1)*(nas-1)/nbins
    while jj < bin*(nas-1)/nbins:
        while ii < nms:
            boostS[ii,jj]= S[ii,jj]
            ii= ii+1
        jj= jj+1
        ii= 0
    restorefile.close()
#also take the last bin
restorefilename= '2dboost'+str(logMh)+'_'+str(nas)
if len(sys.argv) > 3:
    restorefilename+= '_'+str(gamma)
if len(sys.argv) > 4:
    restorefilename+= '_'+str(n)
restorefilename+='_'+str(nbins+1)+'.sav'
restorefile=open(restorefilename,'r')
data=pickle.load(restorefile)
S= data['S']
#Load the boosts from S
ii= 0
jj= nas-1
while ii < nms:
    boostS[ii,jj]= S[ii,jj]
    ii= ii+1
restorefile.close()


#Then plot
plotfilename= '2dboost'+str(logMh)+'_'+str(nas)
if len(sys.argv) > 4:
    plotfilename+= '_'+str(n)
plotfilename+='.eps'

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
im= imshow(boostS,origin='lower',cmap=cm.gist_yarg,
       extent=[-3,-1,-5,-2],aspect=2./3.)
axis([-3,-1,-5,-2])
ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
xlabel(r'$\log_{10}\alpha$')#,fontsize=16)
text(-1.75,-4.75,r'$n$ = '+str(n),color='w')#,fontsize=16)
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



