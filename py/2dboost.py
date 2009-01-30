#Make density plot of total boost for a halo of mass Mh
#function of m_V/m and alpha
#argument1: log_10 Mh
#argument2: grid (default: 64)

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
vdisp= 2.5*10**-2*pow(Mh,1./3.)/(3*10**5)
m=1.0
if len(sys.argv) > 2:
    nms= int(sys.argv[2])
    nas= int(sys.argv[2])
else:
    nms= 64
    nas= 64
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a
    
#construct filenames
savefilename= '2dboost'+str(logMh)+'_'+str(nas)+'.sav'
plotfilename= '2dboost'+str(logMh)+'_'+str(nas)+'.eps'


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
    print 'restarting at ii= '+str(ii)+', jj= '+str(jj)+', at iteration '+str(ii*nms+jj+1)+'/'+str(nms*nas)
else:
    S= zeros((nas,nms))
    ii= 0
    jj= 0
    savefile=open(savefilename,'w')

while ii < nms:
    while jj < nas:
        boo= boost(Mh,mv[ii],a[jj])
        so= avg_enhance(mv[ii],m,a[jj],vdisp)
#        print mv[ii]
#        print a[jj]
#        print boo
#        print so
        S[ii,jj]= log10(so+boo)
#        print S[ii,jj]
        sys.stdout.write('\r'+str(ii*nms+jj+1)+'/'+str(nms*nas))
        sys.stdout.flush()
        jj= jj+1
        data={'mv':mv,'a':a,'S':S,'nms':nms,'nas':nas,'ii':ii,'jj':jj}
        savefile.seek(0,0)
        pickle.dump(data,savefile)
        savefile.flush()
    ii= ii+1
    jj= 0
sys.stdout.write('\n')



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
im= imshow(S,origin='lower',#cmap=cm.gray,
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
