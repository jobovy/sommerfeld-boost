#Makes density plot of total flux from a dsph
#function of m_V/m and alpha
#argument1: dsph
#argument2: grid (default: 64)
#argument3: Einasto alpha, or gamma NFW (based on the value)
#First calculate this density matrix using calc_dsph.py!!!

from math import *
from numpy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
from matplotlib.pyplot import title as pytitle
from enhance import avg_enhance
from boost import boost
import os, sys
import pickle
import matplotlib.cm as cm
from matplotlib import rc


#parameters
if len(sys.argv) > 1:
    if sys.argv[1] == "fornax":
        rs= 1.85
        rhos= 10**7
        M= 1.54*10**8
        D= 138.
        vdisp= 10.5
        title="Fornax"
    elif sys.argv[1] == "leoI":
        rs= 2.77
        rhos= 10**7
        M= 5.16*10**8
        D= 250.
        vdisp= 8.8
        title="Leo I"
    elif sys.argv[1] == "sculptor":
        rs= 2.77
        rhos= 10**7
        M= 5.16*10**8
        D= 79.
        vdisp= 6.6
        title="Sculptor"
    elif sys.argv[1] == "leoII":
        rs= 1.85
        rhos= 10**7
        M= 1.54*10**8
        D= 205.
        vdisp= 6.7
        title="Leo II"
    elif sys.argv[1] == "sextans":
        rs= 0.92
        rhos= 10**7
        M= 1.9*10**7
        D= 86.
        vdisp= 6.6
        title="Sextans"
    elif sys.argv[1] == "carina":
        rs= 1.39
        rhos= 1.13*10**7
        M= 7.37*10**7
        D= 101.
        vdisp= 6.8
        title="Carina"
    elif sys.argv[1] == "ursaminor":
        rs= 2.77
        rhos= 10**7
        M= 5.15*10**8
        D= 66.
        vdisp= 9.3
        title="Ursa Minor"
    elif sys.argv[1] == "draco":
        rs= 3.7
        rhos= 10**7
        M= 1.23*10**8
        D= 82.
        vdisp= 9.5
        title="Draco"
    elif sys.argv[1] == "coma":
        rs= 0.3
        rhos= 3.*10**8
        M= 1.97*10**8
        D= 44.
        vdisp= 4.6
        title="Coma Berenices"
    elif sys.argv[1] == "umaii":
        rs= 0.3
        rhos= 3*10**8
        M= 1.97*10**8
        D= 32.
        vdisp= 6.7
        title="Ursa Major II"
else:
    print "No dsph given, returning..."
    sys.exit(-1)
if len(sys.argv) > 1:
    print "Plotting "+sys.argv[1]
else:
    print "No Mh given, returning..."
    sys.exit(-1)
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
nbins= 16#Hard-coded!!
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a

flux= zeros((nms,nas))

#Loop over bins, restore and put in the right matrix
for bb in range(nbins):
#    if (bb > 0 and bb < 8):
#        continue
    #construct filenames
    bin= bb+1
    restorefilename= 'flux_'+sys.argv[1]+'_'+str(nas)
    if len(sys.argv) > 3:
        restorefilename+= '_'+str(gamma)
    restorefilename+='_'+str(bin)+'.sav'
    restorefile=open(restorefilename,'r')
    data=pickle.load(restorefile)
    S= data['S']
    #Load the boosts from S
    ii= 0
    jj= (bin-1)*(nas-1)/nbins
    while jj < bin*(nas-1)/nbins:
        while ii < nms:
            flux[ii,jj]= log10(0.5/17.)+S[ii,jj]#This is some post-processing because I changed the DM annihilation spectrum
            ii= ii+1
        jj= jj+1
        ii= 0
    restorefile.close()
#also take the last bin
restorefilename= 'flux_'+sys.argv[1]+'_'+str(nas)
if len(sys.argv) > 3:
    restorefilename+= '_'+str(gamma)
restorefilename+='_'+str(nbins+1)+'.sav'
restorefile=open(restorefilename,'r')
data=pickle.load(restorefile)
S= data['S']
#Load the boosts from S
ii= 0
jj= nas-1
while ii < nms:
    flux[ii,jj]= log10(0.5/17.)+S[ii,jj]#This is some post-processing because I changed the DM annihilation spectrum
    ii= ii+1
restorefile.close()


#also calculate the background in this range
back= log10(2.69156*10**-7*2*pi*(1.-cos(rs/D)))

#Then plot
plotfilename= 'flux_'+sys.argv[1]+'_'+str(nas)
if len(sys.argv) > 3:
    plotfilename+= '_'+str(gamma)
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
im= imshow(flux,origin='lower',cmap=cm.gist_yarg,
       extent=[-3,-1,-5,-2],aspect=2./3.)
axis([-3,-1,-5,-2])
ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
xlabel(r'$\log_{10}\alpha$')#,fontsize=16)
CB1= colorbar(im,shrink=0.87)
levels=linspace(0,5,11)
#cont= contour(data['S'],levels,origin='lower',linewidths=1,colors='k',
#               extent=[-3,-1,-5,-2],aspect=2./3.)
#clabel(cont, fontsize=9, inline=1)
#CB2= colorbar(cont,shrink=0.8)
#Move the 2 colorbars on top of each other
#ll,bb,ww,hh = CB2.ax.get_position()
#CB1.ax.set_position([ll, bb, ww, hh])
CB1.set_label(r'$\log_{10}\frac{\mathrm{d} N_{\gamma}}{\mathrm{d} \mathrm{A}\, \mathrm{d}t} [\mathrm{cm}^{-2}\mathrm{ s}^{-1}]$')#,fontsize=16)
pytitle(title,x=0.5,y=1.1,fontsize=10)
suptitle(r'{\scriptsize $\log_{10}\frac{\mathrm{d} N_{\gamma,\mbox{{\tiny background}}}}{\mathrm{d} \mathrm{A}\, \mathrm{d}t} = '+'%4.1f' % back+r' \log_{10} [\mathrm{cm}^{-2}\mathrm{ s}^{-1}]$}',x=0.45,y=0.93)

savefig(plotfilename,format='eps')



