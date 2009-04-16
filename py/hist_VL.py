#Make relevant histograms of VL luminosities etc.
#argument1: log_10 mv
#argument2: log_10 alpha
#argument3: number of 8.5 kpc observers

from read_VL import read_VL
import numpy
from numpy import *
from scipy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
from matplotlib.axes import *
from matplotlib import rc
from math import log, cos, sin, sqrt, log10
from enhance import avg_enhance
from boost import boost
import pickle
import os, sys

#Read arguments
if len(sys.argv) > 1:
    mv= 10**double(sys.argv[1])
else:
    mv= 10**-2
if len(sys.argv) > 2:
    a= 10**double(sys.argv[2])
else:
    a= 10**-5
if len(sys.argv) > 3:
    nobs= int(sys.argv[3])
else:
    nobs= 3
obsR= 8.5

#Read VL data
VLdata= read_VL()

#main host halo
main= array([245452, 0., 207.218, 201.033, 56.7728, 1.93596e+12, 462.274, 0., 0., 0., 0., 0., 0., 1.32317e+08, 5.55052e+08])


#construct filenames
savefilename= 'hist_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'.sav'

if os.path.exists(savefilename):
    savefile=open(savefilename,'r+')
    data=pickle.load(savefile)
    sublum=data['sublum']
    sublumobs=data['sublumobs']
    totallum=data['totallum']
else:
    #first calculate the masses of the subhalos and sort by ascending halo mass
    masses= zeros(VLdata.shape[0])
    massfactor= log(2.)-1./2.
    for ii in range(VLdata.shape[0]):
        rs= VLdata[ii,4]/2.16258
        rhos= 8.56*10**4/pow(rs,2.)*pow(VLdata[ii,3],2.)
        masses[ii]= 4.*pi*rhos*pow(rs,3.)*massfactor
    #sort
    sort_indx= masses.argsort()
    for ii in range(VLdata.shape[1]):
        VLdata[:,ii]= VLdata[sort_indx,ii]
    #Now calculate D^2*luminosities for bound subhalos (that is, those within the virial radius
    lumdata=[]#for each subhalo holds: position(3) subhalo_luminosity total_boost somm_boost sub_boost rs rhos mass
    natt= 10
    nlums= 0
    massfactor= log(2.)-1./2.
    mo= 10**-6#starting value for boost mass
    bo= 0.
    for ii in range(VLdata.shape[0]):
        sys.stdout.write('\r'+str(ii+1)+'/'+str(VLdata.shape[0]))
        sys.stdout.flush()
        if VLdata[ii,1] > main[6] or VLdata[ii,1] == 0.0:#Use 800 kpc?
            continue
        nlums+= 1
        lumdata.append(VLdata[ii,7])
        lumdata.append(VLdata[ii,8])
        lumdata.append(VLdata[ii,9])
        rs= VLdata[ii,4]/2.16258
        if rs > VLdata[ii,6]:
            print "Warning: rs outside of tidal radius"
        rhos= 8.56*10**4/pow(rs,2.)*pow(VLdata[ii,3],2.)
        lumdata.append(4.63*10**-24*pow(rhos,2.)*pow(rs,3.)*7.*pi/6.)
        #calculate total mass within rs
        mass= 4.*pi*rhos*pow(rs,3.)*massfactor
        boo= boost(mass,mv,a,0,1,0,mo,bo)
        mo= mass
        bo= boo
        lumdata.append(boo)
        lumdata.append(boost(mass,mv,a,1))
        vdisp= 2.5*10**-2*pow(mass,1./3.)/(3*10**5)
        lumdata.append(avg_enhance(mv,1.,a,vdisp))
        lumdata.append(rs)
        lumdata.append(rhos)
        lumdata.append(mass)
    sys.stdout.write('\n')
    #Now put them in a nice matrix form
    print "Number of bound subhalos:"
    print nlums
    sublum= zeros((nlums,natt))
    for ii in range(nlums):
        for jj in range(natt):
            sublum[ii,jj]= lumdata[ii*natt+jj]
    #Calculate the actual luminosities and distances etc for different fiducial observers
    #sublumobs holds: distance_from_earth, luminosity, total_boost, reduced_boost, rs, rhos, mass
    natt2= 7
    sublumobs= zeros((nobs,nlums,natt2))
    for ii in range(nobs):
        #Generate random observer
        phi= numpy.random.uniform(0.,2.*pi)
        theta= numpy.random.uniform(0.,pi)
        pos= array([obsR*sin(theta)*sin(phi), obsR*sin(theta)*cos(phi), obsR*cos(theta)])
        for jj in range(nlums):
            sublumobs[ii,jj,0]= sqrt( pow((sublum[jj,0]-pos[0]),2.)+pow((sublum[jj,1]-pos[1]),2.)+pow((sublum[jj,2]-pos[2]),2.))
            if sublumobs[ii,jj,0] < 10.*sublum[jj,7]:
                print "Warning: subhalo is very close"
            sublumobs[ii,jj,1]= log10(sublum[jj,3]*sublum[jj,4]/pow(sublumobs[ii,jj,0],2))
            sublumobs[ii,jj,2]= log10(sublum[jj,4])
            sublumobs[ii,jj,3]= log10(sublum[jj,4]/(sublum[jj,6]*(1.+sublum[jj,5])))
            sublumobs[ii,jj,0]= log10(sublumobs[ii,jj,0])
            sublumobs[ii,jj,4]= sublum[jj,7]
            sublumobs[ii,jj,5]= sublum[jj,8]
            sublumobs[ii,jj,6]= sublum[jj,9]
    #Also generate 1,000 random observers and histogram the total luminosity
    nobs2= 10000
    totallum= zeros(nobs2)
    for ii in range(nobs2):
        #Generate random observer
        phi= numpy.random.uniform(0.,2.*pi)
        theta= numpy.random.uniform(0.,pi)
        pos= array([obsR*sin(theta)*sin(phi), obsR*sin(theta)*cos(phi), obsR*cos(theta)])
        for jj in range(nlums):
            dist2= pow((sublum[jj,0]-pos[0]),2.)+pow((sublum[jj,1]-pos[1]),2.)+pow((sublum[jj,2]-pos[2]),2.)
            totallum[ii]+= sublum[jj,3]*sublum[jj,4]/dist2
        totallum[ii]= log10(totallum[ii])
    #Save
    savefile=open(savefilename,'w')
    data={'sublum':sublum,'sublumobs':sublumobs,'totallum':totallum}
    pickle.dump(data,savefile)

#Calculate log masses
for ii in range(sublum.shape[0]):
    sublum[ii,9]= log10(sublum[ii,9])






savefilename2= 'hist_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_2.sav'

if os.path.exists(savefilename2):
    savefile2=open(savefilename2,'r+')
    data=pickle.load(savefile2)
    totallumnosomm=data['totallumnosomm']
else:
    #calculate the signal for the 1,000 observers without Sommerfeld enhancement
    nobs2= 1001
    totallumnosomm= zeros(nobs2)
    for ii in range(nobs2):
        #Generate random observer
        phi= numpy.random.uniform(0.,2.*pi)
        theta= numpy.random.uniform(0.,pi)
        pos= array([obsR*sin(theta)*sin(phi), obsR*sin(theta)*cos(phi), obsR*cos(theta)])
        for jj in range(sublum.shape[0]):
            dist2= pow((sublum[jj,0]-pos[0]),2.)+pow((sublum[jj,1]-pos[1]),2.)+pow((sublum[jj,2]-pos[2]),2.)
            totallumnosomm[ii]+= sublum[jj,3]*(1.+sublum[jj,6])/dist2
        totallumnosomm[ii]= log10(totallumnosomm[ii]) +log10(0.5/17.)
    #Save
    savefile2=open(savefilename2,'w')
    data={'totallumnosomm':totallumnosomm}
    pickle.dump(data,savefile2)

print numpy.sum(totallumnosomm)/numpy.sum(totallum)


#Post-process the fluxes to incorporate the new spectrum

#First do totallum
for ii in range(totallum.shape[0]):
    totallum[ii]= totallum[ii]+log10(0.5/17.)
#Then do the sublumobs
for ii in range(sublumobs.shape[0]):
    for jj in range(sublumobs.shape[1]):
        sublumobs[ii,jj,1]= sublumobs[ii,jj,1]+log10(0.5/30)

    
#Histogram all of this
plotfilename1= 'hist_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_1.eps'
plotfilename2= 'hist_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_2.eps'
plotfilename3= 'hist_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'_3.eps'

#Plotting parameters
fig_width = 3#3.25  # width in inches
fig_height =4.5      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'legend.fontsize': 8,
          'xtick.labelsize':10,
          'ytick.labelsize':10,
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)
#rc('text',usetex=True)
fig= figure(1)



nbins= 20

subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=None, hspace=0.4)
subplot(2,1,1)
hist(sublumobs[0,:,0],bins=nbins,histtype='step',ec='black')
hist(sublumobs[1,:,0],bins=nbins,histtype='step',ec='black',linestyle='dashed')
hist(sublumobs[2,:,0],bins=nbins,histtype='step',ec='black',linestyle='dashdot')
hist(sublumobs[3,:,0],bins=nbins,histtype='step',ec='black',linestyle='dotted')
axis([0,3,0,2000])
ylabel(r'$N_{\mbox{{\footnotesize sub}}}$')#,fontsize=16)
xlabel(r'$\log_{10}D$ $[\mathrm{kpc}]$')#,fontsize=16)

ax=subplot(2,1,2)
fig.subplots_adjust(hspace=.3)
hist(sublum[:,9],bins=nbins,histtype='step',ec='black')
axis([3,8,0,2000])
ylabel(r'$N_{\mbox{{\footnotesize sub}}}$')#,fontsize=16)
xlabel(r'$\log_{10}M(<r_{-2}) [M_\odot]$')#,fontsize=16)
#ax.set_yticklabels([])
savefig(plotfilename1,format='eps')

#Plotting parameters
fig_width = 3#3.25  # width in inches
fig_height = 4.5      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'legend.fontsize': 8,
          'xtick.labelsize':10,
          'ytick.labelsize':10,
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)
#rc('text',usetex=True)
fig= figure(2)
subplot(2,1,1)
hist(sublumobs[0,:,1],bins=nbins,histtype='step',ec='black')
hist(sublumobs[1,:,1],bins=nbins,histtype='step',ec='black',linestyle='dashed')
hist(sublumobs[2,:,1],bins=nbins,histtype='step',ec='black',linestyle='dashdot')
hist(sublumobs[3,:,1],bins=nbins,histtype='step',ec='black',linestyle='dotted')
axis([-16,-5,0,2000])
#ylabel(r'$N_{\mbox{{\footnotesize sub}}}$')#,fontsize=16)
ylabel(r'$N_{\mbox{{\footnotesize sub}}}$')#,fontsize=16)
xlabel(r'$\log_{10}\frac{\mathrm{d} N_{\gamma}}{\mathrm{d} \mathrm{A}\, \mathrm{d}t} [\mathrm{cm}^{-2}\mathrm{ s}^{-1}]$')#,fontsize=16)


ax=subplot(2,1,2)
fig.subplots_adjust(hspace=.4)
hist(sublumobs[0,:,2],bins=nbins,histtype='step',ec='black')
#axis([.5,2.5,0,2000])
#axis([6.25,7.25,0,2000])
axis([3.25,4.25,0,2000])
ylabel(r'$N_{\mbox{{\footnotesize sub}}}$')#,fontsize=16)
xlabel(r'$\log_{10}\mathcal{B}$')#,fontsize=16)
#ax.set_yticklabels([])


#ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
#xlabel(r'$\log_{10}\alpha$')#,fontsize=16)


savefig(plotfilename2,format='eps')

#Plotting parameters
fig_width = 3.25  # width in inches
fig_height = 3.4      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'legend.fontsize': 8,
          'xtick.labelsize':8,
          'ytick.labelsize':8,
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)

fig= figure(3)

#Set up different axes
ax1 = fig.add_subplot(111)
ax1.hist(totallum[0:1000],bins=nbins,histtype='step',ec='black')
#ax1.set_xlim(-9,-8)
#ax1.set_xlim(-3.5,-2.5)
ax1.set_xlim(-6.5,-4.5)
ax1.set_ylim(0,250)
ax1.set_ylabel(r'$N_{\mbox{{\footnotesize obs}}}$')#,fontsize=16)
ax1.set_xlabel(r'$\log_{10}\mathrm{Total }\frac{\mathrm{d} N_{\gamma}}{\mathrm{d} \mathrm{A}\, \mathrm{d}t} [\mathrm{cm}^{-2}\mathrm{ s}^{-1}]$')#,fontsize=16)

ax2=ax1.twiny()

ax2.hist(totallumnosomm[0:1000],bins=nbins,histtype='step',ec='black',linestyle='dashed')
ax2.set_xlabel(r'$\log_{10}\mbox{Total }\frac{\mathrm{d} N_{\gamma}}{\mathrm{d} \mathrm{A}\, \mathrm{d}t} [\mathrm{cm}^{-2}\mathrm{ s}^{-1}]$'+'(No Sommerfeld enhancement)')#,fontsize=16)
ax2.set_xlim(-10,-9)


savefig(plotfilename3,format='eps')
