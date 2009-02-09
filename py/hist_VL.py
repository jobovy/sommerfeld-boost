#Make relevant histograms of VL luminosities etc.
#argument1: log_10 mv
#argument2: log_10 alpha
#argument3: number of 8.5 kpc observers

from read_VL import read_VL
from numpy import *
from scipy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
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
        lumdata.append(4.63*10**-9*pow(rhos,2.)*pow(rs,3.)*7.*pi/6.)
        #calculate total mass within rs
        mass= 4.*pi*rhos*pow(rs,3.)*massfactor
        boo= boost(mass,mv,a,0,1,0,mo,bo)
        mo= mass
        bo= boo
        lumdata.append(boo)
        lumdata.append(boost(mass,mv,a,1))
        vdisp= 2.5*10**-2*pow(mass,1./3.)/(3*10**5)
        lumdata.append(1.)#avg_enhance(mv,1.,a,vdisp))
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
        phi= random.uniform(0.,2.*pi)
        theta= random.uniform(0.,pi)
        pos= array([obsR*sin(theta)*sin(phi), obsR*sin(theta)*cos(phi), obsR*cos(theta)])
        for jj in range(nlums):
            sublumobs[ii,jj,0]= sqrt( pow((sublum[jj,0]-pos[0]),2.)+pow((sublum[jj,1]-pos[1]),2.)+pow((sublum[jj,2]-pos[2]),2.))
            if sublumobs[ii,jj,0] < 10.*sublum[jj,7]:
                print "Warning: subhalo is very close"
            sublumobs[ii,jj,1]= log10(sublum[jj,3]*sublum[jj,4]/pow(sublumobs[ii,jj,0],2))
            sublumobs[ii,jj,2]= log10(sublum[jj,4])
            sublumobs[ii,jj,3]= log10(sublum[jj,4]/(sublum[jj,5]*(1.+sublum[jj,6])))
            sublumobs[ii,jj,0]= log10(sublumobs[ii,jj,0])
            sublumobs[ii,jj,4]= sublum[jj,7]
            sublumobs[ii,jj,4]= sublum[jj,8]
            sublumobs[ii,jj,4]= sublum[jj,9]
    #Save
    savefile=open(savefilename,'w')
    data={'sublum':sublum,'sublumobs':sublumobs}
    pickle.dump(data,savefile)
  
#Histogram all of this
plotfilename= 'hist_'+sys.argv[1]+'_'+sys.argv[2]+'_'+sys.argv[3]+'dist.eps'

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
hist(sublumobs[0,:,0],histtype='step')
#axis([-3,-1,-5,-2])
#ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
#xlabel(r'$\log_{10}\alpha$')#,fontsize=16)


savefig(plotfilename,format='eps')

