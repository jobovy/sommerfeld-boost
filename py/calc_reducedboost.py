#Calculate the reduced boost for a halo, i.e., totalboost/(somm(Mh)*substboost)
#Makes matrix of reduced boost for a halo of mass Mh, and plots it
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
    nms= int(sys.argv[2])/2*3+1
    nas= int(sys.argv[2])+1
else:
    nms= 97
    nas= 65
mv= linspace(-5,-2,nms)
mv= 10**mv
a= linspace(-3,-1,nas)
a= 10**a

boostS= zeros((nms,nas))

#########################################################
#
# First get the total boosts from the right files
#
#########################################################

#Loop over bins, restore and put in the right matrix
nbins= 16
for bb in range(nbins):
    #construct filenames
    bin= bb+1
    restorefilename= '2dboost'+str(logMh)+'_'+str(nas)+'_'+str(bin)+'.sav'
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
restorefilename= '2dboost'+str(logMh)+'_'+str(nas)+'_'+str(17)+'.sav'
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


#Then calculate the substructure boost for this halo mass, without the Sommerfeld enhancement
subboost= boost(Mh,1,1,1)

#Then calculate the reduced Sommerfeld enhancement for each grid point
#filename
savefilename= '2dredboost'+str(logMh)+'_'+str(nas)+'.sav'

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
    redboostS= data['redboostS']
    print 'restarting at ii= '+str(ii)+', jj= '+str(jj)+', at iteration '+str(ii*nas+jj+1)+'/'+str(nms*nas)
else:
    redboostS= zeros((nms,nas))
    ii= 0
    jj= 0
    savefile=open(savefilename,'w')

while ii < nms:
    while jj < nas:
        somm= avg_enhance(mv[ii],m,a[jj],vdisp)
        redboostS[ii,jj]= boostS[ii,jj]-log10(somm*(1.0+subboost))
        #print boostS[ii,jj], log10(somm*(1+subboost))
        sys.stdout.write('\r'+str(ii*nas+jj+1)+'/'+str(nms*nas))
        sys.stdout.flush()
        jj= jj+1
        data={'mv':mv,'a':a,'redboostS':redboostS,'nms':nms,'nas':nas,'ii':ii,'jj':jj}
        savefile.seek(0,0)
        pickle.dump(data,savefile)
        savefile.flush()
    ii= ii+1
    jj= 0
sys.stdout.write('\n')


#Now plot the result
plotfilename= '2dredboost'+str(logMh)+'_'+str(nas)+'.eps'

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
im= imshow(redboostS,origin='lower',cmap=cm.gist_yarg,
       extent=[-3,-1,-5,-2],aspect=2./3.)
axis([-3,-1,-5,-2])
ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
xlabel(r'$\log_{10}\alpha$')#,fontsize=16)
text(-2,-4.75,r'$M = 10^{'+str(logMh)+r'} M_{\odot}$')#,color='w')#,fontsize=16)
CB1= colorbar(im,shrink=0.87)
levels=linspace(0,5,11)
#cont= contour(data['S'],levels,origin='lower',linewidths=1,colors='k',
#               extent=[-3,-1,-5,-2],aspect=2./3.)
#clabel(cont, fontsize=9, inline=1)
#CB2= colorbar(cont,shrink=0.8)
#Move the 2 colorbars on top of each other
#ll,bb,ww,hh = CB2.ax.get_position()
#CB1.ax.set_position([ll, bb, ww, hh])
CB1.set_label(r'$\log_{10}\left[\frac{\mathcal{B}}{S(M)(1+B_0(M))}\right]$')#,fontsize=16)

savefig(plotfilename,format='eps')


