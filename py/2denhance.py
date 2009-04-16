#Make density plot of Sommerfeld enhancement as a
#function of m_V/m and alpha, for fixed beta
#argument1: log_10 beta
#argument2: grid (default: 64)


from math import *
from numpy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
from enhance import enhance
import os, sys
import pickle
import matplotlib.cm as cm
from matplotlib import rc


#construct filenames
if len(sys.argv) > 1:
    logb= int(sys.argv[1])
else:
    print "No beta given, returning..."
    sys.exit(-1)

savefilename= '2d'+str(logb)+'.sav'
plotfilename= '2d'+str(logb)+'.eps'

if os.path.exists(savefilename):
    savefile=open(savefilename,'r+')
    data=pickle.load(savefile)
else:
    b=10**logb
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
    
    S= zeros((nms,nas))
    
    for ii in range(nms):
        for jj in range(nas):
            S[ii,jj]= log10(enhance(mv[ii],m,a[jj],b))
            sys.stdout.write('\r'+str(ii*nas+jj+1)+'/'+str(nms*nas))
            sys.stdout.flush()
    sys.stdout.write('\n')

    savefile=open(savefilename,'w')
    data={'mv':mv,'a':a,'S':S}
    pickle.dump(data,savefile)




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
im= imshow(data['S'],origin='lower',cmap=cm.gist_yarg,
       extent=[-3,-1,-5,-2],aspect=2./3.)
axis([-3,-1,-5,-2])
ylabel(r'$\log_{10}(m_{\phi}/m_\chi)$')#,fontsize=16)
xlabel(r'$\log_{10}\alpha$')#,fontsize=16)
text(-1.75,-4.75,r'$\beta = 10^{'+sys.argv[1]+'}$',color='w')#,fontsize=16)
CB1= colorbar(im,shrink=0.87)
levels=linspace(0,5,11)
#cont= contour(data['S'],levels,origin='lower',linewidths=1,colors='k',
#               extent=[-3,-1,-5,-2],aspect=2./3.)
#clabel(cont, fontsize=9, inline=1)
#CB2= colorbar(cont,shrink=0.8)
#Move the 2 colorbars on top of each other
#ll,bb,ww,hh = CB2.ax.get_position()
#CB1.ax.set_position([ll, bb, ww, hh])
CB1.set_label(r'$\log_{10}S$')#,fontsize=16)

savefig(plotfilename,format='eps')
