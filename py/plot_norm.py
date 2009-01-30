#Plot the functional dependence of the lum normalization
from math import *
from numpy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
from matplotlib import rc

nsamples= 100

gammas= linspace(0,1.5,nsamples,False)
Ns= zeros(nsamples)
for ii in range(nsamples):
    x= gammas[ii]
    Ns[ii]= 4*pi*pow(2,2*x-5)*(-16+11*x-2*pow(x,2))/(-30+47*x-24*pow(x,2)+4*pow(x,3))

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
plot(gammas,Ns,color='black')
axis([0,1.5,0,25])
xlabel(r'$\gamma$')#,fontsize=16)
ylabel(r'$\mathcal{N}(\gamma)$')#,fontsize=16)

savefig("normgamma.eps",format='eps')
