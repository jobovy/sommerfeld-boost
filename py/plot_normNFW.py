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
Fs= zeros(nsamples)

for ii in range(nsamples):
    x= gammas[ii]
    Ns[ii]= 4*pi*pow(2,2*x-5)*(-16+11*x-2*pow(x,2))/(-30+47*x-24*pow(x,2)+4*pow(x,3))
    Fs[ii]=Ns[ii]/(4*pi)*(30-47*x+24*pow(x,2)-4*pow(x,3))
    
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


figBG   = 'w'        # the figure background color
axesBG  = 'w'  # the axies background color

def get_locator():
        """
            the axes cannot share the same locator, so this is a helper
                function to generate locators that have identical functionality
                    """

        return IndexLocator(10, 1)

left, width = 0.1, 0.8
rect1 = [left, 0.3, width, 0.6]
rect2 = [left, 0.1, width, 0.2]
axUpper      = axes(rect1, axisbg=axesBG)  #left, bottom, width, height
axLower     = axes(rect2, axisbg=axesBG)#, sharex=axUpper)

axUpper.xaxis.set_major_locator( get_locator() )

axUpper.plot(gammas,Ns,color='black')
axUpper.set_ylabel(r'$\mathcal{N}_1(\gamma)$')#,fontsize=16)
axUpper.set_xticks( (0,0.5,1,1.5) )
axUpper.set_xticklabels( () )

axLower.plot(gammas,Fs,color='black')
axLower.set_xlabel(r'$\gamma$')#,fontsize=16)
axLower.set_ylabel(r'$f(<r_s)$')#,fontsize=16)
axLower.set_xticks( axUpper.get_xticks() )
axLower.set_xticklabels( axUpper.get_xticks() )
axLower.set_yticks( (0.4,0.6,.8,1) )
#axLower.set_xticklabels( (0,0.5,1,1.5) )

#Make sure they have the same x-range
axUpper.axis([0,1.5,0,25])
axLower.axis([0,1.5,0.4,1])


savefig("normgammaNFW.eps",format='eps')
