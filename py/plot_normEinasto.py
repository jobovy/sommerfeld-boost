#Plot the functional dependence of the lum normalization, Einasto profile
from math import *
from numpy import *
import matplotlib
matplotlib.use('PS')
from pylab import *
from matplotlib.pyplot import *
from matplotlib import rc
from scipy.special import *

nsamples= 100

gammas= linspace(0.1,0.2,nsamples)
Ns= zeros(nsamples)
Fs= zeros(nsamples)

for ii in range(nsamples):
    x= gammas[ii]
    Ns[ii]= 4*pi*pow(64,-1./x)*pow(x,-1.+3./x)*exp(4./x)*gammainc(3./x,4./x)*gamma(3./x)
    Fs[ii]=4*pi*pow(64,-1./x)*pow(x,-1.+3./x)*exp(4./x)*gamma(3./x)
    

#gammas= array( [0.1, 0.101, 0.102, 0.103, 0.104, 0.105, 0.106, 0.107,
#0.108, 0.109,  0.11, 0.111, 0.112, 0.113, 0.114, 0.115, 0.116, 0.117,
#0.118, 0.119,  0.12, 0.121, 0.122, 0.123, 0.124, 0.125, 0.126, 0.127,
#0.128, 0.129,  0.13, 0.131, 0.132, 0.133, 0.134, 0.135, 0.136, 0.137,
#0.138, 0.139,  0.14, 0.141, 0.142, 0.143, 0.144, 0.145, 0.146, 0.147,
#0.148, 0.149,  0.15, 0.151, 0.152, 0.153, 0.154, 0.155, 0.156, 0.157,
#0.158, 0.159,  0.16, 0.161, 0.162, 0.163, 0.164, 0.165, 0.166, 0.167,
#0.168, 0.169,  0.17, 0.171, 0.172, 0.173, 0.174, 0.175, 0.176, 0.177,
#0.178, 0.179,  0.18, 0.181, 0.182, 0.183, 0.184, 0.185, 0.186, 0.187,
#0.188, 0.189,  0.19, 0.191, 0.192, 0.193, 0.194, 0.195, 0.196, 0.197,
#0.198, 0.199,  0.2] )

#Ns= array( [ 217.038, 212.896, 208.898, 205.038, 201.31, 197.707,
#194.224,  190.854, 187.593, 184.437, 181.379, 178.417, 175.546,
#172.761,  170.061, 167.44, 164.895, 162.424, 160.024, 157.691,
#155.424,  153.218, 151.073, 148.986, 146.955, 144.977, 143.05,
#141.173,  139.344, 137.562, 135.823, 134.128, 132.474, 130.861,
#129.285,  127.748, 126.246, 124.779, 123.346, 121.946, 120.577,
#119.239,  117.93, 116.651, 115.398, 114.173, 112.974, 111.801,
#110.652,  109.526, 108.424, 107.344, 106.286, 105.25, 104.234,
#103.237,  102.261, 101.303, 100.363, 99.4414, 98.5371, 97.6497,
#96.7787,  95.9238, 95.0844, 94.2602, 93.4508, 92.6558, 91.8748,
#91.1075,  90.3535, 89.6126, 88.8843, 88.1684, 87.4645, 86.7724,
#86.0919,  85.4225, 84.764, 84.1163, 83.4789, 82.8518, 82.2346,
#81.6272,  81.0292, 80.4406, 79.861, 79.2903, 78.7283, 78.1748,
#77.6296,  77.0925, 76.5633, 76.0419, 75.5282, 75.0219, 74.5229,
#74.031, 73.5461, 73.068, 72.5967] )

#Fs= array( [226.844, 222.685, 218.67, 214.794, 211.049, 207.43,
#203.93, 200.544,  197.267, 194.094, 191.021, 188.043, 185.156,
#182.356, 179.639,  177.003, 174.443, 171.957, 169.541, 167.193,
#164.911, 162.69,  160.531, 158.429, 156.382, 154.39, 152.449,
#150.557, 148.714,  146.917, 145.165, 143.455, 141.787, 140.16,
#138.571, 137.019,  135.503, 134.023, 132.576, 131.163, 129.78,
#128.429, 127.107,  125.814, 124.549, 123.31, 122.098, 120.912,
#119.75, 118.611, 117.496,  116.404, 115.333, 114.284, 113.255,
#112.247, 111.258, 110.287,  109.336, 108.402, 107.485, 106.586,
#105.703, 104.836, 103.984,  103.148, 102.327, 101.52, 100.728,
#99.9486, 99.1831, 98.4306,  97.6908, 96.9634, 96.2482, 95.5448,
#94.8529, 94.1723, 93.5026, 92.8437, 92.1953, 91.5572, 90.929,
#90.3106, 89.7018, 89.1023, 88.512,  87.9305, 87.3578, 86.7936,
#86.2378, 85.6902, 85.1505, 84.6186,  84.0945, 83.5778, 83.0684,
#82.5662, 82.0711, 81.5828, 81.1013] )

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
axUpper.set_ylabel(r'$\mathcal{N}_2(\alpha)$')#,fontsize=16)
axUpper.set_xticks( (0.1,.12,.14,.16,.18,.2) )
axUpper.set_xticklabels( () )

axLower.plot(gammas,Ns/Fs,color='black')
axLower.set_xlabel(r'$\alpha$')#,fontsize=16)
axLower.set_ylabel(r'$f(<r_{-2})$')#,fontsize=16)
axLower.set_xticks( axUpper.get_xticks())
axLower.set_xticklabels( (0.1,.12,.14,.16,.18,.2) )
axLower.set_yticks( (.85,.9,.95,1) )
#axLower.set_xticklabels( (0,0.5,1,1.5) )

#Make sure they have the same x-range
axUpper.axis([0.1,.2,65,220])
axLower.axis([0.1,.2,.85,1])


savefig("normgammaEinasto.eps",format='eps')
