from math import exp, cos, sin, sqrt
from scipy import integrate
from numpy import *
from pylab import plot, show
from matplotlib.pyplot import *

def deriv(y,t,p1,p2):
    """Implements the right hand side of the Sommerfeld enhancement set of differential equations (dy/dt)

    Input:
       y    - current dependent variable  value
       t    - current independent variable
       p1   - first parameter = b**2/a**2
       p2   - second parameter = mv/(a*m)

    Output:
       derivative
    """
    return array([y[1],
                    (-p1-1.0/t*exp(-p2*t))*y[0],
                    y[3],
                    (-p1-1.0/t*exp(-p2*t))*y[2]])

def enhance(mv,m,a,b):
    """Calculates the Sommerfeld enhancement for a given set of parameters

    Input:
       mv    - mass of the force carrier
       m     - DM mass
       a     - new force fine structure constant
       b     - v/c

    Output:
       Sommerfeld enhancement
    """
    #Integration parameters
    par1= b**2/a**2
    par2= mv/(a*m)

    #Determine start
    if (10/par1 < 1/par2):#Coulomb like
        start= 10**2/par1
    else:
        start= 10**2/par2

    end=10**-4
    y0=array([cos(b/a*start),-b/a*sin(b/a*start),
              sin(b/a*start),b/a*cos(b/a*start)])
    time=[start,end]
    y=integrate.odeint(deriv,y0,time,args=(par1,par2),mxstep=10**8,mxhnil=0)

    return 1/(y[1,0]**2+y[1,2]**2)

def avg_enhance(mv,m,a,b,nsamples=100,vesc=0):
    """Calculates the Sommerfeld enhancement averaged over a truncated Maxwell-Boltzman distribution

    Input:
       mv       - mass of the force carrier
       m        - DM mass
       a        - new force fine structure constant
       b        - v/c
       nsamples - number of equidistant samples to use in the averaging; default 100
       vesc     - escape velocity (truncation velocity); default sqrt(2)*b

    Output:
       Averaged Sommerfeld enhancement
    """
    avg= 0.
    norm= 0.

    if vesc == 0.:
        vesc= 8.*b#2.*sqrt(2.)*b

    bs= linspace(vesc/(nsamples-1.),vesc,nsamples-1)

    for bb in bs:
        temp_norm= bb**2*exp(-bb**2/(4.*b**2))
        avg= avg+ temp_norm*enhance(mv,m,a,bb)
        norm= norm+temp_norm

    return avg/norm
        
