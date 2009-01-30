from pylab import plot, show, semilogx, loglog
from enhance import deriv
from scipy import integrate
from numpy import array,linspace
from math import cos,sin, exp, log, log10, fabs
from matplotlib.pyplot import *

#parameters
a=1/30.
b=.0001
mv=.00001#9
m=68

par1= b**2/a**2
par2= mv/(a*m)


#Determine start
if (10/par1 < 1/par2):#Coulomb like
    start= 10**2/par1
    print "Coulomb-like"
else:
    start= 10**2/par2
    print "Finite range"

print start
end=10**-4
y0=array([cos(b/a*start),-b/a*sin(b/a*start),
          sin(b/a*start),b/a*cos(b/a*start)])
numsteps=10000
time=linspace(log10(start),log10(end),numsteps)
time=10**time

y=integrate.odeint(deriv,y0,time,args=(par1,par2),mxstep=10**8)

print log10(1/(y[numsteps-1,0]**2+y[numsteps-1,2]**2)), 1/(y[numsteps-1,0]**2+y[numsteps-1,2]**2)
print log10(1/(y[numsteps-1,0]**2)), 1/(y[numsteps-1,0]**2)
print log10(1/(y[numsteps-1,2]**2)), 1/(y[numsteps-1,2]**2)
print (y[numsteps-1,3]**2)*a**2/b**2
print (y[numsteps-1,1]**2)*a**2/b**2


#extrapolate?
print (2/(log(10.)*fabs(y[numsteps-1,0]))*y[numsteps-1,1]*end-2*log10(fabs(y[numsteps-1,0])))
print (2/(log(10.)*fabs(y[numsteps-1,2]))*y[numsteps-1,3]*end-2*log10(fabs(y[numsteps-1,2])))


loglog(time,(y[:,0])**2)
loglog(time,(y[:,2])**2)
loglog(time,(y[:,1])**2)
loglog(time,(y[:,3])**2)
axis([end,a**2/b**2,b**2,1.5])
show()
