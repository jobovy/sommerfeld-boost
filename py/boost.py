from scipy import special
from numpy import *
from math import exp, log
from enhance import avg_enhance, enhance
from scipy import integrate

def deriv(y,t,mv,a,gamma,einasto,n,q,A,C,K1,K2,vdisp_sat,somm_sat):
    """Implements the right hand side of the boost differential equation

    Input:
       y    - current dependent variable  value
       t    - current independent variable
       mv      - force carrier mass (/mchi)
       a       - new force fine structure constant
       gamma   - GNFW gamma or Einasto gamma
       einasto - if 0 GNFW, 1 Einasto
       n       - mass function exponent
       q       - mass fraction until which to integrate
       A       - normalization of mass function
       C       - normalization of sigma(m) relation
       K1      - normalization of concentration-mass relation
       K2      - exponent of concentratio-mass relation
       vdisp_sat - velocity dispersion at which the Sommerfeld enhancement saturates
       somm_sat  - saturated value of the Sommerfeld enhancement

    Output:
       derivative
    """
    #Calculate each factor/term in turn
    #First calculate the concentration
    conc= K1*pow(t,K2)/pow(10.,K2*12)
    concq= K1*pow(q*t,K2)/pow(10,K2*12)
    
    #First up, L(qM)/L(M)
    LL= pow(q,3.*K2+1.)
    if einasto == 1:
        fc= special.gammainc(3./gamma,2./gamma*pow(conc,gamma))
        LL*= pow(fc,2)#we don't need to add gamma(3/gamma) since this would get divided out
        LL/= pow(special.gammainc(3./gamma,2./gamma*pow(concq,gamma)),2)#we don't need to add gamma(3/gamma) since this would get divided out
        fc*= special.gamma(3./gamma)
        #dlnLdlnM= 1.+3.*K2+2.*conc/fc*K2*(pow(2,3./gamma)*pow(gamma,1.-3./gamma)*pow(conc,2.)*exp(-2./gamma*pow(conc,gamma)))
        dlnLdlnM= 1.+3.*K2-2.*conc/fc*K2*(pow(conc,2.)*exp(2./gamma*(1.-pow(conc,gamma))))
        print -2.*conc/fc*K2*(pow(conc,2.)*exp(2./gamma*(1.-pow(conc,gamma))))
        print special.gammainc(3./gamma,2./gamma*pow(conc,gamma))
        print fc
        print (pow(conc,2.)*exp(2./gamma*(1.-pow(conc,gamma)))), (pow(conc,2.-gamma)*pow(1.+conc,gamma-3.))
        print dlnLdlnM
    else:
        if gamma == 1:
            fc= log(1.+conc)-conc/(1.+conc)
            LL*= pow(fc,2)#This is only valid for NFW
            LL/= pow(log(1+concq)-concq/(1.+concq),2)#This is only valid for NFW
        else:
            fc= special.hyp2f1(3.-gamma,3.-gamma,4.-gamma,-conc)*pow(conc,3.-gamma)/(3.-gamma)
            LL*= pow(fc,2)
            LL/= pow(special.hyp2f1(3.-gamma,3.-gamma,4.-gamma,-concq)*pow(concq,3.-gamma)/(3.-gamma),2)
        dlnLdlnM= 1.+3.*K2-2.*conc/fc*K2*(pow(conc,2.-gamma)*pow(1.+conc,gamma-3.))
#        print (pow(conc,2.)*exp(2./gamma*(1.-pow(conc,gamma)))), (pow(conc,2.-gamma)*pow(1.+conc,gamma-3.))
#        print -2.*conc/fc*K2*(pow(conc,2.-gamma)*pow(1.+conc,gamma-3.))
#        print dlnLdlnM

    vdisp= 2.5*10**-2*pow(t,1./3.)/(3*10**5)
    if vdisp <= vdisp_sat:
        somm= somm_sat
        #somm= 1
    else:
        somm= avg_enhance(mv,1.,a,vdisp,10)
        #somm= 1
        #somm= enhance(mv,1.,a,vdisp)

    #Now that we have everything, calculate dB/dM
    return 1./(1.+(1.-q)*LL) * 1./t * (A*pow(q,n)*(somm+y)*LL - n*y-y*dlnLdlnM)  


def boost(Mh,mv=10**-2,a=10**-1,nosomm=0,gamma=1,einasto=0,mo=10**-6,bo=0.,n=-0.9,q=0.1,A=2.3*10**-2,C=2.5*10**-2,K1=10.5,K2=-.11):
    """Calculates the boost factor for a given halo mass

    Input:
       Mh      - Halo mass
       mv      - force carrier mass (/mchi)
       a       - new force fine structure constant
       gamma   - GNFW gamma or Einasto gamma
       einasto - if 0 GNFW, 1 Einasto
       mo      - free-streaming limit
       bo      - boost at mo
       n       - mass function exponent
       q       - mass fraction until which to integrate
       A       - normalization of mass function
       C       - normalization of sigma(m) relation
       K1      - normalization of concentration-mass relation
       K2      - exponent of concentration-mass relation
       nosomm  - DON'T include Sommerfeld enhancement

    Output:
       boost
    """

    if nosomm == 0:
        #To speed up the calculation, first find where the Sommerfeld enhancement saturates, and pass that value to the differential equation to use for smaller velocities.
        vdisp_Mh= 2.5*10**-2*pow(Mh,1./3.)/(3*10**5)
        vdisp_mo= 2.5*10**-2*pow(mo,1./3.)/(3*10**5)
        vdisps=[vdisp_Mh,vdisp_Mh*10**-1,vdisp_Mh*10**-2,vdisp_Mh*10**-3,vdisp_Mh*10**-4,vdisp_Mh*10**-5,vdisp_Mh*10**-6,vdisp_Mh*10**-7,vdisp_Mh*10**-8,vdisp_Mh*10**-9,vdisp_Mh*10**-10]
        #vdisps=[10**-2,10**-3,10**-4,10**-5,10**-6,10**-7,10**-8,10**-9,10**-10,10**-11,10**-12]
        somms=zeros(11)
        indx=0
        for ii in range(11):
            if vdisps[ii] < vdisp_mo:
                break
            somms[ii]= avg_enhance(mv,1.,a,vdisps[ii])
            if ii > 0:
                if somms[ii] <= 1.1*somms[ii-1]:
                    break
                else:
                    indx+= 1
        vdisp_sat= vdisps[indx]
        somm_sat= somms[indx]
    else:
        vdisp_sat=10*Mh
        somm_sat= 1


    start= mo
    end= Mh
    y0= bo
    time=[start,end]
    y= integrate.odeint(deriv,y0,time,args=(mv,a,gamma,einasto,n,q,A,C,K1,K2,vdisp_sat,somm_sat),mxstep=10**8)

    return y[1]
