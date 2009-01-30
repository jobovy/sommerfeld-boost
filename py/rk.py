#Runge-Kutta routines

def rk4_step(y,t,step,deriv):
    """rk4_step: one Runge_Kutta step

    Input:
       y   - initial value
       t   - dependent variable value
       step- stepsize
       deriv(y,t)- callable function which evaluates to the derivatives of y at t

    Output:
       Updated value y
    """
    k1= step*deriv(y,t)
    k2= step*deriv(y+0.5*k1,t+0.5*step)
    k3= step*deriv(y+0.5*k2,t+0.5*step)
    k4= step*deriv(y+k3,t+step)

    #print k1, k2, k3, k4

    return y+1.0/6.0*(k1+2*k2+2*k3+k4)

def rk4(y,start,end,nsteps,deriv):
    """rk4: Runge-Kutta ODE solver

    Input:
       y     - Initial value
       start - dependent variable starting point
       end   - endpoint
       nsteps- number of steps to take
       deriv(y,t) - callable function which evaluates the derivative of y at t
       
    Output:
       Updated value of y at end, i.e., y(end)
    """
    #Determine the stepsize
    step= (end-start)/nsteps
    
    #Run rk4_step
    t=start
    import time

    i=0
    while i <nsteps:
        #print y
        y= rk4_step(y,t,step,deriv)
        t-= step
        i+=1
        #time.sleep(1)

    return y
