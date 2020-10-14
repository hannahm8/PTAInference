import numpy as np
import mpmath
from scipy.integrate import quad


def constants():
    Gu     = 6.6738e-11
    cu     = 2.9979e8
    Mpcm   = 3.0857e22
    Gyrs   = 60.0*60.0*24.0*365.25*1.0e9
    Msolkg = 1.9891e30
	
    c = cu * Gyrs / Mpcm
    G = Gu * (Msolkg*(Gyrs**2.0)) / (Mpcm**3.0)
    #c,G=1.,1.
    fgw = 1.0e9  # 
    Ho = 1.0/13.8
    Mscale = 1.0e7

    constants = (4.0*(G**(5.0/3.0))*(Mscale**(5.0/3.0))) / (3.0*((np.pi)**(1.0/3.0))*(c**2.0)*(fgw**(4.0/3.0))*Ho)   
	
    return (constants**(1.0/2.0))






def hmodel(theta,Mlow_choice,Mhigh_choice,zlow_choice,zhigh_choice):

    logno, beta, gamma, alpha, delta = theta

    no = 10.0**(logno)
	
    # Msol cancels from t = M/(10^7Msol)
    Mlow=float(Mlow_choice)  # 10.0**6.0 #Msol
    Mhigh=float(Mhigh_choice) #10.0**11.0 #Msol
    tlow=Mlow/(10.0**(7.0+delta))
    thigh=Mhigh/(10.0**(7.0+delta))

    Mincomplete_gamma_fns = ( mpmath.gammainc( ((5.0/3.0)-alpha), tlow )  - mpmath.gammainc( ((5.0/3.0)-alpha), thigh )  )
	

    # zintegral from 0 to 5
    OmegaM = 0.3
    OmegaL = 0.7
    zlow=float(zlow_choice) #0.0 
    zhigh=float(zhigh_choice) #5.0


    zfn = lambda zi: (  ((1.0+zi)**(beta-(4.0/3.0)) * np.exp(-zi/gamma)) / ( (OmegaM*(1.0+zi)**3.0 + OmegaL)**(1.0/2.0) )  )
    Zint = quad(zfn,zlow,zhigh, epsabs=1.49e-12, epsrel=1.49e-12)

    h2 = no * ( 10.0**(delta*((5.0/3.0)-alpha)) ) * Mincomplete_gamma_fns * Zint[0]  * (1.0/np.log(10))


    h_no_constants = h2**(1.0/2.0)
    h = h_no_constants*constants()

    return h


