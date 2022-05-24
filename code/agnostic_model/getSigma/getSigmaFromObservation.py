"""
Estimates a mu / sigma given a 5%-95% confidence region. 
The values are used to estimate a Gaussian distribution 
which is consistent with the NANOGrav 12.5 yr observation 
at 1/1yr. 
"""


import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def findSigma(low,high,sigmas,mu,credibleRange):

    diffLow,diffHigh=10.,10.
    ratioLow,ratioHigh=100,100
    
    frac=0.5*(100.-credibleRange)/100. #0.05
    print(frac)
    tol = 0.01

    lowerLimit = low-10.
    upperLimit = high+10

    # loop over sigma values to find a place where 0.05 of the distribution
    # is above and below the 5-95% region from NanoGrav
    for sig in sigmas:


        totalInt = quad(Gaussian,lowerLimit,upperLimit,args=(mu,sig))

        lowInt   = quad(Gaussian,lowerLimit,low,args=(mu,sig))

        highInt  = quad(Gaussian,high,upperLimit,args=(mu,sig)) 
    
        #print(totalInt,lowInt)
        try:
          ratioLow  = lowInt[0]/totalInt[0]
          ratioHigh = highInt[0]/totalInt[0]
        except: pass

        if ratioLow<(frac+tol)  and \
           ratioHigh<(frac+tol) and \
           ratioLow>(frac-tol)  and \
           ratioHigh>(frac-tol):
            if abs(ratioLow-0.05)<diffLow or abs(ratioHigh-0.05)<diffHigh:
                sigSave = sig
    return sigSave





def Gaussian(x,mu,sigma):

    y = np.exp(-(x-mu)**2. / (2.*sigma*sigma))

    return y



# unused - now using logNorm for better fit
"""
def getGaussianLikeParams(low,high):

    #######################
    # Gaussian likelihood #       
    #######################

    # 5-95% from NANOGrav result
    #low  = 1.37
    #high = 2.67
    mu = (low + high)/2.

    # sigma range
    sigmas = np.arange(0.01,.5,0.0001)


    sigmaResult = findSigma(low,high,sigmas,mu)

    print('sigma is ',sigmaResult)

    # plot results 
    xs=np.arange(low*0.1,high*2.,0.01)
    x = np.atleast_1d([np.log10(xi*1E-15) for xi in xs])
    y=Gaussian(xs,mu,sigmaResult)

    plt.clf()
    plt.axvline(np.log10(low*1E-15))
    plt.axvline(np.log10(high*1E-15))
    plt.plot(x,y)
    plt.show()


    return None
"""




def getLogNormLikeParams(lowIn,highIn,CR,label=''):


    low  = np.log10(lowIn*(10.**-15))
    high = np.log10(highIn*(10.**-15))
    mu = (low+high)/2.

    #print(low,high,mu)
    print("mu is ", mu )

    sigmas = np.arange(0.01,0.5,0.0001)
   
    sigmaResult = findSigma(low,high,sigmas,mu,CR)

    print('sigma is ', sigmaResult )

    xs=np.arange(low-0.1,high+0.1,0.0001)
    print(xs)
    y=Gaussian(xs,mu,sigmaResult)
    print (y)

    plt.axvline(low)
    plt.axvline((high))
    plt.plot(xs,y,label=label)
    



    return None




plt.clf()

lowNANOGrav = 1.37
highNANOGrav = 2.67

getLogNormLikeParams(lowNANOGrav,highNANOGrav,90,label='N')


lowPPTA = 2.2-0.3
highPPTA = 2.2+0.4

getLogNormLikeParams(lowPPTA,highPPTA,68,label='P')



plt.legend()
plt.show()


    
