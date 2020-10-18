import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def findSigma(low,high,sigmas,mu):

    diffLow,diffHigh=10.,10.
    ratioLow,ratioHigh=100,100
    
    frac=0.05
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
                #print(sig,ratioLow, ratioHigh)
                sigSave = sig
    return sigSave





def Gaussian(x,mu,sigma):

    y = np.exp(-(x-mu)**2. / (2.*sigma*sigma))

    return y





def getGaussianLikeParams():

    #######################
    # Gaussian likelihood #       
    #######################

    # 5-95% from NANOGrav result
    low  = 1.37
    high = 2.67
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





def getLogNormLikeParams():


    low  = np.log10(1.37E-15)
    high = np.log10(2.67E-15)
    mu = (low+high)/2.

    print(low,high,mu)

    sigmas = np.arange(0.01,0.5,0.0001)
   
    sigmaResult = findSigma(low,high,sigmas,mu)

    print('sigma is ', sigmaResult )

    xs=np.arange(low-0.1,high+0.1,0.0001)
    print(xs)
    y=Gaussian(xs,mu,sigmaResult)
    print (y)

    plt.clf()
    plt.axvline(low)
    plt.axvline((high))
    plt.plot(xs,y)
    plt.show()



    return None





getLogNormLikeParams()



    
