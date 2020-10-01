import numpy as np
import cpnest.model

import hModelCircular


# ranges on the priors
def priorBounds():
	
    logno = [-20.0,3.0]
    beta  = [-2.0,7.0]
    gamma = [0.2,5.0]
    alpha = [-3.0,3.0]
    delta = [-1.0,2.0]

    priors = [logno,beta,gamma,alpha,delta]

    return priors




class NANOGravResult(cpnest.model.Model):

    # parameter names and get flat prior bounds
    names=['logn','beta','gamma','alpha','delta']
    bounds=priorBounds()


    def log_likelihood(self, param):

        # based on NANOGrav result
        hCentre = 2.02 #1.9E-15
        hSigma  = 0.42 #1.E-16

        # integration limits
        MlowIntLimit  = 10.0**6.0	
        MhighIntLimit = 10.0**11.0
        zlowIntLimit  = 0.0 
        zhighIntLimit = 5.0

        # calculate h
        theta = param['logn'],param['beta'],param['gamma'],param['alpha'],param['delta']
        hModel = hModelCircular.hmodel(theta,\
                                         MlowIntLimit,MhighIntLimit,\
                                         zlowIntLimit,zhighIntLimit)
        
        # scale h for likelihood
        hModelOver1EMinus15 = hModel/(1.E-15)

        # likelihood 
        log_like = -0.5 *( (hModelOver1EMinus15-hCentre)**2. ) / (hSigma*hSigma)

        return log_like



mymodel = NANOGravResult()
nest = cpnest.CPNest(mymodel,maxmcmc=1000,nlive=10000,verbose=3)
nest.run()
cpnest.CPNest.get_posterior_samples(nest)

