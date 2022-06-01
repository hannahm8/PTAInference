import numpy as np
import cpnest.model
from scipy import stats

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


def getKDE(sampleFile):

    # first part just counts the number of lines (samples) in the file
    count = 0
    with open(sampleFile) as fp:
        for line in fp:
            if line.strip():
                count += 1
       
    # then open the file again to read in the data
    A = np.zeros(count)
    burnIn = int(0.25*count)    
    
    # read in the lines and get the posterior samples for A
    for i in range(lenData):
        line  = f.readline().split()
        A[i] = line[215]
     
    # discard the early samples for burnin    
    ABurntIn = A[burnIn:]
     
    # and get the kde        
    kdeIPTA = stats.gaussian_kde(ABurntIn)
    return kdeIPTA



class IPTAResult(cpnest.model.Model):

    # parameter names and get flat prior bounds
    names=['logn','beta','gamma','alpha','delta']
    bounds=priorBounds()


    def log_likelihood(self, param):


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
      

        # use the kde likelihood here
        log10HModel = np.log10(float(hModel))
        log_like = kdeForIPTA(log10HModel)

        return log_like



posteriorFile = '/home/hannahm/Documents/PTA/IPTADR2Chains/iptadr2-fix_gamma_HD/chain_fg_HD.txt'
kdeForIPTA = getKDE(posteriorFile)


mymodel = IPTAResult()
nest = cpnest.CPNest(mymodel,maxmcmc=1000,nlive=10000,verbose=3)
nest.run()
cpnest.CPNest.get_posterior_samples(nest)

