import numpy as np
import matplotlib.pyplot as plt
 
import sys

sys.path.insert(1,'../sampling')
import hModelCircular


"""
compute h for the posterior samples
"""


def calculateH(posteriorSampes):

  MlowIntLimit  = 10.0**6.0	
  MhighIntLimit = 10.0**11.0
  zlowIntLimit  = 0.0 
  zhighIntLimit = 5.0

  hs = np.zeros(posteriorSamples)

  for i in range(len(hs)):

    theta = posteriorSamples[i][0],\
            posteriorSamples[i][1],\
            posteriorSamples[i][2],\
            posteriorSamples[i][3],\
            posteriorSamples[i][4]

    hs[i] = hModelCircular.hmodel(theta,\
                                  MlowIntLimit,MhighIntLimit,\
                                  zlowIntLimit,zhighIntLimit)


    return hs




# read in posterior samples

runLoc = '../../../../runDirs/run5/cpnest/'#'../../runDir/nlive10000/'
data = np.genfromtxt('{0}posterior.dat'.format(runLoc),names=True)
hout = open('{0}/hposterior.dat'.format(runLoc),'w')


hs = np.zeros(len(data['beta']))

for i in range(len(data['beta'])):

    theta = data['logn'][i],\
            data['beta'][i],\
            data['gamma'][i],\
            data['alpha'][i],\
            data['delta'][i]

    MlowIntLimit  = 10.0**6.0	
    MhighIntLimit = 10.0**11.0
    zlowIntLimit  = 0.0 
    zhighIntLimit = 5.0

    hs[i] = hModelCircular.hmodel(theta,\
                                  MlowIntLimit,MhighIntLimit,\
                                  zlowIntLimit,zhighIntLimit)

    hout.write('{}\n'.format(hs[i]))
                  
hout.close()
 



