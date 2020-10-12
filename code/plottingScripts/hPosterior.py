import numpy as np
import matplotlib.pyplot as plt
 
import sys

sys.path.insert(1,'../sampling')
import hModelCircular
import readPosterior


"""
compute h for the posterior samples
saves the h posterior in the run directory 
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






def main():

    # read in posterior samples
    pathToRuns = "../../runs/simpleModelPosteriors/"
    data = readPosterior.readPosterior(pathToRuns, nRuns=5)
    #runLoc = '../../../../runDirs/run5/cpnest/'#'../../runDir/nlive10000/'
    #data = np.genfromtxt('{0}posterior.dat'.format(runLoc),names=True)

    # where to save the hposterior? 
    hout = open('{0}/combined/hposterior.dat'.format(pathToRuns),'w')


    # to get the right length of the samples
    hs = np.zeros(len(data['beta']))
    print(len(hs))

    # compute h(f=1yr) for each posterior combo
    for i in range(len(hs)):

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
    print(hs)
    print(len(hs))





if __name__ == "__main__":
    main()
