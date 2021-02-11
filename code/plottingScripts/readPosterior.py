"""
for reading in posterior samples from the agnostic model
"""

import numpy as np


def readPosterior(runPath, nRuns=5):

    # takes however many runs were done and reads in all the posteriors together

    for i in range(5):

       runLoc = '{}run{}/'.format(runPath,int(i+1))
       
       if i==0: 
         data  = np.atleast_1d(np.genfromtxt('{0}posterior.dat'.format(runLoc),names=True))
         length = np.shape(data)[0]
       else: 
         dataNew  = np.atleast_1d(np.genfromtxt('{0}posterior.dat'.format(runLoc),names=True))
         data  = np.concatenate((data,dataNew))

         length+=np.shape(dataNew)[0]


    return data







