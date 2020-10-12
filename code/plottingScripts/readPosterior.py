import numpy as np


def readPosterior(runPath, nRuns=5):

    # takes however many runs were done and reads in all the posteriors together

    for i in range(5):

       runLoc = '{}run{}/'.format(runPath,int(i+1))
       
       if i==0: 
         data  = np.atleast_1d(np.genfromtxt('{0}posterior.dat'.format(runLoc),names=True))
         #hData = np.atleast_1d(np.genfromtxt('{0}hposterior.dat'.format(runLoc)))
         #print(np.shape(data),np.shape(hData))
         length = np.shape(data)[0]
       else: 
         dataNew  = np.atleast_1d(np.genfromtxt('{0}posterior.dat'.format(runLoc),names=True))
         #hDataNew = np.atleast_1d(np.genfromtxt('{0}hposterior.dat'.format(runLoc)))
         #print(np.shape(dataNew),np.shape(hDataNew))

         data  = np.concatenate((data,dataNew))
         #hData = np.concatenate((hData,hDataNew))

         length+=np.shape(dataNew)[0]

         #print(np.shape(data),np.shape(hData))     

    #print(data['logn'])
    #print(hData)

    return data







