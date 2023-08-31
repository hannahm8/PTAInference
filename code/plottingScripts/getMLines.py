import numpy as np
from scipy.integrate import quad

import readPosterior

def zInt(beta,gamma,zlow,zhigh):

    OmegaM = 0.3
    OmegaL = 0.7
    Ho = 1.0/13.8

    zfn = lambda zi: ( (1.0 + zi)**(beta-1.0) * np.exp(-zi/gamma) ) / \
                     ( Ho * ( OmegaM*(1.0+zi)**3. + OmegaL )**(1./2.) )
    zInt = quad(zfn,zlow,zhigh, epsabs=1.49e-12, epsrel=1.49e-12)
    
    return zInt[0]





def dndlog10M(theta,log10Ms):

    logno, beta, gamma, alpha, delta = theta

    no = 10.0**(logno)

    zIntegral = zInt(beta,gamma,0.0,5.0)

    Ms = [ 10.0**l10Ms for l10Ms in log10Ms ] 

    dndlogMs = [ no * zIntegral \
                 * (Mi / (10.**7.))**(-alpha) \
                 * np.exp(-Mi / 10.**(delta+7.0)) \
                 for Mi in Ms ]

    return dndlogMs
    






def main():

    pathToRuns = '../../runs/agnosticModel/logNormLikeMstar6to10/'
    nRuns=5
    simpleModelData = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)

    # set up M range and write to file
    log10MRange = np.linspace(6.0, 11.0, 40)
    writeMs = open('{}combined/ms_fix.dat'.format(pathToRuns),'w')
    for logM in log10MRange:
        writeMs.write('{}\n'.format(10.0**logM))
    writeMs.close()

    # place to write results 
    writeLines = open('{}combined/mLines_fix.dat'.format(pathToRuns),'w')    

    # compute lines
    for i in range(len(simpleModelData)):
        thetas = simpleModelData['logn'][i],  \
                 simpleModelData['beta'][i],  \
                 simpleModelData['gamma'][i], \
                 simpleModelData['alpha'][i], \
                 simpleModelData['delta'][i]
        ys = dndlog10M(thetas,log10MRange)   
        
        for yi in ys:        
            writeLines.write('{}\t'.format(yi))
        writeLines.write('\n')
    
    writeLines.close()
    

 
if __name__ == "__main__":
    main()

