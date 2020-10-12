from matplotlib import pylab
import numpy as np
from math import ceil
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})


import readPosterior





def histOutline(dataIn, *args, **kwargs):
    """
    code copied from http://www.scipy.org/Cookbook/Matplotlib/UnfilledHistograms?action=AttachFile&do=get&target=histNofill.py

    Make a histogram that can be plotted with plot() so that
    the histogram just has the outline rather than bars as it
    usually does.

    Example Usage:
    binsIn = numpy.arange(0, 1, 0.1)
    angle = pylab.rand(50)

    (bins, data) = histOutline(binsIn, angle)
    plot(bins, data, 'k-', linewidth=2)

    """

    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)

    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)    
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0
    
    return (bins, data)







def plotOneDPosterior(pathToRuns,nRuns,paramName,paramLabel):

    plt.clf()


    #pathToRuns = "../../runs/simpleModelPosteriors/"
    data = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)

    setBins = np.arange(int(min(data[paramName])), ceil(max(data[paramName])),0.5)

    p05 = np.percentile(data[paramName],5)
    p50 = np.percentile(data[paramName],50)
    p95 = np.percentile(data[paramName],95)


    (bins,dat) = histOutline(data[paramName],setBins)
    pylab.plot(bins,dat,'k-',linewidth=2,alpha=0.7)
    plt.axvline(p05,ls=':',color='k')
    plt.axvline(p95,ls=':',color='k')
    plt.xlabel(paramLabel)
    plt.ylabel('Number of counts')
    plt.tight_layout()
    #plt.savefig('{}lognoHist.pdf'.format(runLoc))
    #plt.savefig('{}lognoHist.png'.format(runLoc))
    pylab.show()


    print("""{}
    05%: {}
    50%: {}
    95%: {}
    """.format(paramName,p05,p50,p95))

    return 0




def plotHDistribution(pathToH):

    plt.clf()

    hData = np.genfromtxt(pathToH)
    log10hData = np.log10(hData)

    ph05 = np.percentile(log10hData,5)
    ph50 = np.percentile(log10hData,50)
    ph95 = np.percentile(log10hData,95)
    print (10**ph05, 10**ph95)


    hSetBins = np.arange(int(min(log10hData)),ceil(max(log10hData)),0.05)
    #print(hSetBins)
    (hbins,hdat) = histOutline(log10hData,hSetBins)
    pylab.plot(hbins,hdat,'k-',linewidth=2,alpha=0.7)
    plt.xlim(-15.5,-14.2)
    plt.axvline(ph05,ls=':',color='k')
    plt.axvline(ph95,ls=':',color='k')
    plt.axvline(ph50,ls="--",color='k')
    plt.ylabel('Number of counts')
    plt.xlabel(r'$A_{\rm yr}$')
    plt.tight_layout()
    #plt.savefig('{}hposterior.png'.format(runLoc))
    #plt.savefig('{}hposterior.pdf'.format(runLoc))
    plt.show()
    

    return 0









def main():

    parameterNames = (['logn','beta','gamma','alpha','delta'])
    parameterLabels = ([r'$\log_{10}\frac{{\dot n}_0}{{\rm Mpc}^{-3}{\rm Gyr}}$',\
                        r'$\beta$', r'$\gamma$FIX', r'$\alpha$', r'$\delta$FIX'])

    path = '../../runs/simpleModelPosteriors/'
    for pName, pLabel in zip(parameterNames, parameterLabels):
        plotOneDPosterior(path,5,pName,pLabel)


    plotHDistribution('../../runs/simpleModelPosteriors/combined/hposterior.dat')





if __name__ == "__main__":
    main()
