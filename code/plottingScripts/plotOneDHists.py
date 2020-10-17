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







def plotOneDPosterior(data,paramLabel,colour='k-',label=''):

    #pathToRuns = "../../runs/simpleModelPosteriors/"
    

    #setBins = np.arange(int(min(data[paramName])), ceil(max(data[paramName])),0.5)
    #setBins = np.linspace(min(data[paramName]), max(data[paramName]),50)
    setBins = np.linspace(-15,3,50)

    p05 = np.percentile(data,5)
    p50 = np.percentile(data,50)
    p95 = np.percentile(data,95)


    (bins,dat) = histOutline(data,setBins,normed=1)
    pylab.plot(bins,dat,colour,linewidth=2,alpha=0.7,label=label)
    plt.axvline(p05,ls=':',color='k')
    plt.axvline(p95,ls=':',color='k')
    plt.xlabel(paramLabel)
    plt.ylabel('Normalised counts')
    plt.tight_layout()
    #plt.savefig('{}lognoHist.pdf'.format(runLoc))
    #plt.savefig('logno.png')
    #pylab.show()


    print("""{}
    05%: {}
    50%: {}
    95%: {}
    """.format(paramLabel,p05,p50,p95))

    return 0




def plotHDistribution(pathToH):

    plt.clf()

    hData = np.genfromtxt(pathToH)
    log10hData = np.log10(hData)

    ph05 = np.percentile(log10hData,5)
    ph50 = np.percentile(log10hData,50)
    ph95 = np.percentile(log10hData,95)
    print ("""
    05%: {}  ({})
    50%: {}  ({})
    95%: {}  ({})
    """.format(10**ph05, ph05, 10**ph50, ph50, 10**ph95, ph95))


    hSetBins = np.arange(int(min(log10hData)),ceil(max(log10hData)),0.05)
    #print(hSetBins)
    (hbins,hdat) = histOutline(log10hData,hSetBins)
    pylab.plot(hbins,hdat,'k-',linewidth=2,alpha=0.7)
    #plt.xlim(-15.5,-14.2)
    #plt.xlim(-15.2,-13.5)
    plt.axvline(ph05,ls=':',color='k')
    plt.axvline(ph95,ls=':',color='k')
    plt.axvline(ph50,ls="--",color='k')
    plt.ylabel('Number of counts')
    plt.xlabel(r'$A_{\rm yr}$')


    # Nanograv numbers
    plt.axvline(np.log10(1.37E-15),color='#ffa500')
    #plt.axvlines(np.log10(1.37E-15),color='#ffa500')
    plt.axvline(np.log10(2.67E-15),color='#ffa500')



    plt.tight_layout()
    #plt.savefig('{}hposterior.png'.format(runLoc))
    #plt.savefig('{}hposterior.pdf'.format(runLoc))
    plt.savefig('hposterior.png')
    #plt.show()
    

    return 0









def main():

    parameterNames = (['logn','beta','gamma','alpha','delta'])
    parameterLabels = ([r'$\log_{10}\frac{{\dot n}_0}{{\rm Mpc}^{-3}{\rm Gyr}}$',\
                        r'$\beta$', r'$\gamma$FIX', r'$\alpha$', r'$\delta$FIX'])


    # one d histograms of three runs
    pathToRuns = '../../runs/simpleModel/logNormLike/'
    nRuns=5
    simpleModelData = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)
    plotOneDPosterior(simpleModelData['logn'],
                      parameterLabels[0],
                      colour='#d95f02',
                      label='Simple model')

    pathToRuns = '../../runs/galaxyModel/n_eff.dat'
    galaxyModel = np.genfromtxt(pathToRuns)
    plotOneDPosterior(galaxyModel,
                      parameterLabels[0],
                      colour='#1b9e77',
                      label='Galaxy model')

    
    pathToRuns = '../../runs/galaxyModel_ext/n_eff.dat'
    galaxyModelExt = np.genfromtxt(pathToRuns)
    plotOneDPosterior(galaxyModelExt,
                      parameterLabels[0],
                      colour='#7570b3',
                      label='Galaxy model ext.')
    
    plt.legend()
    plt.savefig('overplotted_logn.png')
    plt.show()

    
    print(simpleModelData)
    exit()
    # not used at the moment 
    path = '../'
    for pName, pLabel in zip(parameterNames, parameterLabels):
        plotOneDPosterior(path,5,pName,pLabel)
        break


    plotHDistribution('./hposterior.dat')





if __name__ == "__main__":
    main()
