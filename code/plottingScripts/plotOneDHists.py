from matplotlib import pylab
import numpy as np
from math import ceil
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 15})


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







def plotOneDPosterior(data,paramLabel,setBins,colour='k-',label='',prior=False):

    #pathToRuns = "../../runs/simpleModelPosteriors/"
    

    #setBins = np.arange(int(min(data[paramName])), ceil(max(data[paramName])),0.5)
    #setBins = np.linspace(min(data[paramName]), max(data[paramName]),50)
    #setBins = np.linspace(-14.5,2.5,50)

    p05 = np.percentile(data,5)
    p50 = np.percentile(data,50)
    p95 = np.percentile(data,95)

    if prior==True: 
        plt.hist(data, bins=setBins, histtype='step', facecolor=colour, alpha=0.3, edgecolor=None,density=1,fill=True,label=label)

    else: 
        (bins,dat) = histOutline(data,setBins,density=1)
        pylab.plot(bins,dat,colour,linewidth=2,alpha=1.0,label=label)
        plt.axvline(p05,ls=':',color=colour)
        plt.axvline(p95,ls=':',color=colour)
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
    parameterLabels = ([r'$\log_{10}\frac{{\dot n}_0}{{\rm Mpc}^{3}{\rm Gyr}}$',\
                        r'$\beta$', r'$\gamma$FIX', r'$\alpha$', r'$\delta$FIX'])

    orange = '#E69F00'
    blue = '#0072B2'
    green = '#009E73'
    red = '#D55E00'
    purple = '#CC79A7'

    colSimp=orange
    colGalExt=green
    colGal = 'k'
    colSimpNegAlpha = red

    """
    # one d histograms of three runs
    pathToRuns = '../../runs/simpleModel/logNormLike/'
    nRuns=5
    simpleModelData = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)
    plotOneDPosterior(simpleModelData['logn'],
                      parameterLabels[0],
                      colour=colSimp,
                      label='Simple model')
    """
    """
    # negative alpha run
    pathToRuns = '../../runs/simpleModel/logNormLikeNegativeAlpha/'
    nRuns=5
    simpleModelData = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)
    plotOneDPosterior(simpleModelData['logn'],
                      parameterLabels[0],
                      colour=colSimpNegAlpha,
                      label=r'Simple model $\alpha<0$')

    """



    # this one!
    pathToRuns = '../../runs/simpleModel/logNormLikeMstar6to10/'
    nRuns=5
    simpleModelData = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)
    low  = int(min(simpleModelData['logn']))
    high = ceil(max(simpleModelData['logn']))
    bins = np.linspace(low,high,80)
    print(bins)
    plotOneDPosterior(simpleModelData['logn'],
                      parameterLabels[0],
                      setBins=bins,
                      colour=colSimp,
                      label=r'Midd16 posterior')

    
    """ 
    # negative alpha run
    pathToRuns = '../../runs/simpleModel/logNormLikeNegAlphaMstar6to10/'
    nRuns=5
    simpleModelData = readPosterior.readPosterior(pathToRuns, nRuns=nRuns)
    plotOneDPosterior(simpleModelData['logn'],
                      parameterLabels[0],
                      colour=purple,
                      label=r'Simple model $\alpha<0$ and $6<\log_{10}\mathcal{M}_*/M_{\odot}<10$')

    """

    pathToRuns = '../../runs/galaxyModel_ext/n_eff.dat'
    galaxyModelExt = np.genfromtxt(pathToRuns)
    plotOneDPosterior(galaxyModelExt,
                      parameterLabels[0],
                      setBins=bins,
                      colour=colGalExt,
                      label=r'Chen19 posterior')


    
    pathToRuns = '../../runs/galaxyModel_ext/n_effprior.dat'
    galaxyModelExt = np.genfromtxt(pathToRuns)
    plotOneDPosterior(galaxyModelExt,
                      parameterLabels[0],
                      setBins=bins,
                      colour=colGalExt,
                      label=r'Chen19 prior',
                      prior=True)

    



    """
    # prior 
    pathToRuns = '../../runs/galaxyModel/n_effprior.dat'
    galaxyModel = np.genfromtxt(pathToRuns)
    plotOneDPosterior(galaxyModel,
                      parameterLabels[0],
                      colour=colGal,
                      label='Galaxy model prior',
                      prior=True)


    pathToRuns = '../../runs/galaxyModel/n_eff.dat'
    galaxyModel = np.genfromtxt(pathToRuns)
    plotOneDPosterior(galaxyModel,
                      parameterLabels[0],
                      colour=colGal,
                      label='Galaxy model')


    """


    
    plt.legend(fontsize=12)
    plt.xlim(-12.,1.5)
    plt.savefig('combinedAnalysisPlots/lognComparison.pdf',dpi=300)
    plt.show()

    """
    print(simpleModelData)
    exit()
    # not used at the moment 
    path = '../'
    for pName, pLabel in zip(parameterNames, parameterLabels):
        plotOneDPosterior(path,5,pName,pLabel)
        break


    plotHDistribution('./hposterior.dat')

    """



if __name__ == "__main__":
    main()
