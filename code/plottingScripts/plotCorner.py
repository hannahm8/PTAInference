import corner
import numpy as np
import matplotlib.pyplot as plt

import readPosterior

import matplotlib
matplotlib.rcParams.update({'font.size': 16})


def fiveDimSimpleModel():
    # check 
    pathToRuns = '../../runs/agnosticModel/logNormLikeMstar6to10/'
    data = readPosterior.readPosterior(pathToRuns, nRuns=5)

    dataMstar = [ np.log10(10.0**(delta+7.0)) for delta in data['delta'] ]


    d2 =  np.atleast_1d([data['logn'],
                         data['beta'],
                         data['gamma'],
                         data['alpha'],
                         dataMstar]).T
    parameterLabels = ([r'$\log_{10}\frac{{\dot n}_0}{{\rm Mpc}^{3}{\rm Gyr}}$',\
                        r'$\beta_z$', r'$z_0$', r'$\alpha_{\mathcal{M}}$',\
                        r'$\log_{10} \frac{\mathcal{M}_*}{ M_{\odot}}$'])

    figure = corner.corner(d2,labels=parameterLabels,
                           quantiles=[0.05, 0.5, 0.95])
    plt.savefig('simple_5DimCorner.pdf'.format(pathToRuns))
   
    plt.show()
    return 




def fourDimGalModel():

    galExt = np.genfromtxt('../../runs/galaxyModel_ext/chain_1.txt')
    print(galExt)
    galParNames = np.loadtxt('../../runs/galaxyModel_ext/pars.txt',dtype=str)
    

    mStarIndex = list(galParNames).index('Mstar')
    mStar      = galExt[:,mStarIndex]

    tauIndex = list(galParNames).index('tau0')
    tau      = galExt[:,tauIndex]

    alphaTauIndex = list(galParNames).index('alphatau')
    alphaTau      = galExt[:,alphaTauIndex]

    betaTauIndex  = list(galParNames).index('betatau')
    betaTau       = galExt[:,betaTauIndex]


    data = np.atleast_1d([mStar,tau,alphaTau,betaTau]).T
    parameterLabels = ([r'$\log_{10} M_{*}$',
                        r'$\tau_0$',
                        r'$\alpha_{\tau}$',
                        r'$\beta_{\tau}$'])

    figure = corner.corner(data,labels=parameterLabels)
    plt.savefig('gal_ext_4d_corner_tmp.pdf')
    plt.show()

    return 


def twoDimGalModel():

    galExt = np.genfromtxt('../../runs/galaxyModel_ext/chain_1.txt')
    print(galExt)
    galParNames = np.loadtxt('../../runs/galaxyModel_ext/pars.txt',dtype=str)
    

    mStarIndex = list(galParNames).index('Mstar')
    mStar      = galExt[:,mStarIndex]

    tauIndex = list(galParNames).index('tau0')
    tau      = galExt[:,tauIndex]


    data = np.atleast_1d([mStar,tau]).T
    parameterLabels = ([r'$\log_{10} M_{*}$',
                        r'$\tau_0$'])

    figure = corner.corner(data,labels=parameterLabels)
    plt.savefig('gal_ext_2d_corner_tmp.pdf')
    plt.show()

    return 


def main():
    fiveDimSimpleModel()
    twoDimGalModel()
    fourDimGalModel()


if __name__ == "__main__":
    main()
