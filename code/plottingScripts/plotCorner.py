import corner
import numpy as np
import matplotlib.pyplot as plt

import readPosterior

import matplotlib
matplotlib.rcParams.update({'font.size': 14})


def fiveDimSimpleModel():
    # check 
    pathToRuns = '../../runs/simpleModel/logNormLikeNegativeAlpha/'
    data = readPosterior.readPosterior(pathToRuns, nRuns=5)

    dataMstar = [ np.log10(10.0**(delta+7.0)) for delta in data['delta'] ]


    d2 =  np.atleast_1d([data['logn'],
                         data['beta'],
                         data['gamma'],
                         data['alpha'],
                         dataMstar]).T
    parameterLabels = ([r'$\log_{10}\frac{{\dot n}_0}{{\rm Mpc}^{-3}{\rm Gyr}}$',\
                        r'$\beta$', r'$z_0$', r'$\alpha$', r'$\log_{10} \frac{\mathcal{M}_*}{ M_{\odot}}$'])

    figure = corner.corner(d2,labels=parameterLabels,
                           quantiles=[0.05, 0.5, 0.95])
    plt.savefig('{}/plots/5DimCorner.png'.format(pathToRuns))
    plt.savefig('{}/plots/5DimCorner.pdf'.format(pathToRuns))
   
    plt.show()
    return 




def twoDimGalModel():

    galExt = np.genfromtxt('../../runs/galaxyModel_ext/chain_1.txt')
    print(galExt)
    galParNames = np.loadtxt('../../runs/galaxyModel_ext/pars.txt',dtype=str)
    
    tauIndex = list(galParNames).index('tau0')
    betaTauIndex = list(galParNames).index('betatau')
    tau = galExt[:,tauIndex]
    betaTau = galExt[:,betaTauIndex]
    data = np.atleast_1d([tau,betaTau]).T
    parameterLabels = ([r'$\tau_0$',r'$\beta_{\tau}$'])
    figure = corner.corner(data,labels=parameterLabels)
    plt.savefig('gal_ext_2d_corner_tmp.pdf')
    plt.show()

    return 



def main():
    
    twoDimGalModel()



if __name__ == "__main__":
    main()
