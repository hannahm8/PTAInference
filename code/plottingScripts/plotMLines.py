import numpy as np

import matplotlib.pyplot as plt

import confBands 

import matplotlib
matplotlib.rcParams.update({'font.size': 15})



orange = '#E69F00'
blue = '#0072B2'
green = '#009E73'
red = '#D55E00'


green2 = '#004D40'
blue2 = '#1E88E5'
orange2 = '#FFC107'
red2 = '#D81B60'

colSimp=orange
colGalExt=green
colGal = 'k'

simpModRunLoc = '../../runs/simpleModel/logNormLikeMstar6to10/'
simpModData = np.genfromtxt('{}combined/mLines.dat'.format(simpModRunLoc))
simpModMs   = np.genfromtxt('{}combined/ms.dat'.format(simpModRunLoc))

_,p05,p25,p50,p75,p95,_ = confBands.getPercentiles(simpModMs,simpModData)

plt.fill_between(simpModMs,p05,p95,alpha=0.3,color=colSimp)
plt.fill_between(simpModMs,p25,p75,alpha=0.3,color=colSimp)
plt.plot(simpModMs,p50,color=colSimp,ls='--')


galExtModRunLoc = '../../runs/galaxyModel_ext/'
galExtModData = np.genfromtxt('{}/m_lines.dat'.format(galExtModRunLoc))
galExtModLogMs = np.genfromtxt('{}/mc.txt'.format(galExtModRunLoc))
galExtModMs = [ 10.0**logM for logM in galExtModLogMs ]

_,p05,p25,p50,p75,p95,_ = confBands.getPercentiles(galExtModMs,galExtModData)

plt.fill_between(galExtModMs,p05,p95,alpha=0.3,color=colGalExt)
plt.fill_between(galExtModMs,p25,p75,alpha=0.3,color=colGalExt)
plt.plot(galExtModMs,p50,color=colGalExt,ls='--')


galExtModPrior = np.genfromtxt('{}/m_linesprior.dat'.format(galExtModRunLoc))
p00p5,_,_,_,_,_,p99p5 = confBands.getPercentiles(galExtModMs,galExtModPrior)
plt.plot(galExtModMs, p00p5, color='k', alpha=0.7, ls=':')
plt.plot(galExtModMs, p99p5, color='k', alpha=0.7, ls=':')



"""
galModRunLoc = '../../runs/galaxyModel/'
galModData = np.genfromtxt('{}/m_lines.dat'.format(galModRunLoc))
galModLogMs = np.genfromtxt('{}/mc.txt'.format(galModRunLoc))
galModMs = [ 10.0**logM for logM in galModLogMs ]

p05,p25,p50,p75,p95 = confBands.getPercentiles(galModMs,galModData)

plt.fill_between(galModMs,p05,p95,alpha=0.3,color=colGal)
plt.fill_between(galModMs,p25,p75,alpha=0.3,color=colGal)
plt.plot(galModMs,p50,color=colGal,ls='--')
"""








plt.ylim(1E-9,1E6)
plt.xlim(1E6,1E11)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$\mathcal{M}~~({\rm M_{\odot}})$')
plt.ylabel(r'${\rm d}n / {\rm d} \log_{10} \mathcal{M}/{\rm M_{\odot}}~~({\rm Mpc}^{-3})$')
plt.tight_layout()
plt.savefig('combinedAnalysisPlots/dndlogM.pdf'.format(simpModRunLoc))
plt.savefig('combinedAnalysisPlots/dndlogM.png'.format(simpModRunLoc))
plt.show()
