import numpy as np
import mergerrate as mr
import prior_check

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

def read_data(datafile):  
    
    """
    the data file contains the following columns:
    frequency  - log noise amplitude - log hc - SNR - log hc (with error) - standard deviation
    """
    x = np.genfromtxt(datafile,skip_footer=1)

    return x

#names=['Phi0', 'PhiI', 'M0', 'alpha', 'alphaI', 'f0', 'beta', 'gamma', 'delta', 't0', 'epsilon', 'zeta', 'eta', 'Ms', 'theta','sigma','e0','rho']
bounds=[(-3.4,-2.4),(-0.6,0.2),(11,11.5),(-1.5,-1.),(-0.2,0.2),(0.01,0.05),(0.,2.),(-0.5,0.5),(-0.2,0.2),(0.1,10.),(-0.5,0.5),(-3.,1.),(-0.2,0.2),(7.75,8.75),(0.9,1.1),(0.2,0.5),(0.01,0.99),(-2.,2.)]
data=read_data('input.dat')

def in_bounds(par):
    return all(bounds[i][0] < par[i] < bounds[i][1] for i in range(len(par)))

def log_prior(par):
    if in_bounds(par):
        param = dict(Phi0 = par[0], PhiI = par[1], M0 = par[2], alpha = par[3], alphaI = par[4], f0 = par[5], beta = par[6], gamma = par[7], delta = par[8], t0 = par[9], epsilon = par[10], zeta = par[11], eta = par[12], Ms = par[13], theta = par[14], sigma = par[15], e0 = par[16], rho = par[17])        
        return prior_check.check_p(param)
    else:
        return -np.inf

def log_likelihood(par):
    M1 = np.linspace(9,12,25)
    q = np.linspace(0.25,1,10)
    z = np.linspace(0.,1.5,5)
    f = np.log10(data[:,0])
    #fbin = data[0,0]/3.
    initpar = dict(Phi0 = par[0], PhiI = par[1], M0 = par[2], alpha = par[3], alphaI = par[4], f0 = par[5], beta = par[6], gamma = par[7], delta = par[8], t0 = par[9], epsilon = par[10], zeta = par[11], eta = par[12], Ms = par[13], theta = par[14], sigma = par[15], e0 = par[16], rho = par[17])
    model = mr.mergerrate(M1,q,z,f,**initpar).hmodelt(fbin=None)[0]
    l = -0.5/(data[:,5]*data[:,5])*(model-data[:,2])*(model-data[:,2])-0.5*np.log10(2*np.pi) #detection
    return np.sum(l)

x0 = np.array([-2.8,-0.2,11.25,-1.25,0.,0.025,0.8,0.,0.,1.,0.,-0.5,0.,8.25,1.,0.4,0.5,0.])
ndim = len(x0)
cov = np.diag(np.ones(ndim) * 0.01**2)
N = 1000000

sampler = ptmcmc(ndim, log_likelihood, log_prior, cov, outDir='./output', resume=False)
sampler.sample(x0, N, SCAMweight=30, AMweight=15, DEweight=50)
