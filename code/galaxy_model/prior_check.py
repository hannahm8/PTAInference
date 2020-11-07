import numpy
from numpy import inf

def fpairtest(delta):
	return 1./(delta+1.)*(1.-0.25**(delta+1.))

def schechterf(M,M0,phi0,alpha0):
	return phi0*numpy.power(M/M0,1.+alpha0)*numpy.exp(-M/M0)

def Mbh(Mstar,alpha,beta):
	return alpha+beta*numpy.log10(Mstar/1.e11)


def check_p(livepoint):

	data = numpy.genfromtxt('mbulge.dat')
	data0 = numpy.genfromtxt('numdenrange0.txt')
	data1 = numpy.genfromtxt('numdenrange1.txt')
	m0 = numpy.logspace(9,12,15)
	z0 = [0.4,0.6,0.8]
	z1 = [1.5,2.,2.5]

	l = 0.0
	for i in z0:
		phi0 = livepoint['Phi0']+livepoint['PhiI']*i
		a0 = livepoint['alpha']+livepoint['alphaI']*i
		model = schechterf(m0,10.**livepoint['M0'],10.**phi0,a0)
		if (model<=data0[:,6]).all()==True and (model>=data0[:,1]).all()==True:
			for j in z1:
				phi0 = livepoint['Phi0']+livepoint['PhiI']*j
				a0 = livepoint['alpha']+livepoint['alphaI']*j
				model = schechterf(m0,10.**livepoint['M0'],10.**phi0,a0)
				if (model<=data1[:,6]).all()==True and (model>=data1[:,1]).all()==True:
					l += 0.0
				else:
					return -inf
		else:
			return -inf

	model = Mbh(m0,livepoint['Ms'],livepoint['theta'])
	if (model<=data[:,6]).all()==True and (model>=data[:,1]).all()==True:
		return l
	else:
		return -inf


