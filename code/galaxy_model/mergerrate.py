import numpy as np
import numpy.random as nr
import scipy.optimize as so
import scipy.integrate as si
import matplotlib.pyplot as plt

#cosmological functions begin

def invE(z):
	OmegaM = 0.3
	Omegak = 0.
	OmegaLambda = 0.7
	return 1./np.sqrt(OmegaM*(1.+z)**3.+Omegak*(1.+z)**2.+OmegaLambda)

def dtdz(z):
	t0 = 14.
	if z == -1.:
		z = -0.99
	return t0/(1.+z)*invE(z)

def dVcdz(z):
	c = 3.e8
	H0 = 70.e3
	return 4.*np.pi*c/H0*invE(z)*DM(z)*DM(z)

def DM(z):
	c = 3.e8
	H0 = 70.e3
	return c/H0*si.quad(invE, 0., z)[0]

def mchirpq(q):
	return np.log10(q**0.6/(1.+q)**0.2)

def mchirp(m1,m2):
	return (m1*m2)**0.6/(m1+m2)**0.2

def mfraction(m):
	if m > 10.:
		return np.sqrt(6.9)*np.exp(-6.9/(2.*(m-10.)))/(m-10.)**1.5+0.615
	else:
		return 0.615

def mfractionp(m):
	if m > 10.:
		return np.exp(-6.9/(2.*(m-10.)))/(m-10.)**3.5/10.**m*(3.94-1.7*(m-10.))
	else:
		return 0.

#cosmological functions end
#spectrum computation functions begin

Gu = 6.673848e-11
cu = 2.997925e8
Msun = 1.98855e30
pc = 3.0856776e16
Gyr = 3.15576e16
G = Gu*Msun/pc**3.0
c = cu/pc

def sig(m):
	return 261.5*(m/1.0e9)**(1.0/4.38)*1000.0/pc

def ms(m):
	return 3.22e3*m**0.862

def at(ms,gamma):
	return 2.95*(ms/1.0e6)**0.596*(2.0**(1.0/(3.0-gamma))-1.0)/0.7549

def ft(H,rho,mc,e):
	return 0.56/(np.pi*G**0.1)*(5.0/64.0*c**5.0*H*rho/sig(2.5*mc)/fe(e))**0.3*mc**(-0.4)

def rhoi(gamma,m):
	return (3.0-gamma)*ms(m)/(4.0*np.pi*at(ms(m),gamma)**3.0)*(2.0*m/ms(m))**(gamma/(gamma-3.0))

def xm(e0):
	return 1293.0/181.0*(e0**(12.0/19.0)/(1.0-e0**2.0)*(1.0+121.0*e0**2.0/304.0)**(870.0/2299.0))**1.5

def g(f):
	return 1.0e-15*(2.913*(f*1.0e8)**0.254*np.exp(-0.807*f*1.0e8)+74.24*(f*1.0e8)**1.77*np.exp(-3.7*f*1.0e8)+4.8*(f*1.0e8)**(-0.676)*np.exp(-0.6/(f*1.0e8)))

def inpl(e,f,f0):
	return g(f*(1.0e-10/f0)*(xm(0.9)/xm(e)))*(1.0e-10*xm(0.9)/f0/xm(e))**(2.0/3.0)

def nm(mc,e0,f,ga,H):
	return (inpl(e0,f,ft(H,10.**ga*rhoi(1.,2.5*10.0**mc),10.0**mc,e0)))**2.0/585080*(500.0/0.7)**3.0*(10.0**mc/4.166e8)**(5.0/3.0)

def nz(z):
	return (1.02/(1.0+z))**(1.0/3.0)

def fe(e):
	return (1.0+(73.0/24.0)*e**2.0+(37.0/96.0)*e**4.0)/(1.0-e**2.0)**3.5

#spectrum computation functions end


class mergerrate(object):

	def __init__(self, M1, q, zp, f, Phi0, PhiI, M0, alpha, alphaI, f0, beta, gamma, delta, t0, epsilon, zeta, eta, Ms, theta, sigma, e0, rho):
		"""
		M1 mass of primary galaxy
		q mass ratio between galaxies
		zp redshift
		f frequency of the spectrum
		Phi0, PhiI, galaxy stellar mass function renormalisation rate
		M0, scale mass
		alpha, alphaI, galaxy stellar mass function slope
		pair fraction: f0 rate, beta redshift power law, gamma mass power law, delta mass ratio power law
		merger time scale: t0 time scale, epsilon mass power law, zeta redshift power law, eta mass ratio power law
		Ms, theta, sigma M*-MBH relation
		e0 eccentricity of the binaries
		rho galaxy density parameter
		"""
		self.M1 = 10.**M1
		self.M1diff = (M1.max()-M1.min())/(len(M1)-1.)/2.
		self._M1red = np.zeros(len(M1))
		for i in range(len(M1)):
			self._M1red[i] = 10.**M1[i]*mfraction(M1[i])
		self.q = q
		self.qdiff = (q.max()-q.min())/(len(q)-1.)/2.
		self.f = 10.**f
		self.f0 = f0
		self.beta = beta
		self.gamma = gamma
		self.delta = delta
		self.t0 = t0
		self.epsilon = epsilon
		self.zeta = zeta
		self.eta = eta
		self.Ms = 10.**Ms
		self.theta = theta
		self.sigma = sigma
		self.twosigma2 = 2.*sigma*sigma
		self.e0 = e0
		self.rho = rho
		self.zp = zp
		self.zpdiff = (zp.max()-zp.min())/(len(zp)-1.)/2.
		self._f0 = f0/self._fpairtest()
		self._Phi0 = Phi0
		self._PhiI = PhiI
		self._M0 = 10.**M0
		self._alpha = alpha
		self._alphaI = alphaI
		self._beta1 = beta - zeta
		self._gamma1 = delta - eta
		#self.MBH1 = self.MBH(self.M1)
		self.MBH1 = self.MBH(self._M1red)
		self.MBH1diff = (self.MBH1.max()-self.MBH1.min())/(len(self.MBH1)-1.)/2.
		self.M2 = np.zeros((len(M1),len(q)))
		self._M2red = np.zeros((len(M1),len(q)))
		self.MBH2 = np.zeros((len(M1),len(q)))
		self.qBH = np.zeros((len(M1),len(q)))
		self.Mc = np.zeros((len(M1),len(q)))
		for i,j in np.ndindex(len(M1),len(q)):
			self.M2[i,j] = self.M1[i]*q[j]
			self._M2red[i,j] = self.M1[i]*q[j]*mfraction(np.log10(self.M1[i]*q[j]))
			#self.MBH2[i,j] = self.MBH(self.M2[i,j])
			self.MBH2[i,j] = self.MBH(self._M2red[i,j])
			self.qBH[i,j] = 10.**self.MBH2[i,j]/10.**self.MBH1[i]
			self.Mc[i,j] = self.MBH1[i]+mchirpq(self.qBH[i,j])
		
	def zprime(self,M1,q,zp):
		"""
		redshift condition, need to improve cut at the age of the universe
		"""
		t0 = si.quad(dtdz,0,zp)[0]
		tau0 = self.tau(M1,q,zp)
		if t0+tau0 < 13.:
			result = so.fsolve(lambda z: si.quad(dtdz,0,z)[0] - self.tau(M1,q,z) - t0,0.)[0]
			#result = so.root(lambda z: si.quad(dtdz,zp,z)[0] - self.tau(M1,q,z) - t0,0.,method='lm').x
			#print M1,q,result.x,result.success
		else:
			result = -1.
		return result

	def Phi(self,M1,Phi0,alpha):
		"""
		galaxy stellar mass function
		"""
		return np.log(10.)*Phi0*(M1/self._M0)**(1.+alpha)*np.exp(-M1/self._M0)

	def curlyF(self,M1,q,z):
		"""
		galaxy pair fraction
		"""
		return self._f0*(1.+z)**self.beta*(M1/1.e11)**(self.gamma)*q**(self.delta)

	def _fpairtest(self,qlow=0.25,qhigh=1.):
		"""
		galaxy pair fraction normalization
		"""
		return (qhigh**(self.delta+1.)-qlow**(self.delta+1.))/(self.delta+1.)

	def tau(self,M1,q,z):
		"""
		merger time scale
		"""
		return self.t0*(M1/5.7e10)**(self.epsilon)*(1.+z)**self.zeta*q**self.eta

	def dnGdM(self,M,n2,alpha1):
		return n2*self._M0*(M/self._M0)**alpha1*np.exp(-M/self._M0)

	def dnGdq(self,q):
		return q**self._gamma1

	def dnGdz(self,z):
		return (1.+z)**self._beta1*dtdz(z)

	def dndM1par(self,M1,q,z,n2,alpha1):
		"""
		d3n/dM1dqdz from parameters, missing b^epsilon/a^gamma
		b = 0.4*h^-1, a = 1
		"""
		return n2*self._M0*(M1/self._M0)**alpha1*np.exp(-M1/self._M0)*q**self._gamma1*(1.+z)**self._beta1*dtdz(z)*(0.4/0.7*1.e11/self._M0)**self.epsilon/(1.e11/self._M0)**self.gamma

	def dndMc(self,M1,q,z,n2,alpha1):
		"""
		d3n/dMcdqdz
		"""
		Mred = M1*mfraction(np.log10(M1))
		return self.dndM1par(M1,q,z,n2,alpha1)/self.dMbdMG(M1)/self.dMBHdMG(Mred)*10.**self.MBH(Mred)*np.log(10.)

	def dMBHdMG(self,M):
		"""
		dMBH/dMG
		"""
		return self.Ms*self.theta/1.e11*(M/1.e11)**(self.theta-1.)

	def dMbdMG(self,M):
		"""
		dMbulge/dMG
		"""
		return mfraction(np.log10(M))+M*mfractionp(np.log10(M))

	def dndz(self,M1,q,z,n2,alpha1):
		"""
		dn/dz
		"""
		Mlow = 10.**(np.log10(M1)-self.M1diff)
		Mhigh = 10.**(np.log10(M1)+self.M1diff)
		qlow = q-self.qdiff
		qhigh = q+self.qdiff
		return n2*self._M0*(1.+z)**self._beta1*dtdz(z)*(qhigh**(self._gamma1+1.)-qlow**(self._gamma1+1.))/(self._gamma1+1.)*si.quad(lambda M: (M/self._M0)**alpha1*np.exp(-M/self._M0),Mlow,Mhigh)[0]*(0.4/0.7*1.e11/self._M0)**self.epsilon/(1.e11/self._M0)**self.gamma

	def dNdtdz(self,M1,q,z,n2,alpha1,zp):
		"""
		dN/dtdlog_10 z
		"""
		return self.dndz(M1,q,z,n2,alpha1)/dtdz(zp)*dVcdz(zp)/1.e9*zp*np.log(10.)

	def MBH(self,M):
		"""
		mass of the black hole Mstar-Mbulge relation without scattering
		"""
		return np.log10(self.Ms*(M/1.e11)**self.theta)

	def dndzint(self,M1,q,z,n2,alpha1):
		"""
		number of mergers per Gyr integrated over M1,q,z
		"""
		scale = (0.4/0.7*1.e11/self._M0)**self.epsilon/(1.e11/self._M0)**self.gamma
		numberM = si.quad(self.dnGdM,10.**(np.log10(M1)-self.M1diff),10.**(np.log10(M1)+self.M1diff),args=(n2,alpha1))[0]
		numberq = si.quad(self.dnGdq,q-self.qdiff,q+self.qdiff)[0]
		numberz = si.quad(self.dnGdz,z-self.zpdiff,z+self.zpdiff)[0]
		return numberM*numberq*numberz*scale

	def output(self,function='zprime'):
		"""
		input 3 x 1d array M1,q,z
		output 3d array (M1,q,z) (galaxy mass, galaxy mass ratio, redshift) of values for function
		"""
		output = np.zeros((len(self.M1),len(self.q),len(self.zp)))
		for i,j,k in np.ndindex(len(self.M1),len(self.q),len(self.zp)):
			z = self.zprime(self.M1[i],self.q[j],self.zp[k])
			#print self.zp[k], z
			#print self.tau(self.M1[i],self.q[j],self.zp[k]), self.tau(self.M1[i],self.q[j],z)
			if z <= 0.:
				output[i,j,k] = 1.e-20
			else:
				Phi0 = 10.**(self._Phi0+self._PhiI*z)
				alpha = self._alpha+self._alphaI*z
				alpha1 = alpha + self.gamma - self.epsilon
				n2 = Phi0*self._f0/self._M0/self._M0/self.t0
				alpha2 = alpha1 - 1. - self._gamma1
				if function=='zprime':
					output[i,j,k] = z
				elif function=='fpair':
					output[i,j,k] = self.curlyF(self.M1[i],self.q[j],z)/self.q[j]**self.delta*self._fpairtest()
				elif function=='curlyF':
					output[i,j,k] = self.curlyF(self.M1[i],self.q[j],z)
				elif function=='tau':
					output[i,j,k] = self.tau(self.M1[i],self.q[j],z)
				elif function=='Phi':
					output[i,j,k] = self.Phi(self.M1[i],Phi0,alpha)
				elif function=='dndM1':
					output[i,j,k] = self.Phi(self.M1[i],Phi0,alpha)*self.curlyF(self.M1[i],self.q[j],z)/self.tau(self.M1[i],self.q[j],z)*dtdz(z)*4.*self.M1diff*self.qdiff
				elif function=='dndM1par':
					output[i,j,k] = self.dndM1par(self.M1[i],self.q[j],z,n2,alpha1)*4.*self.M1diff*self.qdiff*np.log(10.)*self.M1[i]
				elif function=='dndMc':
					output[i,j,k] = self.dndMc(self.M1[i],self.q[j],z,n2,alpha1)*4.*self.MBH1diff*self.qdiff
				elif function=='dndz':
					output[i,j,k] = self.dndz(self.M1[i],self.q[j],z,n2,alpha1)
				elif function=='dNdtdz':
					output[i,j,k] = self.dNdtdz(self.M1[i],self.q[j],z,n2,alpha1,self.zp[k])
				elif function=='dndzint':
					output[i,j,k] = self.number(self.M1[i],self.q[j],z,n2,alpha1)
				else:
					raise UserWarning("output function not defined")
		return output

	def grid(self,n0=None,M1=None,M2=None,function='dndMc'):
		"""
		input 3d array n0, 1d array MBH1, 2d array MBH2
		output 3d array (Mcbh,qbh,z) (black hole chirp mass, black hole mass ratio, redshift) of values for function
		"""
		if n0 is None:
			n0 = self.output(function)
		if M1 is None:
			M1 = 10.**self.MBH1
		if M2 is None:
			M2 = 10.**self.MBH2
		#print np.sum(n0)
		Mcbh = np.linspace(5,11,30)
		qbh = np.linspace(0,1,10)
		Mcbhdiff = (Mcbh.max()-Mcbh.min())/(len(Mcbh)-1.)/2.
		qbhdiff = (qbh.max()-qbh.min())/(len(qbh)-1.)/2.
		output = np.zeros((len(Mcbh),len(qbh),len(self.zp)))
		Mc = np.zeros((len(M1),len(M2[0,:])))
		q = np.zeros((len(M1),len(M2[0,:])))
		for i,j in np.ndindex(len(M1),len(M2[0,:])):
			Mc[i,j] = np.log10(mchirp(M1[i],M2[i,j]))
			if M2[i,j] > M1[i]:
				q[i,j] = M1[i]/M2[i,j]
			else:
				q[i,j] = M2[i,j]/M1[i]
		#counter = 0
		for i,j in np.ndindex(len(M1),len(M2[0,:])):
			for i0,j0 in np.ndindex(len(Mcbh),len(qbh)):
				if abs(Mc[i,j]-Mcbh[i0]) < Mcbhdiff and abs(q[i,j]-qbh[j0]) < qbhdiff:
					for k in range(len(self.zp)):
						output[i0,j0,k] += n0[i,j,k]/1.3
						#counter += 1
				else:
					pass
		#print counter
		return output

	def dispersion(self,function='dndMc'):
		"""
		input 3d array n0, 1d array MBH1, 2d array MBH2
		output 3d array (MBH1,MBH2,z) (black hole mass 1, black hole mass 2, redshift) of dispersed values for function
		"""
		n0 = self.output(function)
		M1 = np.linspace(5,11,30)
		qbh = np.linspace(0,1,10)
		M2 = np.zeros((len(M1),len(qbh)))
		for i0,j0 in np.ndindex(len(M1),len(qbh)):
			M2[i0,j0] = np.log10(10.**M1[i0]*qbh[j0])
		output = np.zeros((len(M1),len(qbh),len(self.zp)))
		A = (M1[1]-M1[0])*(qbh[1]-qbh[0])
		for i,j in np.ndindex(len(self.MBH1),len(self.MBH2[0,:])):
			vol = np.zeros((len(M1),len(qbh)))
			for i0,j0 in np.ndindex(len(M1),len(qbh)):
				n1 = np.exp(-(self.MBH1[i]-M1[i0])*(self.MBH1[i]-M1[i0])/self.twosigma2)
				n2 = np.exp(-(self.MBH2[i,j]-M2[i0,j0])*(self.MBH2[i,j]-M2[i0,j0])/self.twosigma2)
				vol[i0,j0] = n1*n2*A/(np.pi*self.twosigma2)
				#print vol[i0,j0]
			for i0,j0 in np.ndindex(len(M1),len(qbh)):
				for k in range(len(self.zp)):
					output[i0,j0,k] += vol[i0,j0]/np.sum(vol)*n0[i,j,k]
			#print np.sum(vol), n0[i,j,0]
		return self.grid(n0=output,M1=10.**M1,M2=10.**M2)

	def realization(self,function='dndMc'):
		n0 = self.output(function)
		M1 = np.zeros(len(self.MBH1))
		for i in range (len(M1)):
			M1[i] = nr.normal(loc=self.MBH1[i],scale=self.sigma)
		M2 = np.zeros(self.MBH2.shape)
		for i,j in np.ndindex(len(self.MBH1),len(self.MBH2[0,:])):
			M2[i,j] = nr.normal(loc=self.MBH2[i,j],scale=self.sigma)
		return self.grid(n0=n0,M1=10.**M1,M2=10.**M2)

	def hmodelt(self,fbin=None):
		ga = self.rho
		H = 15.0
		Mcbh = np.linspace(5,11,30)
		result = np.zeros(len(self.f))
		rate = self.dispersion(function='dndMc')
		out = np.copy(rate)
		for k in range(len(self.f)):
			n0 = self.hdrop(rate,self.f[k],fbin)
			n0 = np.sum(n0, axis=1)
			n = np.zeros((len(n0[:,0]),len(self.zp)))
			for i,j in np.ndindex(len(n0[:,0]),len(self.zp)):
				n[i,j] = n0[i,j]*nm(Mcbh[i],self.e0,self.f[k],ga,H)*nz(self.zp[j])*2.*self.zpdiff
			n = np.sum(n)
			result[k] = 0.5*np.log10(n)
		#print result, out
		return result, out

	def hdrop(self,n0,f,fbin):
		if fbin is None:
			return n0
		c = 3.e8
		G = 1.33e20
		fhigh = f + fbin
		flow = f - fbin
		Mcbh = np.linspace(5,11,30)
		nin = np.zeros((n0.shape))
		for i,j,k in np.ndindex(n0.shape):
			nin[i,j,k] = n0[i,j,k]*dVcdz(self.zp[k])/dtdz(self.zp[k])*5./96./(2.*np.pi)**(8./3.)*c**5./(G*10.**Mcbh[i])**(5./3.)*3./8.*(flow**(-8./3.)-fhigh**(-8./3.))/Gyr*2.*self.zpdiff/(0.5+self.zp[k]/2.)**(11./3.)
		nout = n0
		nm = np.sum(nin,axis=(1,2))
		number = 0
		for cut in range(len(nm)):
			if number < 1:
				number += nm[-cut-1]
			else:
				break
		#print cut,number
		for i,j,k in np.ndindex(n0.shape):
			if i in range(cut-1):
				nout[-i-1,j,k] = 0.

		return nout

