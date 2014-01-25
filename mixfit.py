import sys
import scipy as sp
import numpy as np
from scipy.stats import norm
from scipy import integrate
from scipy.optimize import minimize,fmin_bfgs

# Complementary error function, modeling observation decline near observation limit
def erfc(x):
	res = 2*norm.cdf(-sp.sqrt(2)*x)
	return res

# This is modeling the number count of objects N(x)=10^(x-m0); ignores the existence of an observation limit
# N0 and m0 redundant as adds only multiplicative term, neutralised by normalization in dmag
def dN(x,alpha):
	pow = alpha*x
	res = np.power(10,pow)
	return res
	
# Combining dN and erfc gives a model of the observed number counts
def dmag_erf(x,minit,mlim,sigma,alpha):
	def temp(x):
		n = np.shape(x)
		if n == ():
			if dN(x,alpha) != np.inf:
				res = 0.5*erfc((x-mlim)/sigma)*dN(x,alpha)*(1-0.5*erfc((x-minit)/sigma))
			if dN(x,alpha) == np.inf:
				res = 0
		elif n > 1:
			n = np.shape(x)[0]
			res = np.zeros((n))
			for i in xrange(0,n):
				if dN(x[i],alpha) < np.inf:
					res[i] = 0.5*erfc((x[i]-mlim)/sigma)*dN(x[i],alpha)*(1-0.5*erfc((x[i]-minit)/sigma))
				if dN(x[i],alpha) == np.inf:
					res[i] = 0
		return res
	reltol = (2.220446e-16)**0.25# Not sure how to get machine numerical characteristics with Python
	normc = integrate.quad(temp,-100,100,epsrel=reltol,epsabs=reltol/100)[0]
	if np.shape(x) == ():
		res = temp(x)/normc
		return res
	elif np.shape(x) > 1:
		n = np.shape(x)[0]
		res = np.zeros((n))
		for i in xrange(0,n):
			if temp(x[i]) == 0.0:
				res[i] = 0
			else:
				res[i] = temp(x[i])/normc
	return res
	
# Fits a model of observed number counts to magnitude data
# x = vector of magnitude
# para = np.array((mlim,sigma,alpha,minit)) initial parameters	
def FitMagInit(x,para):
	def NegLog(para):
		y = dmag_erf(x,para[3],para[0],para[1],para[2])
		if np.shape(y) == ():
			if y < 10**(-12):
				y = 10**(-12)
		elif np.shape(y) > 1:
			n = np.shape(y)[0]
			for i in xrange(0,n):
				if y[i] < 10**(-12):
					y[i] = 10**(-12)
				else:
					y[i] = y[i]
		y = -np.sum(np.log(y))
		return y		
	pars = fmin_bfgs(NegLog, x0=para)
	return pars
	
# dat[0:(n-1)] = data vector
# dat[n:(2*n-1)] = binary vector w/ length same as x indicating if a coordinate is (1) or is not (0) missing, i.e. not taken (1) or taken (0)
def SelectiveMean(dat):
	n = np.shape(dat)[0]/2
	x = dat[0:n]
	missidx = dat[n:(2*n)]
	if n==0 | n-np.around(n) != 0 | np.sum(missidx[missidx != 0] != 1) > 0:
		sys.exit("Input data in bad format!")
	res = np.mean(x[missidx==0])
	return res
	
