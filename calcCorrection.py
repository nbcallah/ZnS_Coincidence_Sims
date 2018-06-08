#!/usr/local/bin/python3

import math
import numpy as np
import sys

class dip:
	def __init__(self, kappa, length, fraction):
		self.kappa = kappa
		self.length = length
		self.fraction = fraction
		
	def contribAlpha(self, totalN):
		return (self.fraction*totalN/self.kappa)*self.kappa*(1-math.exp(-self.length/self.kappa))
	
	def contribBeta(self, totalN):
		return (self.fraction*totalN/self.kappa)**2*(self.kappa/2.0)*(1-math.exp(-2.0*self.length/self.kappa))

def lt(deltaT, alpha, beta, tau, dips, totalN):
	#oneDip = [dip(kappa, length, n0)]
	alphaSumS = 0.0
	alphaSumL = 0.0
	betaSumS = 0.0
	betaSumL = 0.0
	for d in dips:
		alphaSumS += d.contribAlpha(totalN)
		alphaSumL += d.contribAlpha(totalN)*math.exp(-deltaT/tau)
		betaSumS += d.contribBeta(totalN)
		betaSumL += d.contribBeta(totalN)*math.exp(-2.0*deltaT/tau)
	tau_val = deltaT/(math.log((alpha[0]*alphaSumS + beta[0]*betaSumS)/(alpha[0]*alphaSumL + beta[0]*betaSumL)))
	
	tau_err = 0.0
	tau_err += (alpha[1]*(deltaT*beta[0]*(alphaSumL*betaSumS - alphaSumS*betaSumL)) / ((alpha[0]*alphaSumL + beta[0]*betaSumL) * (alpha[0]*alphaSumS + beta[0]*betaSumS) * math.log((alpha[0]*alphaSumS + beta[0]*betaSumS) / (alpha[0]*alphaSumL + beta[0]*betaSumL))**2))**2
	tau_err += (beta[1]*(deltaT*alpha[0]*(alphaSumL*betaSumS - alphaSumS*betaSumL)) / ((alpha[0]*alphaSumL+beta[0]*betaSumL) * (alpha[0]*alphaSumS+beta[0]*betaSumS) * math.log((alpha[0]*alphaSumS + beta[0]*betaSumS) / (alpha[0]*alphaSumL + beta[0]*betaSumL))**2))**2
	tau_err = math.sqrt(tau_err)
	return (tau_val, tau_err)

def alphabeta(points):
	sumInvErr = sum([1/(e**2) for (x,y,e) in points])
	sumXInvErr = sum([x/(e**2) for (x,y,e) in points])
	sumXXInvErr = sum([x**2/(e**2) for (x,y,e) in points])
	sumYInvErr = sum([y/(e**2) for (x,y,e) in points])
	sumXYInvErr = sum([x*y/(e**2) for (x,y,e) in points])
	#points is [(x,y,yerr), ...]
	delta = np.linalg.det(np.matrix([[sumInvErr, sumXInvErr],[sumXInvErr, sumXXInvErr]]))
	alpha = (1/delta)*np.linalg.det(np.matrix([[sumYInvErr, sumXInvErr],[sumXYInvErr, sumXXInvErr]]))
	beta = (1/delta)*np.linalg.det(np.matrix([[sumInvErr, sumYInvErr],[sumXInvErr, sumXYInvErr]]))
	alphaErr = math.sqrt((1/delta)*sumXXInvErr)
	betaErr = math.sqrt((1/delta)*sumInvErr)
	
	return (alpha, alphaErr), (beta, betaErr)

data = np.loadtxt(sys.stdin, dtype=np.float32)

if(data.shape[1] != 3):
    sys.exit("Error! Malformed input")

alpha, beta = alphabeta(data)
oneDip = [dip(7.8, 100, 1.0)]

ltMeas = lt(1020, alpha, beta, 877.6, oneDip, 25000)
correction = (ltMeas[0]-877.6, ltMeas[1])
print(correction[0], correction[1])
