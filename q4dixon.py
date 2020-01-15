import numpy as np
import time
import sys
from sage.all import *

'''
	Idea fox dixon is to generate (1+epsilon)r  Congruences 
	such that in each X_i**2 is B smooth.
	Then we use Sage's linear system solver mod N, which is not that
	efficient but does the job

'''

def get_factor_base(B):
	'''
		Returns a list with primes < B
		denoted the factor base fb
	'''
	P = Primes()
	pr = P.first()
	fb = []
	while(pr <= B):
		fb.append(pr)
		pr = P.next(pr)
	return fb

def smooth_factorization_trdiv(X, fb):
	'''
		Not optimal but works!
		Computes by repeated division the factorization of X into primes
		in the factor base, if X has a prime greater than largest prime in
		factor base return -1

		Input x is in Z/NZ, fb is a list
	'''
	comp = ZZ(X**2)
	pp = [0]*len(fb) # prime powers
	for i in range(len(fb)):
		j = 0
		pr = fb[i]
		while(comp % pr == 0):
			j = j + 1
			comp = comp / pr 
		pp[i] = j
		if comp == 1:
			return pp
	if comp > 1:
		return -1 # its not B smooth
	return pp
def smooth_factorization_sieved(X, fb):
	'''
		Faster method to check if B smooth and get
		factor base factorization
	'''
	comp = ZZ(X**2)
	P = Primes()
	pp = [0]*len(fb)
	F = list(factor(comp))
	if F[-1][0] > fb[-1]:
		return -1
	for pri, power in F:
		pp[P.rank(pri)] = power
	return pp
	
def get_X(Xarr, eps):
	'''
		Computes  product of X_i ^ {epsilon_i} (mod N)
		Input: Xarr, when each element is already in Z/NZ
	'''
 	prod = 1
 	for i in range(len(Xarr)):
 		if eps[i] == 1:
 			prod = prod * Xarr[i]
 	return prod

def get_Y(fb, alpha, N):
	'''
		Computes product of p_j ^ alpha_j (mod N)
		
		fb elements are not in Z/NZ so they need to be coerced
	'''
	prod = 1
	fb = vector(Integers(N),fb)
	for i in range(len(alpha)):
 		prod = prod * fb[i]**alpha[i]
 	return prod

def coerce_vec(ele):
	'''	
		explicitly coerces sage vector
		resolves some errors that come with automatic coercion in sage
	'''
	eZ = []
	for g in ele:
		eZ.append(int(g))
	return np.array(eZ)

def getFac(N, B, frac = 0.1):
	'''
		Dixons Algorithm
	'''
	fb = get_factor_base(B)
	r = len(fb) #size of factor basis
	R = ceil((1 + frac) * r)# the number of congruences we require
	ZnZ = Integers(N)
	Xarr = []
	delmat = []
	runs = 0
	for i in range(R):
		sf = -1
		while(sf == -1): # keep iterating till we get sqaure free
			X = ZnZ.random_element()
			sf = smooth_factorization_trdiv(X, fb)
		Xarr.append(X)
		delmat.append(sf)

	M = Matrix(Integers(2), delmat) # solving the congruences, finding epsilon_i
	for ele in M.kernel().basis():
		eZ = coerce_vec(ele)
		dij_tr = np.array(delmat).transpose()
		alpha = np.matmul(dij_tr, eZ) / 2
		Xf = get_X(Xarr, eZ)
		Yf = get_Y(fb, alpha, N)
		if Xf != Yf and Xf != -Yf:
			return gcd(ZZ(Xf-Yf), N), gcd(ZZ(Xf+Yf),N) #return the non trivial factors
	return "Fail"

def B_heuristic(N):
	'''
		A heuristic for choosing B, so that congruene generation is fast, 
		However in this q4dixon.py we just use a smaller value as Sage's Linear system of congruences solver is slow
	'''
	a = log(N*1.0)
	b = log(a)
	return floor(exp((a*b)**0.5))

def main():
	#val = 8591966237
	val = 2251802665812493
	#val = 73786976659910426999
	B_practical = 5000
	start_time = time.time()
	print getFac(val, B_practical)
	print "Runnning time: " + str(time.time() - start_time) + "s"

if __name__ == "__main__":
    main()