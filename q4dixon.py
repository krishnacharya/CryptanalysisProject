import numpy as np
import time
from tqdm import tqdm
import sys
from sage.all import *

def get_factor_base(B):
	P = Primes()
	pr = P.first()
	fb = []
	while(pr <= B):
		fb.append(pr)
		pr = P.next(pr)
	return fb

def smooth_factorization_trdiv(X, fb):
	'''
		X is originally is Z/nZ
		Early abort doesnt help much the last example
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
		# if comp > fb[-1] and comp in Primes():
		# 	return -1
	if comp > 1:
		return -1 # its not B smooth
	return pp
def smooth_factorization_sieved(X, fb):
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
 	prod = 1
 	for i in range(len(Xarr)):
 		if eps[i] == 1:
 			prod = prod * Xarr[i]
 	return prod

def get_Y(fb, alpha, N):
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

def getFac(N, B):
	fb = get_factor_base(B)
	#print(fb)
	r = len(fb)
	frac = 0.1
	R = ceil((1 + frac) * r) + 1 
	ZnZ = Integers(N)
	Xarr = []
	delmat = []
	runs = 0
	for i in tqdm(range(R)):
		sf = -1
		while(sf == -1):
			X = ZnZ.random_element()
			sf = smooth_factorization_trdiv(X, fb)
		Xarr.append(X)
		delmat.append(sf)
		# epsilon = lsolve_ZnZ(transpose-delij, mod 2)
		# alphas  = getalphas(1/2epsilon)
		#Xf = prod of X**epsiloni mod N
		#Y = prod p_j**\alpha_j mod N
	# print "A non trivial factor is" gcd(Xf+Y,N)
	print "Started" 
	M = Matrix(Integers(2), delmat) # in sage this is right kernel for Mtranspose 
	for ele in M.kernel():
		eZ = coerce_vec(ele)
		dij_tr = np.array(delmat).transpose()
		alpha = np.matmul(dij_tr, eZ) / 2
		Xf = get_X(Xarr, eZ)
		Yf = get_Y(fb, alpha, N)
		if Xf != Yf and Xf != -Yf:
			return gcd(ZZ(Xf-Yf), N), gcd(ZZ(Xf+Yf),N)	
# fb = get_factor_base(20)
	return "Fail"

def B_heuristic(N):
	a = log(N*1.0)
	b = log(a)
	return floor(exp((a*b)**0.5))
#val = 8591966237
val = 2251802665812493
#val = 73786976659910426999
B_practical = 1000 
start_time = time.time()
#print len(get_factor_base(B_heuristic(val)))
print getFac(val, B_practical)
print time.time() - start_time
# X = r.random_element()
# print smooth_factorization(X, fb)