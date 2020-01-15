import sys
from sage.all import *
import time
def get_product(budget, prime_bound):
	'''
		Computes X = product of primes^floor(log budget/ log prime)
	'''
	X = 1
	P = Primes()
	p = P.first()
	while(p <= prime_bound):
		X = X * (p**floor(ZZ(budget).log(p)))
		p = P.next(p)
	return X

def multiply(tri1,tri2):
	'''
		Multiplication of two group elements in G_d(N)
		Follows the rule (a1,b1) (a2,b2) = (a1a2 + db1b2, a1b2+a2b1)
	'''
	d = tri1[0]
	a1 = tri1[1]
	b1 = tri1[2]
	a2 = tri2[1]
	b2 = tri2[2]
	return (d, a1*a2 + d*b1*b2, a1*b2 + a2*b1)

def get_random_Gd(N):
	'''
		Get random group element, First choose 'a' and 'b' randomly such that
		b has an inverse mod N, so that 'd' = (a**2 - 1) / b**2  (mod N) can be
		calculated.
	'''
	r = Integers(N)
	a = r.random_element()
	while(True):
		b = r.random_element()
		if(gcd(int(b),N) == 1):
			break
	d = (a**2 -1) / (b**2)
	return (d,a,b)

def brute_mult(X,tri):
	'''
		Brute force multiplication for X * alpha, not recommneded
	'''
	gen_prod = (tri[0],1,0)
	for i in range(X):
		gen_prod = multiply(gen_prod,tri)
	return gen_prod

def sq_mult(X, tri):
	'''
		square and multiply algorithm for computing X * alpha
	'''
	binarr = X.digits(2)
	d = tri[0]
	a = tri[1]
	b = tri[2]
	gen_prod = (d,1,0)
	if binarr[0] == 1:
		gen_prod = tri
	temp_sq = tri
	for i in range(1,len(binarr)):
		temp_sq = multiply(temp_sq, temp_sq)
		if binarr[i] == 1:
			gen_prod = multiply(temp_sq, gen_prod)
	return gen_prod

def pplusone_alg(budget, prime_bound, N):
	'''
		repeat till we get non trivial factor of N
	'''
	while(True):
		re = get_random_Gd(N)
		#print(re)
		X = get_product(budget, prime_bound)
		(d,u,v) = sq_mult(X, re)
		g = gcd([int(u-1), int(v), int(N)])
		if 1 < g < N:
			print "A non trivial factor is:", g
			break

def main():
	#Ex1
	# Pt = Primes()
	# pb = Pt.unrank(1000)
	# bg = 40000000000
	# budget = prime_bound = 2100
	# N = 95853544864250299111409

	#Ex2
	N = 746482824012238308661619491135773503333385064366762059957618554835738449567418578817253229
	Pt = Primes()
	prime_bound = Pt.unrank(1000)
	budget = 40000000000
	start_time = time.time()
	pplusone_alg(budget, prime_bound, N)
	print "Running time: " + str(time.time() - start_time) + "s"

if __name__ == "__main__":
    main()