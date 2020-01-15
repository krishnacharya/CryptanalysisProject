import sys
from sage.all import *
import time
'''
	Idea:
	If the factorization of N is known, dlog in the original group
	reduces to solving dlog in subgroups of smaller order p_i^{e^i}
	and then recombining using CRT
'''

def pohlig_hellman(h,g,N):
	'''
		g^x = h
		group order is N => g^N = identity
	'''
	bp =list(factor(N)) #factors must be known
	y  = []
	m  = []
	for ele in bp:
		m1 = ele[0]**ele[1]
		nexp = (N / m1)
		g1 = g**nexp
		h1 = h**nexp
		y.append(brute_force(h1,g1,m1)) #we can even use shank_steps here
		m.append(m1)
	return CRT_list(y,m)

def brute_force(h,g,m):
	'''
		Find the smallest power i such that g^i = h
		Brute force also works well if the factors of N are small
	'''
	val = g**0
	i = 0
	while(i < m):
		if val == h:
			return i
		i = i+1
		val = val*g
	return "Error, possibly incorrect generator"

def shank_steps(h,g,m):
	'''
		An optimization:
		We can be smarter and instead of brute force use shanks baby step, giant step.
	'''
	return bsgs(g,h,(0,m))

def main():	
	p = 13827821670227353601 #A prime number
	N = p-1 # order of the group Z/pZ
	r = Integers(p)
	h = r(10780909174164501009) # element whose dlog we wish to find
	g = r(3) # the generator for the group
	start_time = time.time()
	print "The dlog is " + str(pohlig_hellman(h,g,N))
	print "Running time: " + str(time.time() - start_time) + "s"

if __name__ == "__main__":
    main()