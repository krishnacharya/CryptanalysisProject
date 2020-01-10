import sys
from sage.all import *
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
	val = g**0
	i = 0
	while(i < m):
		if val == h:
			return i
		i = i+1
		val = val*g
	return "Error, possibly incorrect generator"

def shank_steps(h,g,m):
	return bsgs(g,h,(0,m))

p = 13827821670227353601
N = p-1
r = Integers(p)
h = r(10780909174164501009)
g = r(3)

# p = 11251
# N = p-1
# r = Integers(p)
# h = r(9689)
# g = r(23)
print pohlig_hellman(h,g,N)
#print discrete_log_generic(h,g)