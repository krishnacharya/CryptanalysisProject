import sys
from sage.all import *
import time

'''
    Implements Pollards rho for factoring
    Pseducode:
        Choose X_0 randomly in Z/NZ
        Y_0 = X_0
        Run for atmost O(B) number of times:
            X_0 = f(X_0)
            Y_0 = f(f(Y_0))
            if(gcd(X_0-Y_0,N) > 1) 
                implies we have a non trivial factor, return it and end.
'''

def f(x, c = 1):
    '''
        f(x) is of the form x**2 + c, c is generally small valued
        The idea is that this sufficiently randomizes when taken mod N, 
        default value for c set to 1
        Input x, c
        Where x is in Z/NZ
    '''
    return x**2 + c

def pollard_rho(N, B):
    '''
        N is the number to factor
        B is the budget
    '''
    R = Integers(N) #Z/NZ
    x = R.random_element()
    y = x 
    b=1
    while b < B:
        x = f(x)
        y = f(f(y))
        d = gcd(int(x-y), N)
        if d > 1:
            return x-y, d
        b = b + 1
    return "Fail", "IncrB"

def main():
    '''
        Makes the call to pollard_rho
    '''
    val = 60331193824455101058028269521753
    #val = 276474933387964773460419532857385928669681
    budget = 10**17
    start_time = time.time()
    rval, d = pollard_rho(val, budget)
    if rval != "Fail":
        print "A non trivial facor is, " + str(d)
    else:
        print "Try increasing budget"
    print "Running time: " + str(time.time() - start_time) + "s"

if __name__ == "__main__":
    main()