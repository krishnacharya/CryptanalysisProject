import sys
from sage.all import *

def f(x):
    return x**2 + 1

def pollardrho(N, B):
    '''
        N is the number to factor
        B is the budget
    '''
    R = Integers(N)
    y = x = R.random_element()
    b=1
    while b < B:
        x = f(x)
        y = f(f(y))
        d = gcd(int(x-y), N)
        if d > 1:
            return x-y
        b = b + 1
    return "Fail"

val1 = 60331193824455101058028269521753
val2 = 276474933387964773460419532857385928669681
diff = pollardrho(val2, 10**17)
print diff
print gcd(val2, diff)