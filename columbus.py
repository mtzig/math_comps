from brute_stratum import *
from sts_curves import *
from math import gcd


def alg1(h,v):
    r_0 = min(*get_cycles(h), *get_cycles(v))

    l_0 = r_0 ** 2

    S_l0 = compute_S_l0(r_0)

def compute_S_l0(r):
    '''
    r is the square root of l0
    
    TODO: might want to only generate a "basis" of vectors
          i.e. don't generate (-i,j), (j,i), (-j,i) 
    '''
    l_0 = r**2
    S_l0 = []

    # we always include (1,0)
    S_l0.append((1,0))
    S_l0.append((0,1))


    for i in range(1,r):
        for j in range(i, r):
            if gcd(i,j) == 1: # relatively prime
                if i**2 + j**2 <= l_0:
                    S_l0.append((i,j))
                    S_l0.append((-i,j))

                    if i != j: # hardcode to ignore (1,1) duplicate
                        S_l0.append((j,i))
                        S_l0.append((-j,i))


    return S_l0

def alg2(v):
    A = computeA(v)

    return 0

def computeA(v):
    v_1 = abs(v[0])
    v_2 = abs(v[1])

    # first handle degenerate cases (we have a one)
    if v_1 == 1:
        d = v_2 + 1
        b = 1
    elif v_2 == 1:
        d = 1
        b = v_1-1
    else:
        for i in range(1, v_2+1):
            if i*v_1 % v_2 == 1:
                d = i
                b = (i*v_1 - 1) // v_2
                break
    
    if v[0] < 0:
        d = -d

    if v[1] < 0:
        b =  -b

    # might want to use numpy here
    A = [[d, -b],[-v[1],v[0]]]

    return A
