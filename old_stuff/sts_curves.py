from brute_stratum import *
from math import lcm

def compose(f,g):
    '''
    returns the composition of f and g i.e. f \circ g
    '''

    h = {}

    for key in f:
        h[key] = f[g[key]]

    return h

def fixed_point(p):

    for key in p:
        if key == p[key]: # exists fixed point
            return True 
    
    return False


def get_cycles(p):
    '''
    returns the cycles of perm p

    (modified from get_stratum)
    '''

    n = len(p)

    visited = set()

    cycles = []

    for i in range(1,n+1):
        if i in visited:
            continue

        j = i
        cycle_length = 0
        while j not in visited:
            visited.add(j)
            cycle_length += 1
            j = p[j]
        
        cycles.append(cycle_length) 

    return cycles

def unit_saddle(hi, c):
    '''
    NOTE: Currently only checks in horizontal direction
    hi is inverse of horizontal perm
    c is commutator

    basically checks if hi c^k has fixed point
    '''

    max_k = lcm(*get_cycles(c))

    comp = hi

    for _ in range(max_k):
        comp = compose(comp, c)
        if fixed_point(comp):
            return True

    return False

    # lcm of cycles

def get_saddle_dist(num_squares, fixed=False, num_samples=1000):
    '''
    get empirical unit saddle prob

    NOTE: currently only checks horizontal
    '''

    saddle_count = 0

    if fixed:
        h, hi = one_cycle_perm(num_squares)

    for _ in range(num_samples):


        if not fixed:
            h, hi = random_perm(num_squares)
        
        v, vi = random_perm(num_squares)

        c =  get_commutator(h, hi, v, vi)

        if unit_saddle(hi, c):
            saddle_count += 1

    return saddle_count