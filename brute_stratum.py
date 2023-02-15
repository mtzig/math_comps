from random import shuffle
from surface_dynamics.all import Origami
from collections import Counter



def get_strat_perm(stratum):
    '''
    given stratum, generates corresonding commutator perm
    '''

    p_com = {}
    num = 1
    for cycle in stratum:
        start = num
        for i in range(cycle):
            p_com[num] = num+1
            num+=1
        p_com[num] = start
        num+=1


    return p_com

def random_perm(perm_size):
    vals = [i+1 for i in range(perm_size)]
    shuffle_vals = vals.copy()
    shuffle(shuffle_vals)

    perm = dict(zip(vals, shuffle_vals))

    perm_i = dict(zip(shuffle_vals, vals))

    return perm, perm_i

def one_cycle_perm(perm_size):
    '''
    returns perm of one size
    '''

    perm = {}
    perm_i = {}

    for i in range(1, perm_size):
        perm[i] = i+1
        perm_i[i+1] = i
    perm[perm_size] = 1
    perm_i[1] = perm_size
    return perm, perm_i
                

def find_strat_sts(stratum):

    strat_counts = Counter(stratum)

    perm_size = sum(stratum) + len(stratum)
    counter = 0

    h, hi = one_cycle_perm(perm_size)

    while True:
        counter += 1
        # h, hi = random_perm(perm_size)
        v, vi = random_perm(perm_size)

        c = get_commutator(h,hi,v,vi)
        c_strat = get_stratum(c)

        if Counter(c_strat) == strat_counts:
            print(f'It took {counter} guesses')
            return h, v
    

def get_commutator(h, hi, v, vi):
    '''
    Given h,v and its inverse

    create the commutator subgroup
    '''

    n = len(h)
    p = {}
    for i in range(1,n+1):
        p[i] = v[h[vi[hi[i]]]]

    return p

def get_stratum(p):
    '''
    returns the "stratum" of a permutation
    i.e. we are assuming p is commutator subgroup
    '''

    n = len(p)

    visited = set()

    stratum = []

    for i in range(1,n+1):
        if i in visited:
            continue

        j = i
        cycle_length = 0
        while j not in visited:
            visited.add(j)
            cycle_length += 1
            j = p[j]
        
        stratum.append(cycle_length - 1) # by def of stratum

    return stratum


def convert_cycle(p):
    '''
    returns permutation p
    in cycle notation as string that can be converted to origami obj
    '''

    n = len(p)

    visited = set()

    perm = ''

    for i in range(1,n+1):
        if i in visited:
            continue

        j = i
        cycle = []
        while j not in visited:
            visited.add(j)
            cycle.append(str(j))
            j = p[j]
        
        perm += f'({",".join(cycle)})' # by def of stratum

    return perm
