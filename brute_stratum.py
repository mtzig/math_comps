from random import shuffle
from surface_dynamics.all import Origami
from collections import Counter
from itertools import permutations


def get_strat_perm(stratum):
    '''
    given stratum, generates corresponding commutator perm
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
    '''
    generates random perm of perm_size
    '''
    vals = [i+1 for i in range(perm_size)]
    shuffle_vals = vals.copy()
    shuffle(shuffle_vals)

    perm = dict(zip(vals, shuffle_vals))

    perm_i = dict(zip(shuffle_vals, vals))

    return perm, perm_i

def one_cycle_perm(perm_size):
    '''
    returns perm of one size
    and its inverse
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
    '''
    finds a pair of permutation for given stratum
    (by random guessing)
    '''

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
    Given h,v and its inverses

    create the commutator permutation
    '''

    n = len(h)
    p = {}
    for i in range(1,n+1):
        p[i] = v[h[vi[hi[i]]]]

    return p

def get_stratum(p):
    '''
    returns the "stratum" of a permutation
    i.e. we are assuming p is commutator
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

def get_genus(p):
    strat = get_stratum(p)

    return int((sum(strat) + 2) / 2)


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

def gen_fixed_sts(num_squares):
    '''
    generates sts of num_squares squares
    fixes h to be a one cycle

    Namemly generates h, hi, v, vi
    '''

    perm = [i for i in range(1, num_squares+1)]

    h = [i % num_squares + 1 for i in range(1,num_squares+1)]
    for v in permutations(perm):
        yield dict(zip(perm, h)), dict(zip(h, perm)), \
                dict(zip(perm, v)), dict(zip(v, perm))


def gen_all_sts(num_squares):
    '''
    generates all sts of num_squares squares

    Namemly generates h, hi, v, vi
    '''

    perm = [i for i in range(1, num_squares+1)]


    for h in permutations(perm):
        for v in permutations(perm):
            yield dict(zip(perm, h)), dict(zip(h, perm)), \
                  dict(zip(perm, v)), dict(zip(v, perm))
    
def gen_cyl_sts(num_squares):
    '''
    assumes an even number squares
    '''
    perm = [i for i in range(1, num_squares+1)]
    perm_half = perm[num_squares/2:]



    h = [i % (num_squares/2) + 1 for i in range(1,num_squares/2+1)] + \
        [i % (num_squares/2) + 1 + num_squares/2 for i in range(1,num_squares/2+1)]

    for v in permutations(perm[:num_squares/2]):
        # print(perm_half)
        yield dict(zip(perm, h)), dict(zip(h, perm)), \
                dict(zip(perm, perm_half + list(v))), dict(zip(perm_half + list(v), perm))



def get_stratums(num_squares, fixed=False, perms=False):
    '''
    Finds the distribution of stratums given we fix n
    '''

    stratums = {}
    
    gen = gen_fixed_sts if fixed else gen_all_sts
    
    # gen =  gen_cyl_sts

    for  sts in gen(num_squares):
        c =  get_commutator(*sts)
        stratum = get_stratum(c)

        # converts stratum to hashable tuple
        stratum_t =  tuple((sorted(Counter(stratum).items())))

        
        if perms:
            _ = stratums.setdefault(stratum_t,[])
            stratums[stratum_t].append((sts[0],sts[2]))
        else:
            stratums[stratum_t] = stratums.setdefault(stratum_t,0) + 1

    return stratums



def get_stratums_sample(num_squares, fixed=False, genus=False, num_samples=1000):
    '''
    Finds distribution of stratums given we fix n
    Uses random sampling to find distriubtion
    '''

    stratums = {}
    
    if fixed:
        h, hi = one_cycle_perm(num_squares)

    for _ in range(num_samples):

        if not fixed:
            h, hi = random_perm(num_squares)
        
        v, vi = random_perm(num_squares)

        c =  get_commutator(h, hi, v, vi)

        # TODO: rename variables
        if genus: # counts genus
            stratum_t = get_genus(c)
        else: # counts stratum
            stratum = get_stratum(c)
            stratum_t =  tuple((sorted(Counter(stratum).items())))

        stratums[stratum_t] = stratums.setdefault(stratum_t,0) + 1

    return stratums

def compare_fixed_unfixed(num_squares, genus=False, num_samples=1000):

    stratums_unfixed = get_stratums_sample(num_squares, fixed=False, genus=genus, num_samples=num_samples) 
    stratums_fixed = get_stratums_sample(num_squares, fixed=True, genus=genus, num_samples=num_samples) 

    # we want to get set of all stratums that were visted
    keys = set(stratums_unfixed.keys()).union(set(stratums_fixed.keys()))

    strat_diff = {}

    for key in keys:
        u_val = stratums_unfixed.setdefault(key, 0)
        f_val = stratums_fixed.setdefault(key, 0)

        diff = abs(u_val - f_val)

        strat_diff[key] = diff

    return stratums_unfixed, stratums_fixed, strat_diff 
