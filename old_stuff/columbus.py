from brute_stratum import *
from sts_curves import *
from math import gcd
import heapq
import re

class STS:

    def __init__(self, h, v):
        '''
        h and v are permutations, represented as a dictionary
        s.t. the key is the input and the value is the output

        e.g. (1,2)(3)(4,5,6) is
        {1:2, 2:1, 3:3, 4:5, 5:6, 6:4}
        '''
        self.h = h
        self.v = v
        self.h_i = STS.get_inverse(h)
        self.v_i = STS.get_inverse(v) 

    def __str__(self):
        return f'h: {STS.cycle_print(self.h)}\n'+ \
               f'v: {STS.cycle_print(self.v)})'
    
    def get_commutator(self):
        '''
        returns the commutator
        '''
        p = {}
        for i in self.h:
            p[i] = self.v[self.h[self.v_i[self.h_i[i]]]]

        return p
    
    def get_singularities(self):
        '''
        returns the singularies as equivalence classes (i.e. tuples)
        '''
        com = self.get_commutator()
        sings = STS.get_cycles(com)
        
        equivs = []
        for sing in sings:
            equivs.append(tuple(sing))


        return equivs    

    @staticmethod
    def get_inverse(p: dict)-> dict:
        '''
        returns the inverse of perm
        '''
        p_i = {}
        for key in p:
            p_i[p[key]] = key
        
        return p_i
    
    @staticmethod
    def get_cycles(p: dict) -> list:
        '''
        returns permutation p
        as a string in cycle notation
        '''

        n = len(p)

        visited = set()
        cycles = []

        for i in range(1,n+1):
            if i in visited:
                continue

            j = i
            cycle = []
            while j not in visited:
                visited.add(j)
                cycle.append(j)
                j = p[j]
            
            # print(cycle)
            cycles.append(cycle)

        return cycles

    @staticmethod
    def cycle_print(p: dict) -> str:
        '''
        returns permutation p
        as a string in cycle notation
        '''

        perm = ''
        cycles = STS.get_cycles(p)

        for cycle in cycles:
            perm += f'({",".join(map(str,cycle))})'

        return perm
    
    @staticmethod
    def get_cycles(p: dict) -> list:
        '''
        returns permutation p
        as a list of cycles
        '''

        n = len(p)

        visited = set()
        cycles = []

        for i in range(1,n+1):
            if i in visited:
                continue

            j = i
            cycle = []
            while j not in visited:
                visited.add(j)
                cycle.append(j)
                j = p[j]
            
            # print(cycle)
            cycles.append(cycle)

        return cycles
    
    @staticmethod
    def get_cycle_dict(perm_string, n):
        visited = set()
        
        regex = r'\(.+?\)|".+?"|\w+' 
        cycles = re.findall(regex, perm_string)

        perm = {}

        for cycle in cycles:
            c = cycle[1:-1].split(',')
            l = len(c)
            for i in range(l):
                visited.add(int(c[i]))
                perm[int(c[i])] = int(c[(i+1) % l])

        for i in range(1, n+1):
            if i not in visited:
                perm[i] = i

        return perm



    # create static method that creates dictionaries that map each value to a set

def alg1(s: STS):
    r_0 = min(*get_cycles(s.h), *get_cycles(s.v))

    l_0 = 4 #r_0 ** 2

    S_l0 = compute_S_l0(r_0)

    graph = {}

    for v in S_l0:
        subgraph = alg2(s, v)

        for vertice in subgraph:
            graph[vertice] = graph.setdefault(vertice, set()).union(subgraph[vertice])

    return graph


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

def alg2(s: STS, v: tuple):
    A = computeA(v)
    factors = factorA(A)

    s_ = s

    # initialize alpha to identity
    alpha = {} # a bijection??
    s_sings = s.get_singularities()
    for sing in s_sings:
        alpha[sing] = sing

    for m in reversed(factors):
        if m == 'R':
            s_ = STS(s_.v_i, s_.h)
            alpha = updateAlpha(alpha, s_, 'R')

        else: # then m is T
            s_ = STS(s_.h, compose(s_.v, s_.h_i))
            alpha = updateAlpha(alpha, s_, 'T')

    # compute dictionary mapping singularities of A\dot S to S
    alpha_i = STS.get_inverse(alpha)

    # construct l
    l = []
    for equiv in alpha_i:
        if len(equiv) > 1: # prevent fixed point
            l.extend(equiv)
    
    l = set(l)

    ############################ construct the graph
    graph = {}
    EquivsDict = getEquivsDict(s_.get_singularities())
    lv = (v[0] ** 2 + v[1] ** 2) ** .5 # get length of v
    for i in l:
        j = s_.h[i]
        k = 1
        while j not in l:
            j = s_.h[j]
            k += 1

        # generate the directed graph here with vertices labeled by alpha_i[equiv]
            # use list representaion of graph

        graph.setdefault(alpha_i[EquivsDict[i]], []).append((alpha_i[EquivsDict[j]], k*lv))

    return graph

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

def factorA(A):
    '''
    Factors A into T and R

    Code taken from 
    https://codegolf.stackexchange.com/questions/198012/decomposition-of-a-matrix-in-sl-2-mathbbz
    '''
    a,b,c,d = A[0][0], A[0][1], A[1][0], A[1][1]
    
    def f(a,b,c,d):
        if (a,b,c,d)==(1,0,0,1): 
            return ''
        elif a*c+b*d>0: 
            return 'T'+f(a-c,b-d,c,d)
        else: 
            return 'R'+f(c,d,-a,-b)

    return f(a,b,c,d)

def getEquivsDict(equivs: dict) -> dict:
    '''
    returns a dictionary that matches all values in a equivalence class to
    its equivalence class
    '''

    EquivsDict = {}
    for equiv in equivs:
        for i in equiv:
            EquivsDict[i] = equiv

    return EquivsDict

def updateAlpha(alpha: dict, s_: STS, transformation: str):
    '''
     updates the alpha bijection after one iteration of decomposition
    '''

    if transformation == 'T':
        return alpha
    
    # Otherwise, the transformation is R
    EquivsDict = getEquivsDict(s_.get_singularities())

    for equiv in alpha:

        currEquiv = alpha[equiv]

        i = currEquiv[0] # take one element from equiv class
        j = s_.h[i]

        alpha[equiv] = EquivsDict[j]



    return alpha

# thanks ChatGPT for this code
def shortest_closed_path(graph):
    # Find the shortest path from each vertex to every other vertex
    shortest_paths = {}
    for start_vertex in graph:
        dists = {vertex: float('inf') for vertex in graph}
        dists[start_vertex] = 0

        queue = [(0, start_vertex)]
        while queue:
            (distance, current_vertex) = heapq.heappop(queue)

            # Shortest path to current_vertex found
            if current_vertex in shortest_paths:
                continue

            shortest_paths[current_vertex] = distance
            for neighbor, weight in graph[current_vertex]:
                if neighbor not in shortest_paths:
                    heapq.heappush(queue, (distance + weight, neighbor))

    # Find the shortest closed path (cycle) in the graph
    shortest_cycle = float('inf')
    for start_vertex in graph:
        for neighbor, weight in graph[start_vertex]:
            cycle_distance = shortest_paths[start_vertex] + weight + shortest_paths[neighbor]
            if cycle_distance < shortest_cycle:
                shortest_cycle = cycle_distance

    return shortest_cycle