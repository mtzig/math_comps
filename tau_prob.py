# Note have to directly paste into notebook

from surface_dynamics.all import Origami


class Tau:
    
    def __init__(self, n):


        self.n = n

        self.G = SymmetricGroup(n)
        self.A = AlternatingGroup(n)

        self.fixed_h = PermutationGroupElement([tuple([i for i in range(1, n+1)])])
        self.ident = self.G.identity()
        self.irr_G = self.G.irreducible_characters()


    
    @staticmethod
    def commutator(x,y):
        return x*y*x^-1*y^-1
    
    @staticmethod
    def ct_to_cs(ct):
        '''
        converts cycle type to cycle string
        '''
        i = 1
        cs = ''

        for n in ct:
            cs += '('
            for _ in range(n):
                cs += f'{i},'
                i+= 1
            
            cs = cs[:-1] + ')'
        
        return cs

    def compute_closed_form(self, tau, model='rand', conj_tau=True):
        '''
        model types: rand, fixed, fixed_conj
        '''

        conj_tau_size = self.G.conjugacy_class(tau).cardinality()

        if model == 'rand':
                prob = sum(chi(tau) / chi(self.ident) for chi in self.irr_G) / factorial(self.n)
        elif model == 'fixed':
            prob = sum(complex(abs(chi(self.fixed_h))^2) * conjugate(chi(tau)) / chi(self.ident) for chi in self.irr_G) / factorial(self.n)
        else: # fixed cong
            prob = sum(complex(abs(chi(self.fixed_h))^2) * chi(tau) / chi(self.ident) for chi in self.irr_G) / factorial(self.n)

        if conj_tau:
            return N(prob * conj_tau_size)
        else:
            return N(prob)


    def compute_brute(self, model='rand'):
        
        ctype_success = {}

        if model=='rand':
            ctype_success = self.__brute_rand()
        elif model == 'fixed':
            ctype_success = self.__brute_fixed()
        else:
            ctype_success = self.__brute_fixed_conj()

        
  

        return ctype_success

    def __brute_rand(self):
        '''
        computes probability for random model (small n)
        '''
        ctype_success = {}

        for h in self.G:
            for v in self.G:

                com = Tau.commutator(h,v).cycle_type()

                ctype_success[com] = ctype_success.get(com,0) + 1
        for com in ctype_success.keys():    
            ctype_success[com] = N(ctype_success[com]/(self.G.cardinality()^2))

        return ctype_success

    def __brute_fixed(self):
        
        ctype_success = {}

        for v in self.G:
            com = Tau.commutator(self.fixed_h,v).cycle_type()

            ctype_success[com] = ctype_success.get(com,0) + 1

        for com in ctype_success.keys():    
            ctype_success[com] = N(ctype_success[com]/self.G.cardinality())

        return ctype_success

    def __brute_fixed_conj(self):

        ctype_success = {}

        for h in self.G.conjugacy_class(self.fixed_h):
            for v in self.G:
                com = Tau.commutator(h,v).cycle_type()

                ctype_success[com] = ctype_success.get(com,0) + 1
       


        for com in ctype_success.keys():    
            ctype_success[com] = N(ctype_success[com]/(self.G.cardinality()*self.G.conjugacy_class(self.fixed_h).cardinality()))

        return ctype_success

    def gen_sts(self, model='rand'):
        if model == 'rand':
            h = self.G.random_element().cycle_string()
        elif model=='fixed':
            h = self.fixed_h.cycle_string()
        else:
            h = self.G.conjugacy_class(self.fixed_h).random_element()

        v = self.G.random_element().cycle_string()

        return Origami(h, v)




