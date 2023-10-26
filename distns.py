import random
import numpy as np
from scipy.stats import norm
import math

def make_neg(a):
    if a<=0:
        return a
    return -a


def make_cov_mat_1(n,p):
    cov_mat = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j:
                cov_mat[i][j]=1
            else:
                cov_mat[i][j] = make_neg(math.pow(p,abs(i-j)*2))
    return cov_mat

def make_cov_mat_2(n,p):
    cov_mat = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j:
                cov_mat[i][j]=1
            else:
                cov_mat[i][j] = make_neg(math.pow(p,abs(i-j)))
    return cov_mat

class MultNorm:
    def __init__(self,n=10000,p=-0.1):
        self.n = n
        self.p = p
        self.null = 'Normal(0,1)'
        self.alt = 'Normal(shift,1)' # specify shift in null_to_alt, hence in generate_p_vals
        self.dep = '1 in diag, cov_mat_ij = make_neg(p^abs(i-j)) for i \\ne j, p -ve'
        self.mean = [0 for i in range(n)]
        self.cov_mat = make_cov_mat_1(self.n,self.p)
        
    def generate(self):
        data = np.random.multivariate_normal(mean=self.mean,cov=self.cov_mat)
        data = [(i+1,data[i]) for i in range(self.n)]
        return data
    
    def null_to_alt(self,data,alt_indices,shift=1):
        for alt in range(1,1+len(data)):
            if alt in alt_indices:
                data[alt-1]=(data[alt-1][0],data[alt-1][1]+shift,1)
            else:
                data[alt-1]=(data[alt-1][0],data[alt-1][1],0)
        return data

    def generate_p_vals(self,pi_0=0.5,shift=1): # right sided p_value, h1>0
        data = self.generate()
        n_0 = round(self.n*pi_0)
        alt_ind = random.sample(range(1,self.n+1),self.n-n_0)
        data = self.null_to_alt(data,alt_ind,shift)
        p_vals = [(i+1,1-norm.cdf(data[i][1]),int((i+1) in alt_ind)) for i in range(self.n)] #(i,u,1) means follows alt
        return data,p_vals


class Multinomial:
    def __init__(self):
        self.null = 'Multnom all same prob'
        self.alt = 'Multnom not all same prob'
        self.dep = 'already neg dep'

        def generate(self,n=10000):
            pass

        def null_to_alt(self):
            pass

        def generate_p_vals(self,n=10000,pi_0=0.5):
            data = self.generate(n=n)
            n_0 = n*pi_0
            alt_ind = random.sample(range(1,n+1),n-n_0)
            data = self.dist.null_to_alt(self.dist.generate(),alt_ind)
            p_vals = [(i+1,0,int((i+1) in alt_ind)) for i in range(n)] # binomial cdf something 
            return data,p_vals