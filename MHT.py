import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import trange

import random
import tests as t
import distns


class MHT:
    def __init__(self,distn='multnorm',param=-0.1,n=100,seed=80085,repl=1000,shift_range=(0,3),resol_shift_prop = 0.01):
        self.n = n
        self.p = param
        self.seed = seed
        self.replicates = repl
        self.shift_range = shift_range
        self.resol_shift_prop = resol_shift_prop
        random.seed(self.seed)
        if distn.lower() == 'multnorm':
            self.dist = distns.MultNorm(n=self.n,p=self.p)
        elif distn.lower() == 'multnom':
            self.dist = distns.Multinomial()
        else:
            print('The distribution is not in our hands, please try "multnorm" for Multivariate Normal, and "multnom" for Multinomial Distribution.')
            self.dist = None
        
    def test_indiv(self,dataa):
        pass

    def generate_type1error_power(self,alltests=True,tests=None, shift=1,alpha=0.05,pi_0=0.5):
        data,p_vals = self.dist.generate_p_vals(pi_0=pi_0,shift=shift)
        #print(p_vals)
        num_of_null = [p_vals[i][2] for i in range(len(p_vals))].count(0)
        num_of_alt = len(p_vals)-num_of_null
        #print('Data generated')
        if alltests:    
            tests = [t.bonferroni_test,t.holms_test,t.sidak_su_test,t.bh_test,t.bh_test_lk_correction,t.bh_test_additive_err_bound_for_neg_dep,t.bh_test_multiplicative_err_bound_for_neg_dep]
            seq = ['bonferroni','holms','sidak','bh','bh_log_corrected','bh_additive_error_bound','bh_mult_error_bound']
            all_rej = [tst(p_vals=p_vals,sigf_val=alpha) for tst in tests]
        
        elif tests:
            seq = [test.__name__ for test in tests]
            all_rej = [tst(p_vals=p_vals,sigf_val=alpha) for tst in tests]
        
        else:
            print('Rejection Function not chosen')

        power = []
        type_1_error = []
        for rej in all_rej:
            t_rej = 0
            f_rej = 0
            for r in rej:
                if p_vals[r-1][2]:
                    t_rej+=1
                else:
                    f_rej+=1
            if num_of_alt:
                pow__ = t_rej/num_of_alt
            else:
                pow__ = 0
            power.append(pow__)
            if num_of_null:
                t1e__ = f_rej/num_of_null
            else:
                t1e__ = 0
            type_1_error.append(t1e__)
        #print('one testing done')
        return power,type_1_error,seq

    
    def test_power_for_alpha_pi(self,alpha=0.05,pi_0=0.5):
        shifts = [(self.shift_range[0]+(self.shift_range[1]-self.shift_range[0])*self.resol_shift_prop*i) for i in range(round(1/self.resol_shift_prop))]
        powers_to_df = []
        type_1errors_to_df = []
        for ii in trange(len(shifts)):
            shift = shifts[ii]
            first_time = True # for setting seq_final and powers_to_df and type_...
            for i in range(self.replicates):
                pow_,typ1err_,seq = self.generate_type1error_power(shift = shift,alpha=alpha,pi_0=pi_0)
                if first_time:
                    seq_final = seq
                    power_add = [0 for i in range(len(pow_))]
                    type_1error_add = [0 for i in range(len(typ1err_))]
                    first_time = False
                power_add = [power_add[i]+pow_[i] for i in range(len(pow_))] # adding for average 
                type_1error_add = [type_1error_add[i]+typ1err_[i] for i in range(len(typ1err_))]
            power = [power_add[i]/self.replicates for i in range(len(power_add))]
            type_1error = [type_1error_add[i]/self.replicates for i in range(len(power_add))]
            #print(power)
            
            powers_to_df.append(power)
            type_1errors_to_df.append(type_1error)
                
        powers = pd.DataFrame(powers_to_df,columns=seq_final)
        type_1errors = pd.DataFrame(type_1errors_to_df,columns=seq_final)
        powers.index = shifts
        type_1errors.index = shifts
        colors = ['r','g','b','k','m','c','y']
        labels = seq_final
                
        for i in range(len(seq_final)):
            plt.plot(powers.index,powers[seq_final[i]],color = colors[i],label = labels[i])
        plt.xlabel(f'Shifts in alternate hyp parameter : {self.dist.__class__.__name__},p = {self.p}')
        plt.ylabel('Power')
        plt.legend()
        #plt.show()
        plt.savefig(f'Plots/Set/testing-power-alpha-{alpha}-pi_0-{pi_0}-seed-{self.seed}-{self.dist.__class__.__name__}-{self.p}-cov_mat1.png')
        plt.clf()
        for i in range(len(seq_final)):
            plt.plot(type_1errors.index,type_1errors[seq_final[i]],color = colors[i],label = labels[i])
        plt.xlabel(f'Shifts in alternate hyp parameter : {self.dist.__class__.__name__},p = {self.p}')
        plt.ylabel('Type 1 error')
        plt.legend()
        #plt.show()
        plt.savefig(f'Plots/Set/testing-type_1_error-alpha-{alpha}-pi_0-{pi_0}-seed-{self.seed}-{self.dist.__class__.__name__}-{self.p}-cov_mat1.png')
        plt.clf()

    def test_power(self,directory='Plots'):
        param=[-0.5,-0.2,0,0.2,0.5]
        alphas = [0.01,0.05,0.1,0.2,0.37]
        pi = [1,0.8,0.6,0.4,0.2,0]
        for p in param:
            self.p = p
            for alpha in alphas:
                for pi_0 in pi:
                    self.test_power_for_alpha_pi(alpha=alpha,pi_0=pi_0)


        

            
        
            