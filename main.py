from MHT import MHT
import tests as t
import distns as d

test = MHT(n=100,repl=200,shift_range=(0,5),param=-0.5,resol_shift_prop=0.05)
#test.test_power_for_alpha_pi(alpha=0.1,pi_0=0.3)
#test.generate_type1error_power(shift=3)
#print(t.bh_test_additive_err_bound_for_neg_dep([(1,0.03,0),(2,0.06,1),(3,0.001,0)],0.05))

test.test_power()