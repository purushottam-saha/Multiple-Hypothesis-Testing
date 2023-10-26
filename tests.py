import math

def global_merging_test(p_value,sigf_val):
    if p_value<=sigf_val:
        return 1
    else:
        return 0

def single_step_procedure(p_vals,sigf_val):
    rejections = []
    for p_val in p_vals:
        if p_val[1]<=sigf_val:
            rejections.append(p_val[0])
    return rejections

def stepdown_procedure(p_vals,sigf_vals):
    n = len(p_vals)
    accepted = False
    rejections = []
    p_vals_ = sorted(p_vals,key=lambda x:x[1])
    for i in range(n):
        if p_vals_[i][1]<=sigf_vals[i] and (not accepted):
            rejections.append(p_vals_[i][0])
        else:
            accepted = True
    return rejections

def stepup_procedure(p_vals,sigf_vals):
    n = len(p_vals)
    rejected = False
    rejections = []
    p_vals_ = sorted(p_vals,key=lambda x:x[1])
    for i in range(n-1,-1,-1):
        if p_vals_[i][1]>=sigf_vals[i] and (not rejected):
            pass
        else:
            rejections.append(p_vals_[i][0])
    return rejections

def global_simes_merge(p_vals): # [(1,p1,e1),(2,p2,e2),...], ei=0 implies ith pvalue is from null
    p_val = sorted(p_vals,key=lambda x: x[1])
    return min([p_val[i]/(i+1) for i in range(len(p_vals))])*len(p_vals)

def global_simes_test(p_vals,sigf_val):
    return global_merging_test(global_simes_merge(p_vals),sigf_val)

def global_sidak_test(p_vals,sigf_val):
    n = len(p_vals)
    sigf_val = 1-math.pow((1-sigf_val),1/n)
    return single_step_procedure(p_vals=p_vals,sigf_val=sigf_val)

def bonferroni_test(p_vals,sigf_val):
    return single_step_procedure(p_vals=p_vals,sigf_val=sigf_val/len(p_vals))

def holms_test(p_vals,sigf_val):
    n = len(p_vals)
    sigf_vals = [sigf_val/(n-i) for i in range(n)]
    return stepdown_procedure(p_vals=p_vals,sigf_vals=sigf_vals)

def bh_test(p_vals,sigf_val):
    n = len(p_vals)
    sigf_vals = [(i+1)*sigf_val/n for i in range(n)]
    return stepdown_procedure(p_vals=p_vals,sigf_vals=sigf_vals)

def sidak_su_test(p_vals,sigf_val):
    n = len(p_vals)
    sigf_vals = [(1-math.pow((1-sigf_val),1/(i+1))) for i in range(n)]
    return stepup_procedure(p_vals=p_vals,sigf_vals=sigf_vals)

def bh_test_lk_correction(p_vals,sigf_val):
    lk = sum([1/(1+i) for i in range(len(p_vals))])
    sigf_val = sigf_val/lk
    return bh_test(p_vals=p_vals,sigf_val=sigf_val)

def bh_test_additive_err_bound_for_neg_dep(p_vals,sigf_val):
    lk = sum([1/(1+i) for i in range(len(p_vals))])
    sigf_val = sigf_val/min(lk,-math.log(sigf_val)+3.18)
    return bh_test(p_vals=p_vals,sigf_val=sigf_val)

def bh_test_multiplicative_err_bound_for_neg_dep(p_vals,sigf_val):
    lk = sum([1/(1+i) for i in range(len(p_vals))])
    sigf_val = sigf_val/min(lk,-math.log(sigf_val)+3.18)
    return bh_test(p_vals=p_vals,sigf_val=sigf_val)

def k_induced_test(p_vals,sigf_val):
    pass

def a_b_optimization_test(p_vals,sigf_val):
    pass



