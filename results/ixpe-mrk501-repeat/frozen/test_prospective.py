"""Synthetic gates, covariance propagation, and nuisance minimization checks."""
import numpy as np
from prospective import evaluate, vector_and_cov

GATES=['input_lineage_disjoint','gtis_valid','response_valid','regions_valid','centroid_rule_valid','spectra_within_family','nuisance_bounds_valid','statistical_approximation_valid']
cfg={'gates':{k:True for k in GATES},'difference_systematic_cov':[[0.,0.],[0.,0.]],'background_multiplier_bounds':[0.,2.],'rmf_recovery_bounds':[0.9636,0.9958],'threshold_chisq':9.210340372}
def summary(q,u):
    d={'source_sum':[10000,10000*q,10000*u],'background_sum':[0,0,0],'source_outer':[[10000,0,0],[0,10000,0],[0,0,10000]],'background_outer':[[0,0,0],[0,0,0],[0,0,0]],'scale':0.1}
    return {'du_sufficient':[d]}
def both(q,u):return {'nominal':summary(q,u),'xy_centroid':summary(q,u)}
a=evaluate(both(0,0),both(0,0),cfg);assert not a['all_alternatives_exceed']
assert a['alternatives']['nominal']['minimum_chisq'] < 1e-10
b=evaluate(both(0,0),both(.2,.2),cfg);assert b['all_alternatives_exceed']
assert b['alternatives']['nominal']['minimum_chisq'] > cfg['threshold_chisq']
rng=np.random.default_rng(772)
lower=b['alternatives']['nominal']['certified_lower_chisq']
for _ in range(1000):
    ab,ah=rng.uniform(0,2,2);fb,fh=rng.uniform(.9636,.9958,2)
    v1,c1=vector_and_cov(both(0,0)['nominal'],ab);v2,c2=vector_and_cov(both(.2,.2)['nominal'],ah)
    delta=v2/fh-v1/fb; value=delta@np.linalg.solve(c1/fb**2+c2/fh**2,delta)
    assert lower <= value + 1e-9
cfg['gates']['response_valid']=False
try:evaluate(both(0,0),both(.2,.2),cfg)
except ValueError as e:assert 'response_valid' in str(e)
else:raise AssertionError('failed to gate unvalidated response')
print('synthetic null, separated alternative, gate pass')
