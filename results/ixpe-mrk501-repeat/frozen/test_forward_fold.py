"""Synthetic forward-model recovery with actual DU1/3 response matrices."""
import numpy as np
from forward_fold import responses,fit

chunks=[]
for d in (1,3):
    for date in ('20260101','20260701'):
        rsp=responses(d,date);e,a,m,r,widths=rsp;v=e**(-2.4)*a*widths
        norm=20.;i=v@r*norm;p=(v*m)@r*norm
        chunks.append({'response':rsp,'obs':np.column_stack((i,0.12*p,-0.20*p)),'whiten':np.array([np.eye(3) for _ in i])})
result=fit(chunks)
print(result)
assert abs(result['gamma']-2.4)<1e-5
assert abs(result['q']-0.12)<1e-5
assert abs(result['u']+0.20)<1e-5
