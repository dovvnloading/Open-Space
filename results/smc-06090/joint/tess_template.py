"""Refit two eclipse centers to sector-stratified, phase-folded TESS flux sums."""
import numpy as np
from scipy.optimize import least_squares

PHASE=(np.arange(120)+.5)/120
GROUPS={'2018_2019':[1,2,13],'2020':[27,28],'2023':[67,68],'2025':[94,95,96]}
FIELDS=['primary_phase','secondary_phase','separation','primary_depth_e_s','secondary_depth_e_s','primary_width','secondary_width']

def wrap(q):return (q+.5)%1-.5

def fit(S,N,wide=False):
    # Each sector supplies its native-cadence phase means; pool sectors by covered
    # orbital-cycle count, avoiding overweighting high-cadence sectors.
    A=[];W=[]
    for s,n in zip(S,N):
        ss=s.sum(axis=0);nn=n.sum(axis=0)
        arr=np.divide(ss,nn,out=np.zeros_like(ss),where=nn>0)
        weight=nn/max(nn.mean(),1)*len(s)
        A.append(arr);W.append(weight)
    w=np.sum(W,axis=0); y=np.average(np.array(A),axis=0,weights=np.array(W))
    def model(z):
        c1,c2,d1,d2,b,a,h=z[:7];w1,w2=z[7:] if wide else (.025,.025)
        return b+a*np.cos(2*np.pi*PHASE)+h*np.sin(2*np.pi*PHASE)-d1*np.exp(-.5*(wrap(PHASE-c1)/w1)**2)-d2*np.exp(-.5*(wrap(PHASE-c2)/w2)**2)
    start=[0,.53,1.2,1.2,0,0,0];lo=[-.09,.41,0,0,-10,-10,-10];hi=[.11,.61,25,25,10,10,10]
    if wide:start.extend([.025,.025]);lo.extend([.010,.010]);hi.extend([.060,.060])
    res=least_squares(lambda z:(model(z)-y)*np.sqrt(w/w.mean()),start,bounds=(lo,hi),max_nfev=120)
    c1,c2,d1,d2=res.x[:4];return np.r_[c1,c2,c2-c1,d1,d2,res.x[7:] if wide else [.025,.025]]

