"""Split a later visit on a fixed UTC date, estimate event phase on its first segment, and test it on its second."""
from pathlib import Path
import csv
import sys
import numpy as np
from locked.image import get,reduce
from locked.score import rot,depth
VISITS = {
    "29489_73": (29489,73,2460293.5,12.35327485665821,2.00740716,3.1679),
    "3177_22": (3177,22,2458912.5,18.559880330703862,2.31998504,3.547245),
}
def clean(rows,p,epoch):
    flux=rows[:,8]
    ok=(flux>0)&np.isfinite(flux)&(rows[:,7]<np.quantile(rows[:,7],.97))
    t=(rows[ok,0]-epoch)*24
    y=-2.5*np.log10(flux[ok])
    y-=np.array([np.median(y[abs(t-ti)<=15]) for ti in t])
    r=rot(t,y,p)
    ok=abs(r-np.median(r)) < 4*max(1e-8,1.4826*np.median(abs(r-np.median(r))))
    return t[ok],rot(t[ok],y[ok],p)
def phase_score(t,r,p,w,c):
    ph=abs((t/p-c+.5)%1-.5)
    return np.mean(r[ph<w/(2*p)])-np.mean(r[ph>.27]) if np.sum(ph<w/(2*p))>=5 else -np.inf
def analyze(key):
    obj,s,split,p,w,spin=VISITS[key]
    frames=get('stamps'+key)
    training=[row for row in frames if row[0]<split]
    test=[row for row in frames if row[0]>split+.5]
    train,_,_,_,_=reduce(training)
    if len(train)<40:return [obj,s,'insufficient training',len(train)]
    t,r=clean(train,spin,split);cs=np.arange(100)/100
    c=cs[np.argmax([phase_score(t,r,p,w,v) for v in cs])]
    dat,_,_,_,_=reduce(test)
    if len(dat)<40:return [obj,s,'insufficient testing',len(dat)]
    tt,rr=clean(dat,spin,split);ev=depth(tt,rr,p,w,c)
    if not len(ev):return [obj,s,'no covered test cycles',len(dat)]
    d=ev[:,1];lo,hi=np.quantile(np.median(np.random.default_rng(3839).choice(d,(10000,len(d))),axis=1),[.025,.975])
    return [obj,s,len(t),len(tt),split,c,len(d),sum(d>0),round(np.median(d),5),round(lo,5),round(hi,5),round(1.4826*np.median(abs(rr-np.median(rr))),5)]
if __name__=='__main__':
    for key in sys.argv[1:] or VISITS:
        print(analyze(key),flush=True)
