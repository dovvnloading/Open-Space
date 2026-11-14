"""Read the two hash-verified OGLE I-band light curves, using the original Open-Space follow-up conventions."""
import hashlib
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
EP=2460000.0369
P=1.1696719
TC=2456000
Y=365.25
SCALE=1.8
HASHES={'ogle3.dat':'d57614117ef140adf1413baab42f92182b9c39752eed39706b6fc040e203a897','ogle4.dat':'429caaf74f2086f3a8f53c71a91bfb3c07049dbc2e81122c06646387c8fa21a0'}
def get():
    data=[]
    for i,n in enumerate(HASHES):
        f=HERE/'inputs'/n;assert hashlib.sha256(f.read_bytes()).hexdigest()==HASHES[n]
        a=np.loadtxt(f)[:,:3]
        if i==0:a[:,0]+=2450000
        a=a[np.isfinite(a).all(axis=1)&(a[:,2]>0)&(a[:,2]<=.1)]
        data.append(a)
    h,m,e=np.concatenate(data).T
    return h,(h-TC)/Y,(h-EP)/P%1,m,np.hypot(e,.012)
def circ(q):return (q+.5)%1-.5
