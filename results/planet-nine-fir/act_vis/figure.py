"""Reproduce the coverage map and blind synthetic-source recovery plot."""
import csv,json
import numpy as np
from scipy.ndimage import uniform_filter,maximum_filter
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from itertools import combinations
from search_patch import sources
from common import CACHE,HERE,N,HALF,MOTION
from dense import score_tile

def main():
    D=np.array([np.load(CACHE/f'data_{i}.npy') for i in range(4)])
    W=np.array([np.load(CACHE/f'weight_{i}.npy') for i in range(4)])
    yy,xx=np.indices((N,N));inside=(abs(xx-HALF)<=600)&(abs(yy-HALF)<=600)
    local=np.array([uniform_filter((w>0).astype('f4'),size=7)>.7 for w in W])
    sep=sources().query(np.column_stack((xx.ravel(),yy.ravel())))[0].reshape(xx.shape)
    field=(local.sum(axis=0)>=3).astype(int)*2
    for i,j in combinations(range(4),2):
        pair=[i,j];sel=(local.sum(axis=0)==2)&local[pair].all(axis=0)
        z=score_tile(D[pair],W[pair],MOTION[pair][None])[0]
        field[sel&(maximum_filter(abs(z),size=25)<8)]=1
    field[(sep<=15)|~inside]=0
    fig,axes=plt.subplots(1,2,figsize=(8,3.65),layout='constrained')
    a,b=axes;clrs=['#e9edf0','#a8cfd1','#276674']
    a.imshow(field,origin='lower',extent=[-62.55,62.55,-62.55,62.55],cmap=ListedColormap(clrs),interpolation='nearest')
    a.scatter(0,0,marker='+',c='#f44b2e',s=95,lw=1.6,label='ACT best fit')
    a.set(xlabel='East offset (arcsec)',ylabel='North offset (arcsec)',xlim=(-60,60),ylim=(-60,60))
    a.legend(handles=[Patch(color=clrs[0],label='Untested'),Patch(color=clrs[1],label='Two frames'),Patch(color=clrs[2],label='Three or four'),a.collections[0]],loc='lower left',fontsize=7)
    tbl=list(csv.DictReader((HERE/'supplementary_injections.csv').open()))
    m=np.array([float(q['VIS_AB']) for q in tbl]);ratio=np.array([float(q['either_fraction']) for q in tbl])
    lo,hi=[np.array([float(q[col]) for q in tbl]) for col in ('wilson95_lo','wilson95_hi')]
    b.fill_between(m,lo,hi,color='#b4d5d6',alpha=.75)
    b.plot(m,ratio,'o-',color='#276674',label='Either test',lw=2,ms=4)
    b.plot(m,[int(q['three_or_four_frames_recovered'])/2000 for q in tbl],'s--',color='#bc6b35',label='Three or four frames',ms=3)
    b.set(xlabel='Injected VIS magnitude (AB)',ylabel='Fraction of 2,000 offsets recovered',ylim=(0,1.015),xlim=(20,27))
    b.grid(alpha=.17);b.legend(fontsize=8)
    fig.savefig(HERE/'act4_06_vis.png',dpi=170)
if __name__=='__main__':main()
