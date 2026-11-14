"""Recover three eclipse centers per TESS visit with local baselines and overlap masks."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
SECTORS=(8,9,34,35,36,61,62,63,88,89,90)
CLOCKS={'A primary':(1.8035,1544.5094,.014),'A secondary':(1.8035,1544.5094+1.8035*.5042,.015),'B primary':(4.3829,1561.2557,.007)}
# Trial phase centers under the rounded 4-decimal literature periods. The
# offsets are recovered by fitting each visit and do not constrain the solution.
def centers(t):
    year=0 if t<1800 else 1 if t<2500 else 2 if t<3300 else 3
    return [ [.0,.002,-.01,-.029][year],[-.009,-.006,-.018,-.038][year],[-.001,0,.009,.02][year] ]

def single_visit(data, sector, aperture, half_range):
    t=data.time.values
    counts=data[aperture].to_numpy()
    baseline=pd.Series(counts,index=pd.to_timedelta(t,unit='D')).rolling('1.1D',center=True,min_periods=5).median().to_numpy()
    y=counts/baseline
    valid=np.isfinite(y)&(abs(y-1)<.55)&(baseline>0)
    t,y=t[valid],y[valid]
    cadence=np.median(np.diff(t))
    approximate=centers(np.median(t))
    results=[]
    for j,(name,(per,epoch,width)) in enumerate(CLOCKS.items()):
        mask=np.ones(len(t),bool)
        for z,(op,oe,ow) in enumerate(CLOCKS.values()):
            if j!=z:
                distance=abs(((t-oe)/op-approximate[z]+.5)%1-.5)
                mask &= distance > max(ow*4,.03)
        first=int(np.floor((t.min()-epoch)/per))
        last=int(np.ceil((t.max()-epoch)/per))
        for cycle in range(first,last+1):
            expected=epoch+cycle*per
            center=expected+approximate[j]*per
            sel=mask&(abs(t-center)<half_range)
            x=t[sel]-center;f=y[sel]
            if min(sum(abs(x)<.035),sum(x<-.06),sum(x>.06))<2:continue
            def residual(a):
                const,slope,depth,delta=a
                sub=x[:,None]+cadence*np.array([-.4,-.2,0,.2,.4])-delta
                prof=np.exp(-.5*(sub/(width*per))**2).mean(axis=1)
                return const+slope*x-depth*prof-f
            fit=least_squares(residual,[1,0,.08,0],bounds=([.65,-1.5,0,-.075],[1.35,1.5,.7,.075]),loss='soft_l1',f_scale=.015)
            depth,shift=fit.x[2:]
            if depth>.02 and abs(shift)<.074:
                results.append(dict(sector=sector,clock=name,aperture=aperture,window=half_range,cycle=cycle,nominal_btjd=expected,residual_min=(center-expected+shift)*1440,depth=depth))
    return results

def run():
    rows=[];coverage=[]
    for s in SECTORS:
        data=pd.read_csv(HERE/'data'/f's{s}.csv.gz')
        eligible=data.box.median()>100
        coverage.append(dict(sector=s,cadences=len(data),median_box=data.box.median(),eligible=eligible))
        if eligible:
            for ap in ['box','central']:
                for hw in [.17,.13]:
                    rows+=single_visit(data,s,ap,hw)
    events=pd.DataFrame(rows)
    events.to_csv(HERE/'events.csv.gz',index=False,compression={'method':'gzip','mtime':0},float_format='%.10f')
    pd.DataFrame(coverage).to_csv(HERE/'coverage.csv',index=False)
    visits=events.groupby(['clock','aperture','window','sector']).agg(median_btjd=('nominal_btjd','median'),dt=('residual_min','median'),events=('residual_min','count')).reset_index()
    visits.to_csv(HERE/'visits.csv',index=False)
if __name__=='__main__':run()
