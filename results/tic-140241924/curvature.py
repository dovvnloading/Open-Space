"""Compare four seasonal epochs against a linear ephemeris for each clock."""
from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
SEASON={9:0,34:1,36:1,61:2,62:2,63:2,88:3,89:3}

def diagnostic(t,y):
    """Middle epochs relative to the line between the endpoints, plus blind last prediction."""
    scale=(t-t[0])/(t[-1]-t[0])
    resid=y-y[0]-(y[-1]-y[0])*scale
    # Linear least squares on first three seasons gives a prediction for #4.
    X=np.stack([np.ones(3),scale[:3]],axis=1)
    pred=np.array([1,scale[-1]])@np.linalg.pinv(X)@y[:3]
    return np.array([resid[1],resid[2],y[-1]-pred])

def main():
    data=pd.read_csv(R/'events.csv.gz')
    data['season']=data.sector.map(SEASON)
    rng=np.random.default_rng(140241924)
    rows=[]; visits=[]
    for (ap,hw,clock),g in data.groupby(['aperture','window','clock']):
        times=[];centers=[];sample=[];sample_time=[]
        for season,sub in g.groupby('season'):
            per=[];boot=[];t=[];times_boot=[]
            for s,a in sub.groupby('sector'):
                residual=a.sort_values('cycle').residual_min.to_numpy()
                n=len(residual)
                if not n:continue
                per.append(np.median(residual));t.append(a.nominal_btjd.median())
                starts=rng.integers(0,n,(4000,int(np.ceil(n/3))))
                idx=((starts[:,:,None]+np.arange(3))%n).reshape(4000,-1)[:,:n]
                boot.append(np.median(residual[idx],axis=1))
                times_boot.append(np.median(a.sort_values("cycle").nominal_btjd.to_numpy()[idx],axis=1))
            times.append(np.mean(t));centers.append(np.mean(per))
            # Nested bootstrap of visits plus contiguous blocks of three observed
            # eclipses per visit; the singleton sector 9 cannot measure visit systematics.
            ns=len(boot); pick=rng.integers(0,ns,(ns,4000))
            sample.append(np.mean(np.array(boot)[pick,np.arange(4000)],axis=0))
            sample_time.append(np.mean(np.array(times_boot)[pick,np.arange(4000)],axis=0))
            visits.append(dict(aperture=ap,window=hw,clock=clock,season=season,
                               mean_btjd=np.mean(t),delta=np.mean(per),visits=len(t),events=len(sub)))
        times=np.array(times);centers=np.array(centers);sample=np.array(sample)
        val=diagnostic(times,centers)
        distribution=np.array([diagnostic(np.array(sample_time)[:,i],sample[:,i]) for i in range(sample.shape[1])])
        for i,metric in enumerate(['2021_minus_chord','2023_minus_chord','2025_minus_forecast']):
            low,high=np.percentile(distribution[:,i], [2.5,97.5])
            rows.append(dict(aperture=ap,window=hw,clock=clock,metric=metric,value=val[i],low=low,high=high))
    pd.DataFrame(rows).to_csv(R/'curvature.csv',index=False,float_format='%.6f')
    pd.DataFrame(visits).to_csv(R/'seasons.csv',index=False,float_format='%.8f')
if __name__=='__main__':main()
