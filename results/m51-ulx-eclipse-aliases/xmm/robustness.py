"""Stress-test the spatial fit and compare diametrically opposite apertures."""
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from analyze import HERE, NEIGHBOR, read_visit, measure, windows, EPOCH, INGRESS, LATE_EGRESS

def diagnostic(v,k):
    prefix,p,gti=v
    P=(LATE_EGRESS-EPOCH)/k
    D=(1 if k<300 else 2)*P-(INGRESS-EPOCH)
    w=windows(gti,P,D,.08,True)
    if not len(w):return
    t=p.mjd_tdb_geo.to_numpy()
    jj=np.clip(np.searchsorted(w[:,0],t)-1,0,len(w)-1)
    use=(t>=w[jj,0])&(t<=w[jj,1])&p.energy_ev.between(900,7000)
    xy=p[['x_arcsec','y_arcsec']].values[use]
    z=[]
    # Same distance from ULX-2 at four rotations of its direction toward ULX-1.
    for degrees in (0,90,180,270):
        a=np.radians(degrees)
        rot=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        c=NEIGHBOR-rot@NEIGHBOR
        z.append(int(np.sum(np.linalg.norm(xy-c,axis=1)<=3)))
    return dict(k=k,visit=prefix[1:11],detector=prefix[11:13],target=z[0],left=z[1],opposite=z[2],right=z[3],symmetric_p=binomtest(z[0],sum(z),.25,alternative='greater').pvalue)

if __name__=='__main__':
    fixtures={('0824450901',499),('0677980701',499),('0883550201',506),('0883550301',506),('0883550301',253),('0883550301',507),('0883550201',507)}
    summary=[]; controls=[]
    axis=NEIGHBOR/np.linalg.norm(NEIGHBOR)
    grid=[('base',{}),('narrow_psf',dict(rc=3.5,alpha=1.65)),('broad_psf',dict(rc=6,alpha=1.40)),('soft',dict(lo=500,hi=2200)),('hard',dict(lo=1800,hi=7000)),('box15',dict(half_size=15)),('box23',dict(half_size=23)),('shift_toward_neighbor',dict(dx=1.5*axis[0],dy=1.5*axis[1])),('shift_away',dict(dx=-1.5*axis[0],dy=-1.5*axis[1]))]
    for f in sorted((HERE/'inputs').glob('*-photons.csv.gz')):
        v=read_visit(f)
        for obs,k in sorted(fixtures):
            if obs!=v[0][1:11]:continue
            controls.append(diagnostic(v,k))
            for label,opts in grid:
                for row in measure(v,k,**opts):
                    if row['state']=='in':
                        summary.append(dict(variation=label,**row))
    pd.DataFrame(summary).to_csv(HERE/'xmm-sensitivity.csv',index=False,float_format='%.7g')
    pd.DataFrame([r for r in controls if r]).to_csv(HERE/'xmm-azimuth.csv',index=False,float_format='%.7g')
