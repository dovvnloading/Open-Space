"""Repeat recovery with limb-darkened transit shapes and verify a full-cadence injection subset."""
import csv
import sys
from pathlib import Path
import numpy as np
from analyze import HERE, KNOWN, filter_matrix, shape


def occultation(dt, impact, ppm):
    """Normalized occulted intensity from 512 stellar annuli, quadratic limb darkening."""
    k = np.sqrt(ppm * 1e-6)
    x = np.linspace(0., 40.5, 811)
    z = np.sqrt(impact**2 + ((1+k)**2-impact**2) * (x/40.5)**2)
    r = (np.arange(512)+.5)/512
    mu = np.sqrt(1-r*r)
    intensity = 1 - .25*(1-mu) - .18*(1-mu)**2
    cosphi = (r[None,:]**2 + z[:,None]**2 - k*k) / (2*np.maximum(z[:,None],1e-12)*r[None,:])
    loss = np.arccos(np.clip(cosphi,-1,1)) @ (r*intensity)
    loss /= loss[0]
    return np.interp(np.abs(dt), x, loss, right=0.)


def run(source=None):
    a = np.genfromtxt(HERE/'one_minute.csv',delimiter=',',names=True)
    bj,y,n=a['bjd_tdb_minus_2400000_5'],a['relative_flux'],a['count']
    t=(bj-bj[0])*1440
    mask=t>t.max()-8
    for mid,duration in KNOWN:
        mask |= abs(bj+2400000.5-mid)*1440 < duration/2+8
    centers = np.arange(0,np.floor(t[-1])+1,2.)
    injected = np.arange(1,np.floor(t[-1]),2.)
    W,valid=filter_matrix(t,n,mask,centers,1,80)
    W = W[valid]
    centers=centers[valid]
    depths=W@(y-1)
    threshold=max(.0005,-depths.min())
    eligible=np.min(abs(centers[:,None]-injected[None,:]),axis=0)<1.1
    with (HERE/'shape_checks.csv').open('w') as f:
        writer=csv.writer(f)
        writer.writerow(['injected_ppm','impact_parameter','covered_recovered','covered_centers','all_recovered','all_centers'])
        for ppm in [1000,1500]:
            for impact in [0.,.7,.95]:
                template=occultation(t[:,None]-injected[None,:],impact,ppm)
                amps=depths[:,None]-ppm*1e-6*(W@(y[:,None]*template))
                best=amps.argmax(axis=0)
                recovered=(amps[best,np.arange(len(best))]>threshold)&(abs(centers[best]-injected)<=15)
                writer.writerow([ppm,impact,int(np.sum(recovered&eligible)),int(eligible.sum()),int(recovered.sum()),len(recovered)])
    if source is None:
        return
    raw=np.loadtxt(line for line in Path(source).open('rb') if line[:1].isdigit())
    rt=(raw[:,1]-bj[0])*1440
    ry=raw[:,2]/np.nanmedian(raw[:,2])
    ids=np.floor((raw[:,1]-raw[0,1])*1440).astype(int)
    starts=np.r_[0,np.flatnonzero(np.diff(ids))+1]
    stops=np.r_[starts[1:],len(ids)]
    full=[]
    for minute in injected[::13]:
        injected_y=ry*(1-.001*shape(rt,minute,81.,11.34))
        measured=np.array([np.median(injected_y[s:e]) for s,e in zip(starts,stops) if e-s>=10])
        actual=W@(measured-1)
        approx=W@(y*(1-.001*shape(t,minute,81.,11.34))-1)
        v1=actual.argmax();v2=approx.argmax()
        full.append([int(minute),round(1e6*(actual[v1]-approx[v2]),3),
                     round(centers[v1]-centers[v2],2),
                     int(actual[v1]>threshold and abs(centers[v1]-minute)<=15),
                     int(approx[v2]>threshold and abs(centers[v2]-minute)<=15)])
    with (HERE/'full_cadence_checks.csv').open('w') as f:
        writer=csv.writer(f)
        writer.writerow(['injected_minute','peak_depth_difference_ppm','peak_center_difference_minutes','raw_recovery','binned_recovery'])
        writer.writerows(full)

if __name__=='__main__':
    run(sys.argv[1] if len(sys.argv)>1 else None)
