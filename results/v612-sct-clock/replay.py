"""Recompute the V612 Sct clocks, positional controls and phase predictions."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
GAIA_ID = "4104113350446549888"
SKY = SkyCoord(277.941193, -14.315539, unit=u.deg)
EPOCH = 2459500.
FREQ_MIN, FREQ_MAX = .3, 22
RNG = np.random.default_rng(94712)


def bjd(time, body):
    # HJD_UTC needs the Sun-barycentre correction. Geocentric MJD_UTC needs Earth's.
    stamp = Time(time, format="jd", scale="utc")
    position = get_body_barycentric(body, stamp)
    correction = (position.dot(SKY.cartesian) / c).to_value(u.day)
    return stamp.tdb.jd + correction


def clean(group):
    group = group[(group.catflags == 0) & (group.magerr < .18) & (group.mag > 12)].copy()
    group = group.sort_values("hjd").drop_duplicates("expid")
    t, y = group.hjd.to_numpy(), group.mag.to_numpy()
    boundaries = zip(np.searchsorted(t, t-25), np.searchsorted(t, t+25))
    group["res"] = y - np.array([np.median(y[a:b]) for a,b in boundaries])
    width = max(.04, np.quantile(np.abs(group.res), .75))
    return group[np.abs(group.res) < 5*width + .1].copy()


def design(t, frequency):
    angle = 2*np.pi*(t-EPOCH)*frequency
    return np.column_stack((np.ones(len(t)), np.cos(angle), np.sin(angle)))


def fit(t, y, w, frequency):
    matrix = design(t, frequency)
    coeff = np.linalg.lstsq(matrix*w[:,None]**.5, y*w**.5, rcond=None)[0]
    return coeff, float(np.sum(w*(y-matrix@coeff)**2))


def refine(t, y, w, f, width):
    return float(minimize_scalar(lambda q: fit(t,y,w,q)[1], bounds=(f-width,f+width),
                                 method="bounded",options={"xatol":1e-12}).x)


def frequency_scan(t, y, w, samples=8):
    ls = LombScargle(t-EPOCH, y, dy=1/np.sqrt(w))
    f, power = ls.autopower(minimum_frequency=FREQ_MIN, maximum_frequency=FREQ_MAX,
                            samples_per_peak=samples, method="fast")
    permitted = abs(f/1.0027379-np.round(f/1.0027379))>.025
    idx=np.argmax(np.where(permitted,power,-np.inf))
    return refine(t,y,w,f[idx],1/np.ptp(t)/samples), len(f)


def phase(coeff):
    return float(np.arctan2(coeff[2],coeff[1])/(2*np.pi))


def circular(x):
    return float((x+.5)%1-.5)


def prediction(t, y, w, f, coeff, invert=False):
    model = design(t,f)@coeff
    observed, _ = fit(t,y,w,f)
    return dict(n=len(t), correlation=float(np.corrcoef(-model if invert else model,y)[0,1]),
                phase_difference_cycles=circular(phase(observed)-phase(coeff)-(.5 if invert else 0)),
                semiamplitude=float(np.hypot(*observed[1:])))


def sector_stats(d, f, coeff):
    t=d.btjd.to_numpy()+2457000; y=d.corrected_e_per_s.to_numpy();w=np.ones(len(t))
    out=prediction(t,y,w,f,coeff,invert=True)
    # In this independent data set the frequency search is allowed to choose any peak.
    ls=LombScargle(t-EPOCH,y)
    freq,power=ls.autopower(minimum_frequency=FREQ_MIN,maximum_frequency=FREQ_MAX,samples_per_peak=10)
    out["free_peak_d_inv"]=float(freq[np.argmax(power)])
    out["strongest_within_0p15_d_inv"]=float(freq[abs(freq-f)<.15][np.argmax(power[abs(freq-f)<.15])])
    out["power_at_ztf"]=float(ls.power(f))
    out["daily_alias_powers"]=[float(ls.power(f+offset)) for offset in (-1.0027379,1.0027379)]
    # Null shifts phase independently by UTC day, retaining within-day covariance and cadence.
    angle=2*np.pi*(t-EPOCH)*f
    labels=np.floor(t).astype(int); unique=np.unique(labels)
    z=np.array([np.sum((y[labels==k]-np.mean(y[labels==k]))*np.exp(-1j*angle[labels==k])) for k in unique])
    # One-sided test of the phase from ZTF, accounting for the magnitude-to-flux sign.
    align=np.exp(2j*np.pi*(phase(coeff)+.5))
    obs=float(np.real(np.sum(z)*align))
    null=np.real((np.exp(2j*np.pi*RNG.random((10000,len(z))))@z)*align)
    out["daily_phase_null"]={"days":len(z),"trials":len(null),"exceedances":int(np.sum(null>=obs))}
    out["without_background_projection"]=prediction(t[d.unprojected_good==1],d.unprojected_e_per_s[d.unprojected_good==1].to_numpy(),w[d.unprojected_good==1],f,coeff,invert=True)
    out["mean_cadence_minutes"]=float(np.median(np.diff(t))*1440)
    return out


def photocenter(rows,t,frequency,coeff):
    """First-moment diagnostic for a variable component inside an unresolved blend."""
    template=(design(t,frequency)@coeff-coeff[0])/np.hypot(*coeff[1:])
    season=np.floor((t-EPOCH)/150)
    offsets=np.column_stack([season==k for k in np.unique(season)])
    model=np.column_stack((offsets,template))
    coords=rows[["ra","dec"]].to_numpy().copy()
    coords-=np.median(coords,axis=0)
    coords*=np.array([np.cos(SKY.dec.to_value(u.rad)),1])*3600
    data=np.column_stack((coords,rows.res))
    def estimate(index):
        signal=np.linalg.lstsq(model[index],data[index],rcond=None)[0][-1]
        return -signal[:2]/(.4*np.log(10)*signal[2])
    months=np.floor((t-EPOCH)/30)
    groups=[np.flatnonzero(months==k) for k in np.unique(months)]
    rng=np.random.default_rng(1537)
    boots=[estimate(np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))]))
           for _ in range(1600)]
    return {"offset_east_north_arcsec":estimate(np.arange(len(rows))).tolist(),
            "marginal_95_bounds_east_north":np.quantile(boots,[.025,.975],axis=0).tolist()}


def main():
    ztf=pd.read_csv(DATA/"ztf-field.csv"); target=ztf[ztf.oid==333114400007906]
    # Original candidate screen used the catalog's 2017.6 eruption year plus 400 days.
    rows=clean(target);rows=rows[rows.hjd>2451544.5+(2017.6-2000)*365.25+400]
    t=bjd(rows.hjd.to_numpy(),"sun");y=rows.res.to_numpy()
    w=1/(rows.magerr.to_numpy()**2+.06**2)
    train=t<np.quantile(t,.65); masks={"early":train,"late":~train,"all":np.ones(len(t),dtype=bool)}
    clocks={}
    for label,mask in masks.items():
        f,grid=frequency_scan(t[mask],y[mask],w[mask])
        v,_=fit(t[mask],y[mask],w[mask],f)
        clocks[label]=dict(n=int(mask.sum()),bjd_range=[float(t[mask].min()),float(t[mask].max())],
                           period_d=1/f,grid_points=grid,semiamplitude_mag=float(np.hypot(*v[1:])))
    ft=1/clocks["early"]["period_d"]; vt=fit(t[train],y[train],w[train],ft)[0]
    clocks["early_predict_late"]=prediction(t[~train],y[~train],w[~train],ft,vt)
    f=1/clocks["all"]["period_d"];coeff=fit(t,y,w,f)[0]
    # Resample entire observing months within the chosen period peak.
    block=np.floor((t-EPOCH)/30).astype(int)
    groups=[np.flatnonzero(block==k) for k in np.unique(block)]
    periods=[]
    for _ in range(600):
        ind=np.concatenate([groups[k] for k in RNG.integers(len(groups),size=len(groups))])
        periods.append(1/refine(t[ind],y[ind],w[ind],f,.4/np.ptp(t)))
    clocks["block_bootstrap_95_d"]=np.quantile(periods,[.025,.975]).tolist()
    clocks["max_brightness_bjd_near_epoch"]=float(EPOCH+circular(phase(coeff)+.5)/f)
    # Nearby sources and individual exposure metadata are controls for a blended variable or systematic.
    controls={}
    for oid,group in ztf.groupby("oid"):
        if len(group)<40 or oid==333114400007906:continue
        g=clean(group);tc=bjd(g.hjd.to_numpy(),"sun")
        controls[str(oid)]=dict(separation_arcsec=float(SKY.separation(SkyCoord(g.ra.median(),g.dec.median(),unit=u.deg)).arcsec),
                                **prediction(tc,g.res.to_numpy(),1/(g.magerr.to_numpy()**2+.06**2),f,coeff))
    gaia=pd.read_csv(DATA/"gaia-field.csv",dtype={"source_id":str})
    pos=SkyCoord(target.ra.median(),target.dec.median(),unit=u.deg)
    sep=pos.separation(SkyCoord(gaia.ra,gaia.dec,unit=u.deg)).arcsec
    ordered=np.argsort(sep)
    nearest=[dict(id=gaia.source_id.iloc[ix],offset_arcsec=float(sep[ix])) for ix in ordered[:2]]
    sectors={s:sector_stats(pd.read_csv(DATA/f"tess-{s}.csv.gz"),f,coeff) for s in (80,92)}
    ps=pd.read_csv(DATA/"panstarrs.csv");ps=ps[(ps.psfFlux>0)&(ps.psfFluxErr<.17*ps.psfFlux)].copy()
    ps["mag"]=-2.5*np.log10(ps.psfFlux/3631)
    py=(ps.mag-ps.groupby("filterID").mag.transform("median")).to_numpy()
    ps1=prediction(bjd(ps.obsTime.to_numpy()+2400000.5,"earth"),py,np.ones(len(py)),f,coeff)
    report={"ztf":clocks,"tess":sectors,"source_position":dict(ra=float(pos.ra.deg),dec=float(pos.dec.deg),gaia=nearest),
            "phase_dependent_photocenter":photocenter(rows,t,f,coeff), "field_controls":controls,"panstarrs_pre_eruption":ps1}
    (ROOT/"measurements.json").write_text(json.dumps(report,indent=2)+"\n")
    return report, rows, t, y, f, coeff


if __name__=="__main__":
    main()
