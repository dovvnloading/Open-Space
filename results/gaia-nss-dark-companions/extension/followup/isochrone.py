"""Fit a giant's Gaia/2MASS SED and independent temperatures to rotating MIST v1.2 tracks."""
import argparse
import csv
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
URL='https://mist.science/data/tarballs_v1.2/MIST_v1.2_vvcrit0.4_UBVRIplus.txz'
# Model columns from MIST v1.2 photometric grid, including Gaia EDR3 passbands.
COLS=[0,1,2,3,4,5,7,30,31,32,14,15,16]
BANDS=['G','BP','RP','J','H','Ks']
EXT=np.array([.82,1.02,.61,.282,.175,.112])
FM=.640  # solar masses, joint SB1+LAMOST fit
METALS=('m0.75','m0.50','m0.25','p0.00','p0.25')


def m2(f,m1):
    lo=np.zeros_like(m1)+f; hi=f+2*m1
    for _ in range(40):
        mid=(lo+hi)/2
        below=mid**3 < f*(m1+mid)**2
        lo,hi=np.where(below,mid,lo),np.where(below,hi,mid)
    return (lo+hi)/2


def grid(cache, compressed):
    if compressed.exists():
        return np.load(compressed)['grid']
    if not cache.exists():
        urllib.request.urlretrieve(URL,cache)
    grids=[]
    with tarfile.open(cache,'r|xz') as arc:
        for member in arc:
            if not member.isfile() or not any('_feh_'+m+'_' in member.name for m in METALS): continue
            with arc.extractfile(member) as f:
                data=np.loadtxt(f,usecols=COLS)
            grids.append(data)
    data=np.concatenate(grids)
    np.savez_compressed(compressed,grid=data)
    return data


def systems(data):
    values=[]
    for met in np.unique(data[:,6]):
        rows=data[data[:,6]==met]
        for age in np.unique(rows[:,1]):
            iso=rows[rows[:,1]==age]
            m=iso[:,2]; weights=np.abs(np.gradient(m))*m**-2.35*10**(age-9)
            good=(iso[:,4]>np.log10(4350))&(iso[:,4]<np.log10(5550))&(iso[:,5]>1.0)&(iso[:,5]<3.4)&(weights>0)&(iso[:,0]<808)
            prim=iso[good]
            if not len(prim): continue
            # The ordinary star alternative is a coeval hydrogen-burning companion.
            ms=iso[(iso[:,0]>=202)&(iso[:,0]<=454)]
            ms=ms[np.argsort(ms[:,3])]
            sec=m2(FM,prim[:,3])
            for n in (0,1,2):
                mags=prim[:,7:].copy()
                if n:
                    if len(ms)<2: continue
                    smag=np.stack([np.interp(sec/n,ms[:,3],ms[:,j],left=np.nan,right=np.nan) for j in range(7,13)],axis=1)
                    mags=-2.5*np.log10(10**(-.4*mags)+n*10**(-.4*smag))
                values.append(np.column_stack([np.full(len(prim),n),prim[:,:7],sec,weights[good],mags]))
    return np.concatenate(values)


def main(cache,compressed):
    s=systems(grid(cache,compressed))
    star=next(csv.DictReader((HERE/'source.csv').open()))
    tm=next(csv.DictReader((HERE/'2mass.csv').open()))
    phot=np.array([float(star[k]) for k in ['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']]
                  +[float(tm[k]) for k in ['j_m','h_m','k_m']])
    plx=float(star['parallax'])
    dm=10-5*np.log10(plx)
    sdm=(5/np.log(10))*float(star['parallax_error'])/plx
    # Systematic photometric floors include passband/extinction approximations.
    errors=np.array([.055,.07,.06,.04,.04,.04])
    lasp=list(csv.DictReader((HERE/'lasp_visits.csv').open(encoding='utf-8-sig')))
    spec=next(csv.DictReader((HERE/'star.csv').open()))
    params={'GSP-Spec':tuple(float(spec[k]) for k in ['teff_gspspec','logg_gspspec','mh_gspspec'])}
    for label,keys in [('LASP',['teff_lasp','logg_lasp','feh_lasp']),('LAMOST CNN',['teff_cnn','logg_cnn','feh_cnn'])]:
        params[label]=tuple(np.median([float(r[k]) for r in lasp]) for k in keys)
    results=[]
    for ir_only in (False,True):
        j=slice(3,None) if ir_only else slice(None)
        A=np.column_stack((EXT[j],np.ones(len(EXT[j]))))
        w=errors[j]**-2
        precision=A.T@(w[:,None]*A)+np.diag([0.,sdm**-2])
        resid=(phot-dm-s[:,10:])[:,j]
        par=((resid*w)@A)@np.linalg.inv(precision)
        av=np.clip(par[:,0],0,6)
        shift=((resid-av[:,None]*EXT[j])*w).sum(1)/(w.sum()+sdm**-2)
        cost=np.sum(((resid-av[:,None]*EXT[j]-shift[:,None])/errors[j])**2,axis=1)+(shift/sdm)**2
        for label,(t,g,z) in params.items():
            # Conservative single-star atmosphere systematics. Evaluate the disagreeing pipelines separately.
            chi=cost+((10**s[:,5]-t)/150)**2+((s[:,6]-g)/.25)**2+((s[:,7]-z)/.25)**2
            for n in (0,1,2):
                valid=(s[:,0]==n)&np.isfinite(chi)
                idx=np.flatnonzero(valid); i=idx[np.argmin(chi[idx])]
                weight=np.exp(-.5*(chi[idx]-chi[i]))*s[idx,9]
                weight/=sum(weight)
                def q(k):
                    v=s[idx,k]; order=np.argsort(v)
                    return np.interp([.16,.5,.84],np.cumsum(weight[order]),v[order]).round(4).tolist()
                results.append(dict(spectra=label,primary_constraints=[t,g,z],bands=BANDS[j],number_of_normal_secondaries=n,
                    min_chi2=round(float(chi[i]),3),best_av=round(float(av[i]),3),best_distance_modulus=round(float(dm+shift[i]),3),
                    best_current_m1=round(float(s[i,4]),3),best_edge_on_m2=round(float(s[i,8]),3),
                    m1_q16_50_84=q(4),m2_edge_on_q16_50_84=q(8),conditional_mass_weight_above3=round(float(weight[s[idx,8]>3].sum()),4),
                    best_band_residual_mag=(resid[i]-av[i]*EXT[j]-shift[i]).round(4).tolist()))
    with cache.open('rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
    (HERE/'isochrone_fits.json').write_text(json.dumps(dict(url=URL,sha256=sha,grid_version='MIST v1.2, v/vcrit=0.4',
        mass_function_fixed_msun=FM,models=results),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('grid_cache',type=Path)
    a=p.parse_args();main(a.archive,a.grid_cache)
