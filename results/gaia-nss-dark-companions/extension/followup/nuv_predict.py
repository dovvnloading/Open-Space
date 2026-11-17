"""Held-out GALEX NUV posterior prediction from the optical/IR MIST fits."""
import argparse
import csv
import hashlib
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from isochrone import HERE, EXT, systems

BC_URL='https://mist.science/BC_tables/v1/GALEX.txz'


def tables(path):
    if not path.exists():urllib.request.urlretrieve(BC_URL,path)
    grids={}
    with tarfile.open(path,'r:xz') as arc:
        for member in arc.getmembers():
            if not member.isfile():continue
            with arc.extractfile(member) as f: a=np.loadtxt(f)
            z=a[0,2]
            if z not in [-.75,-.50,-.25,0.,.25]:continue
            t,g,av=map(np.unique,[a[:,0],a[:,1],a[:,3]])
            shape=(len(t),len(g),len(av)); mat=np.full(shape,np.nan)
            mat[np.searchsorted(t,a[:,0]),np.searchsorted(g,a[:,1]),np.searchsorted(av,a[:,3])]=a[:,6]
            assert np.isfinite(mat).all() and np.all(a[:,4]==3.1)
            grids[z]=RegularGridInterpolator((np.log10(t),g,av),mat,bounds_error=True)
    return grids


def fit_data(s):
    star=next(csv.DictReader((HERE/'source.csv').open()))
    ir=next(csv.DictReader((HERE/'2mass.csv').open()))
    obs=np.array([float(star[k]) for k in ('phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag')]+
                 [float(ir[k]) for k in ('j_m','h_m','k_m')])
    plx=float(star['parallax']); dm=10-5*np.log10(plx)
    sdm=(5/np.log(10))*float(star['parallax_error'])/plx
    error=np.array([.055,.07,.06,.04,.04,.04]);w=error**-2
    A=np.column_stack([EXT,np.ones(len(EXT))])
    cov=np.linalg.inv(A.T@(w[:,None]*A)+np.diag([0,sdm**-2]))
    residual=obs-dm-s[:,10:]
    par=((residual*w)@A)@cov
    av=np.clip(par[:,0],0,6)
    ddm=((residual-av[:,None]*EXT)*w).sum(1)/(w.sum()+sdm**-2)
    chi=np.sum(((residual-av[:,None]*EXT-ddm[:,None])/error)**2,axis=1)+(ddm/sdm)**2
    return chi,av,dm+ddm,cov


def nuv_bc(grids,t,g,z,av):
    result=np.zeros(len(t))
    for metal in np.unique(z):
        k=z==metal
        result[k]=grids[metal](np.column_stack((t[k],g[k],av[k])))
    return result


def mag(grids,m,t,g,z,av):
    # Solar g in cgs, nominal bolometric magnitude and nominal effective temperature.
    log_lum=np.log10(m)+4.438-g+4*(t-np.log10(5772))
    return 4.74-2.5*log_lum-nuv_bc(grids,t,g,z,av)


def predict(s,draw,grids,tracks,n,av,dm):
    rows=s[draw]
    m,t,g,z=rows[:,4],rows[:,5],rows[:,6],rows[:,7]
    flux=10**(-.4*mag(grids,m,t,g,z,av))
    if n:
        sm=rows[:,8]/n
        st=np.empty(len(sm)); sg=np.empty(len(sm))
        for metal,age in set(zip(z,rows[:,2])):
            k=(z==metal)&(rows[:,2]==age)
            zz=tracks[metal,age]
            st[k]=np.interp(sm[k],zz[:,3],zz[:,4])
            sg[k]=np.interp(sm[k],zz[:,3],zz[:,5])
        flux+=n*10**(-.4*mag(grids,sm,st,sg,z,av))
    return dm-2.5*np.log10(flux)


def run(star_path,bc_path):
    data=np.load(star_path)['grid']
    s=systems(data); grids=tables(bc_path)
    tracks={}
    for z in np.unique(data[:,6]):
        met=data[(data[:,6]==z)&(data[:,0]>=202)&(data[:,0]<=454)]
        for age in np.unique(met[:,1]):
            rows=met[met[:,1]==age]
            tracks[z,age]=rows[np.argsort(rows[:,3])]
    chi0,av,dm,cov=fit_data(s)
    opts=json.loads((HERE/'isochrone_fits.json').read_text())['models']
    fitpars={r['spectra']:r['primary_constraints'] for r in opts if len(r['bands'])==6}
    obs=json.loads((HERE/'galex.json').read_text())['data'][0]
    rng=np.random.default_rng(456091)
    out=[]
    for label,(t,g,z) in fitpars.items():
        chi=chi0+((10**s[:,5]-t)/150)**2+((s[:,6]-g)/.25)**2+((s[:,7]-z)/.25)**2
        for n in [0,1,2]:
            k=np.flatnonzero((s[:,0]==n)&np.isfinite(chi))
            w=np.exp(-.5*(chi[k]-np.min(chi[k])))*s[k,9];w/=w.sum()
            draw=rng.choice(k,20000,p=w)
            shift=rng.multivariate_normal([0,0],cov,len(draw))
            v=av[draw]+shift[:,0];d=dm[draw]+shift[:,1]
            assert np.all((v>0)&(v<6))
            y=predict(s,draw,grids,tracks,n,v,d)
            intrinsic=predict(s,draw,grids,tracks,n,v*0,d)
            attenuation=y-intrinsic
            def q(x):return np.quantile(x,[.025,.16,.5,.84,.975]).round(4).tolist()
            out.append(dict(spectra=label,n_normal_dwarfs=n,nuv_predicted=q(y),
                    A_nuv=q(attenuation),required_delta_A_nuv=q(obs['nuv_mag']-y),
                    prediction_if_A_nuv_minus20pct=q(y-.2*attenuation),
                    prediction_if_A_nuv_plus20pct=q(y+.2*attenuation)))
    digests={}
    for p in (star_path,bc_path,HERE/'galex.json',HERE/'isochrone_fits.json'):
        with p.open('rb') as f:digests[p.name]=hashlib.file_digest(f,'sha256').hexdigest()
    return dict(uv_source=obs['IAUName'],nuv_ab=obs['nuv_mag'],err=obs['nuv_magerr'],nuv_exptime_s=obs['nuv_exptime'],
                fuv_exptime_s=obs['fuv_exptime'],bolometric_correction_url=BC_URL,sha256=digests,
                draws_per_model=20000,quantiles=[.025,.16,.5,.84,.975],models=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mist_rows',type=Path);p.add_argument('bc_archive',type=Path)
    a=p.parse_args()
    (HERE/'nuv_predictions.json').write_text(json.dumps(run(a.mist_rows,a.bc_archive),indent=2)+'\n')
