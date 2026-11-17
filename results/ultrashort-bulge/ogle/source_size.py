"""Conditional angular-source-size sensitivity for two rho-only KMT flags.

Use published VVV J/Ks catalogs and the released 2024/2025 pySIS light curves.
No single-lens angular Einstein radius or calibrated event color is inferred.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tarfile
import urllib.parse
import urllib.request

import numpy as np
from scipy.optimize import curve_fit

from counterparts import HERE,coords
sys.path.insert(0,str(HERE.parent))
import fit
sys.path.insert(0,str(HERE.parent/"kmt"))
from neighbours import query as query_virac

TARGETS=('KMT-2024-BLG-1516','KMT-2025-BLG-2688')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def cadence_colour(path,geom):
    vb=fit.engine()
    out=[]
    with tarfile.open(path) as tar:
        for m in tar:
            if not m.name.endswith('.pysis') or not m.isfile():continue
            t,df,err=np.loadtxt(tar.extractfile(m),usecols=(0,1,2),ndmin=2).T
            t=t-geom['t0']
            ok=(err>0)&np.isfinite(df)&np.isfinite(err)&np.isfinite(t)
            if not ok.any():continue
            t,y,e=t[ok],-df[ok],err[ok]
            u,te,rho=geom['u0_or_fall_d'],geom['te_or_rise_d'],geom['rho']
            a=np.asarray(vb.ESPLLightCurve([np.log(u),np.log(te),0.,np.log(rho)],t.tolist())[0])-1
            near=abs(t)<max(1.,rho)*te
            off=abs(t)>max(3*rho*te,3*te)
            if sum(off)<5:continue
            scale=max(1.,np.median(np.abs((y[off]-np.median(y[off]))/e[off]))/.67449)
            w=1/(scale*e)**2
            centered=a-np.average(a,weights=w)
            v=(w*centered**2).sum()
            amp=float((w*centered*y).sum()/v)
            out.append(dict(stream=m.name,cat_band=m.name.split('_')[-1][0],
                            n=int(len(t)),n_in_source_crossing=int(sum(near)),
                            source_flux=amp,formal_error=float(1/np.sqrt(v)),noise_scale=float(scale)))
    pairs=[]
    for i in out:
        if i['cat_band']!='I':continue
        v=next((b for b in out if b['stream']==i['stream'].replace('_I.','_V.')),None)
        if v:
            pairs.append(dict(stream=i['stream'].split('_')[0],nV_inside_crossing=v['n_in_source_crossing'],
                              source_V_I_count_ratio=round(v['source_flux']/i['source_flux'],4),
                              V_amplitude_significance=round(v['source_flux']/v['formal_error'],2)))
    return pairs


def clump(phot):
    radius,j,k=phot.T
    records=[]
    for rad in [60.,90.,120.]:
        local=(radius>2)&(radius<=rad)&(j-k>.8)&(j-k<1.8)
        edges=np.arange(12.,14.71,.125); bins=(edges[1:]+edges[:-1])/2
        hist=np.histogram(k[local],edges)[0]
        def density(m,a,b,h,mean,sd):
            return a*np.exp(b*(m-13.3))+h*np.exp(-.5*((m-mean)/sd)**2)
        p,_=curve_fit(density,bins,hist,p0=[max(1,np.median(hist)),.3,max(3,np.ptp(hist)),13.25,.22],
                      bounds=([0,-1,0,12.65,.12],[np.inf,1,np.inf,13.8,.5]),
                      sigma=np.sqrt(np.maximum(1,hist)),maxfev=3000)
        sel=local&(np.abs(k-p[3])<.22)
        c=j[sel]-k[sel]
        records.append(dict(radius_arcsec=rad,n=int(len(c)),K_peak=round(float(p[3]),3),
                            J_K_centroid=round(float(np.median(c)),3),
                            color_p16_p84=np.percentile(c,[16,84]).round(3).tolist()))
    return records


def sensitivity(k,clump_centroids):
    # A sensitivity envelope, not a sampled posterior. Allow a differential
    # color-excess offset of ±0.05 mag, baseline K bandpass offset ±0.05 mag,
    # and the entire valid V-K range of the giant-star calibration.
    res=[]
    for source_fraction in [1.,.5,.1]:
        vals=[]
        for reddening in (min(clump_centroids)-.68-.05,max(clump_centroids)-.68+.05):
            for coeff in [.428,.528]:
                for dk in [-.05,.05]:
                    k0=k-2.5*np.log10(source_fraction)-coeff*reddening+dk
                    for c in [2.049,3.839]:
                        for scatter in [-.019,.019]:
                            vals.append(500*10**(.562+.051*c-.2*k0+scatter))
        res.append(dict(hypothetical_K_source_fraction=source_fraction,
                        radius_uas_min=round(min(vals),2),radius_uas_max=round(max(vals),2)))
    return res


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--archives',type=Path,required=True)
    args=parser.parse_args();cache=args.cache;cache.mkdir(exist_ok=True,parents=True)
    out=[]
    for event in TARGETS:
        yr=event.split('-')[1]; path=HERE.parent/'kmt'/f'catalog-{yr}.tsv'
        r=next(r for r in csv.DictReader(path.open(),delimiter='\t') if r['Event']==event)
        ra,dec=coords(r['RA'],r['Dec'])
        pars={'-source':'II/348/vvv2','-c':f'{ra:.9f} {dec:+.9f}','-c.rs':120,'-out.add':'_r','-sort':'_r',
              '-out.max':30000,'-out':'srcid,raJ2000,deJ2000,Jmag3,Hmag3,Ksmag3,Jperrbits,Ksperrbits'}
        url='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urllib.parse.urlencode(pars)
        field=cache/(event+'.tsv')
        if not field.exists():field.write_bytes(urllib.request.urlopen(url,timeout=100).read())
        raw=field.read_bytes();lines=raw.decode().splitlines()
        idx=next(i for i,l in enumerate(lines) if l.startswith('_r\t'))
        cols=lines[idx].split('\t'); rows=[dict(zip(cols,l.split('\t'))) for l in lines[idx+3:] if l.strip() and not l.startswith('#')]
        valid=[r for r in rows if r['Jmag3'].strip() and r['Ksmag3'].strip() and int(r['Jperrbits'])<256 and int(r['Ksperrbits'])<256]
        phot=np.array([[float(r['_r']),float(r['Jmag3']),float(r['Ksmag3'])] for r in valid])
        c=clump(phot)
        virac_cache=cache/'virac2'
        virac_cache.mkdir(exist_ok=True)
        query_virac(dict(event=event,ra=r['RA'],dec=r['Dec']),'VIRAC2',virac_cache)
        virac=json.loads((virac_cache/(event+'-VIRAC2.json')).read_text())
        if virac.get('status')!=200:raise RuntimeError(virac)
        colors=[x['J_K_centroid'] for x in c]
        for obj in virac['matches']:
            if obj['Ksmag'].strip():
                obj['conditional_giant_radii']=sensitivity(float(obj['Ksmag']),colors)
        archive=args.archives/(event+'.tar.gz')
        geom=json.loads((HERE.parent/'peaks'/(event+'.json')).read_text())['models']['fs']
        out.append(dict(event=event,url=url,field_sha256=digest(raw),total_rows=len(rows),valid_JK=len(phot),
                        clump=c,virac2=virac,nearest_2arcsec=[r for r in rows if float(r['_r'])<=2],
                        geometry_condition='FSPL at the shared peak-fit parameters',
                        pySIS_sha256=digest(archive.read_bytes()),
                        differential_color=cadence_colour(archive,geom),source_radius_cases=sensitivity(phot[0,2],[x['J_K_centroid'] for x in c])))
    p=cache/'bundle.tgz'
    u='https://www.astrouw.edu.pl/ogle/ogle4/ews/2024/blg-0859.tar.gz'
    if not p.exists():p.write_bytes(urllib.request.urlopen(u,timeout=90).read())
    with tarfile.open(p) as tf:raw=tf.extractfile(next(m for m in tf if m.name.endswith('map.dat'))).read()
    table=np.loadtxt(io.BytesIO(raw));star=table[table[:,0]==214764][0]
    out[0]['OGLE_map']=dict(url=u,bundle_sha256=digest(p.read_bytes()),map_sha256=digest(raw),field_stars=len(table),
                            star_ID=214764,baseline_V=star[3],baseline_V_I=star[4],baseline_I=star[5])
    (HERE/'source-size-cases.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
