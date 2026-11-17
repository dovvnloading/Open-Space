"""Measure narrow cube residuals in separately dated SPHEREx detector cutouts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import io
import json
import logging
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import map_coordinates
from index_visits import approximate_xy, cut_url, read_url, sha

logging.getLogger('astropy').setLevel(logging.ERROR)
SOURCE = 1 << 21
OFFSETS = [(9*np.cos(a),9*np.sin(a)) for a in np.arange(8)*np.pi/4]


def choose_images(index, calibration, ra, dec, center, scale, per_bin, all_options=False):
    wl = RegularGridInterpolator((calibration['y'],calibration['x']),calibration['wave_um'])
    choices = {}
    with index.open() as f:
        for row in csv.DictReader(f):
            uri = row['uri']
            if '/dr1/' not in uri and '/qr3/' not in uri:
                continue
            release = 'qr3' if '/qr3/' in uri else 'dr1'
            x,y = approximate_xy(row['poly'],ra,dec)
            if not (22 < x < 2016 and 22 < y < 2016):
                continue
            prediction=float(wl([[y+1,x+1]])[0])
            distance=(prediction-center)/scale
            for label,offset in (('blue',-2.25),('core',0.),('red',2.25)):
                tolerance=.75 if label!='core' else .55
                if abs(distance-offset)<=tolerance:
                    choices.setdefault((release,label),[]).append((float(row['time_bounds_lower']),abs(distance-offset),uri,prediction))
    records=[]
    for (release,label),entries in sorted(choices.items()):
        # First maximize independent dates, then proximity to the requested wavelength.
        if all_options:
            for t,delta,uri,pred in sorted(entries,key=lambda a:(a[1],a[0],a[2])):
                records.append(dict(release=release,bin=label,mjd=t,uri=uri,pred_um=round(pred,6)))
            continue
        days={}
        for t,delta,uri,pred in sorted(entries,key=lambda a:(a[1],a[0],a[2])):
            days.setdefault(int(t),(t,delta,uri,pred))
        candidates=sorted(days.values())
        n=min(len(candidates),per_bin)
        for j in np.round(np.linspace(0,len(candidates)-1,n)).astype(int):
            t,delta,uri,pred=candidates[j]
            records.append(dict(release=release,bin=label,mjd=t,uri=uri,pred_um=round(pred,6)))
    return records


def aperture(img,zodi,var,flags,weights,sky):
    good=np.isfinite(img)&np.isfinite(zodi)&np.isfinite(var)&(var>0)&((flags&~SOURCE)==0)
    ap=weights>0
    sky=sky&good
    good_ap=np.sum(weights[good])
    total=np.sum(weights)
    flag_or=int(np.bitwise_or.reduce(flags[ap]))
    if good_ap < .9999*total or sky.sum()<28:
        return None,flag_or,total,good_ap
    bkg=float(np.median((img-zodi)[sky]))
    signal=np.sum(weights[ap]*(img[ap]-zodi[ap]-bkg))
    noise=np.sqrt(np.sum(var[ap]*weights[ap]**2)+total**2*np.pi/2*np.median(var[sky])/sky.sum())
    zod=np.sum(weights[ap]*(zodi[ap]-np.median(zodi[sky])))
    return (signal,noise,zod,bkg),flag_or,total,good_ap


def measure(hd,ra,dec,maps):
    h=hd['IMAGE'].header
    x,y=WCS(h).all_world2pix([[ra,dec]],0)[0]
    sp=WCS(h,hd,key='W')
    sp.sip=None
    quick_wave,quick_band=sp.all_pix2world([[x,y]],0)[0]
    yy,xx=np.indices(hd['IMAGE'].data.shape)
    arrays=[hd[name].data for name in ('IMAGE','ZODI','VARIANCE','FLAGS')]
    xx,yy=xx-x,yy-y
    sub=(np.arange(5)+.5)/5-.5
    estimates=[]
    detx,dety=x+1-h['CRPIX1W'],y+1-h['CRPIX2W']
    wavelength,bandwidth=[float(map_coordinates(a,[[dety],[detx]],order=1)[0]) for a in maps]
    result=dict(quick_wave_um=round(float(quick_wave),6),quick_bandwidth_um=round(float(quick_band),6),date_obs=h['DATE-OBS'],obs_id=h['OBSID'],detector=int(h['DETECTOR']),
                det_x=round(detx,2),det_y=round(dety,2),wave_um=round(float(wavelength),6),
                bandwidth_um=round(float(bandwidth),6),flag_or='',coverage='',n_blank=0,
                flux='',sigma='',blank_median='',zodi_aperture='',sky_median='',status='aperture_flags')
    if min(x,y,xx.shape[1]-1-x,xx.shape[0]-1-y)<17:
        result['status']='edge'
        return result
    source=None
    for i,(dx,dy) in enumerate([(0,0),*OFFSETS]):
        r=np.hypot(xx-dx,yy-dy)
        weights=sum((xx-dx+ox)**2+(yy-dy+oy)**2<=2.25**2 for ox in sub for oy in sub)/25.
        value,flags,total,good=aperture(*arrays,weights,(r>=4)&(r<=7))
        if i==0:
            result.update(flag_or=flags,coverage=round(good/total,4))
            source=value
        elif value is not None:
            estimates.append(value[0])
    result['n_blank']=len(estimates)
    if source is not None:
        flux,formal,zod,bg=source
        result.update(flux=round(float(flux),6),zodi_aperture=round(float(zod),6),sky_median=round(float(bg),6))
        if len(estimates)>=3:
            med=np.median(estimates)
            err=max(formal,1.4826*np.median(np.abs(estimates-med)))
            result.update(blank_median=round(float(med),6),sigma=round(float(err),6),status='ok')
        else:
            result['status']='few_blanks'
    return result


def fetch_one(record,ra,dec,center,scale,cache,maps):
    rec=dict(record)
    url=cut_url(rec.pop('uri'),ra,dec)
    rec['url']=url
    local=cache/(sha(url.encode())+'.fits')
    try:
        if local.exists():
            data=local.read_bytes()
        else:
            data=read_url(url)
            if not data.startswith(b'SIMPLE  '):
                raise ValueError(data[:100])
            local.write_bytes(data)
        rec['sha256']=sha(data)
        with fits.open(io.BytesIO(data)) as hd:
            rec.update(measure(hd,ra,dec,maps))
        distance=(rec['wave_um']-center)/scale
        if abs(distance-{'blue':-2.25,'core':0.,'red':2.25}[rec['bin']])>(.55 if rec['bin']=='core' else .75):
            rec['status']='off_bin'
    except Exception as exc:
        rec['status']=f'read_error: {exc}'
    return rec


def analyse(rows):
    result={}
    for release in ('qr3','dr1'):
        selected=[r for r in rows if r['release']==release and r['status']=='ok']
        for subgroup in ('all','early','late'):
            # Split by observation time across all bins within the release, before looking at flux.
            mjd=[r['mjd'] for r in rows if r['release']==release]
            split=np.median(mjd) if mjd else np.nan
            part=[r for r in selected if subgroup=='all' or ((r['mjd']<split)==(subgroup=='early'))]
            bands={}
            for label in ('blue','core','red'):
                rr=[r for r in part if r['bin']==label]
                if not rr:
                    continue
                vals=np.array([r['flux'] for r in rr])
                # Single-bin median with empirical inter-date scatter, blank/formal noise,
                # and a 3% floor on the flux to allow for uncorrected aperture variation.
                sigma=max(1.4826*np.median(abs(vals-np.median(vals))),np.median([r['sigma'] for r in rr]))*1.253/np.sqrt(len(rr))
                bands[label]=dict(n=len(rr),dates=len({r['date_obs'][:10] for r in rr}),wave=float(np.median([r['wave_um'] for r in rr])),
                                  median=float(np.median(vals)),err=float(np.hypot(sigma,.03*abs(np.median(vals)))))
            summary=dict(bands=bands)
            if len(bands)==3:
                l,c,r=[bands[k] for k in ('blue','core','red')]
                t=(c['wave']-l['wave'])/(r['wave']-l['wave'])
                continuum=(1-t)*l['median']+t*r['median']
                diff=c['median']-continuum
                err=np.sqrt(c['err']**2+((1-t)*l['err'])**2+(t*r['err'])**2)
                summary.update(fractional_contrast=diff/continuum if continuum>0 else None,
                               contrast=float(diff),err=float(err),signed_score=float(diff/err))
            result[release+'_'+subgroup]=summary
    return result


def run(args):
    top=list(csv.DictReader(args.targets.open()))
    sources={r['source_id']:r for r in csv.DictReader(args.sources.open())}
    args.out.mkdir(parents=True,exist_ok=True)
    cache=args.index_cache/'cutouts'
    cache.mkdir(exist_ok=True)
    for target in top:
        name=target['source_id']
        src=sources[name]
        ra=float(src['extraction_ra']);dec=float(src['extraction_dec'])
        center=float(target['wave_um']);scale=float(target['scale_um']);d=int(target['detector'])
        suffix=target.get('label','')
        stem=name+('_'+suffix if suffix else '')
        if args.replay_dir:
            with (args.replay_dir/f'{stem}.csv').open() as f:
                receipts=list(csv.DictReader(f))
            records=[dict(release=r['release'],bin=r['bin'],mjd=float(r['mjd']),
                          uri=r['url'].split('ipac.caltech.edu/',1)[1].split('?',1)[0],
                          pred_um=float(r['pred_um'])) for r in receipts]
        else:
            cal=json.loads((args.calibration/f'calibration_D{d}.json').read_text())
            records=choose_images(args.index_cache/f'footprints-D{d}.csv',cal,ra,dec,center,scale,args.per_bin)
        print(name,len(records),flush=True)
        with fits.open(args.science_cache/f'sci_D{d}.fits') as science, ThreadPoolExecutor(args.workers) as pool:
            maps=[science[k].data for k in ('CWAVE','CBAND')]
            rows=list(pool.map(lambda r:fetch_one(r,ra,dec,center,scale,cache,maps),records))
        if args.replay_dir:
            for row,original in zip(rows,receipts):
                expected=original.get('sha256')
                if expected and row.get('sha256')!=expected:
                    row['status']='archive_bytes_changed_or_unavailable'
        keys=list(dict.fromkeys(k for row in rows for k in row)) or ['release','bin','mjd','pred_um','url','sha256','status']
        with (args.out/f'{stem}.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        summary=dict(source_id=name,ra=ra,dec=dec,center_um=center,scale_um=scale,
                     attempted=len(rows),status_counts={s:sum(r['status']==s for r in rows) for s in sorted({r['status'] for r in rows})},groups=analyse(rows))
        (args.out/f'{stem}.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        print(name,summary['status_counts'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('targets','sources','index_cache','calibration','out','science_cache'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--replay-dir',type=Path,help='read the archived selection from this directory and check hashes')
    parser.add_argument('--per-bin',type=int,default=6)
    parser.add_argument('--workers',type=int,default=6)
    run(parser.parse_args())
