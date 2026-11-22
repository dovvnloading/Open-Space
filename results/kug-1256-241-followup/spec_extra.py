"""Targeted optical follow-up and offline replay of deposited science pixels."""
import csv
import io
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urljoin

import astropy.units as u
import numpy as np
from astropy.io.votable import parse
from bs4 import BeautifulSoup
from scipy.linalg import lstsq

import spectra as sp
from common import OUT, HOST_RA, HOST_DEC, fetch, status, event

PIX = OUT / 'optical_pixels'
ZREF = .07468937  # Pass-one BOSS redshift, spec-5991-56076-0484.fits.
META = []
INDEX = []


def store(survey, ident, wave, flux, ivar, z, mjd, end, unit, sep, note):
    PIX.mkdir(exist_ok=True)
    label=sp.keyname(ident)
    w=np.asarray(wave,dtype=float).ravel()
    f=np.asarray(flux,dtype=float).ravel()
    iv=np.asarray(ivar,dtype=float).ravel() if ivar is not None else np.full(w.shape,np.nan)
    if w.shape!=f.shape or w.shape!=iv.shape or len(w)<50:
        raise ValueError('unusable pixel arrays')
    np.savez_compressed(PIX/(label+'.npz'),wave=w,flux=f,ivar=iv)
    m=dict(archive=survey,product=ident,file=label+'.npz',z=float(z),mjd=float(mjd) if mjd!='' else '',
           end_mjd=float(end) if end!='' else '',unit=unit,separation_arcsec=round(sep,3),note=note)
    META.append(m)
    event(survey,'optical spectrum (pixel record)',mjd,kind='spectrum',source=ident,end_mjd=end,
          note=f'{len(w)} wavelength pixels deposited; {note}')


def boss():
    for p,m,fi,z in [(5991,56076,484,.07468937),(5989,56312,965,.07488265)]:
        ident=f'spec-{p:04d}-{m}-{fi:04d}.fits'
        url=f'https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/{p:04d}/{ident}'
        with sp.unpack(fetch('pixel-copy-'+ident,url)) as h:
            t=h['COADD'].data
            iv=t['ivar'].copy()
            iv[t['and_mask']!=0]=0
            store('BOSS DR17',ident,10**t['loglam'],t['flux'],iv,z,m,'',sp.STD,.449,
                  'Optical nuclear spectrum; MAST serves duplicate reductions of these dates.')


def dirs(key,url):
    raw=fetch(key,url,timeout=18)
    soup=BeautifulSoup(raw,'html.parser')
    return sorted({a['href'].strip('/') for a in soup.select('a[href]')
                   if re.fullmatch(r'[A-Za-z0-9_-]+/',a['href'])})


def desi_direct():
    from astropy_healpix import HEALPix
    root='https://data.desi.lbl.gov/public/dr1/spectro/redux/iron'
    hp=HEALPix(nside=64,order='nested')
    pixel=int(hp.lonlat_to_healpix(HOST_RA*u.deg,HOST_DEC*u.deg))
    edge=hp.cone_search_lonlat(HOST_RA*u.deg,HOST_DEC*u.deg,10*u.arcsec)
    pixels=sorted(set(map(int,edge))|{pixel})
    base=root+'/healpix/'
    try:
        surveys=dirs('desi-directory-surveys',base) or ['main','sv1','sv2','sv3','special','cmx']
    except Exception:
        surveys=['main','sv1','sv2','sv3','special','cmx']
    jobs=[]
    for survey in surveys:
        try:
            programs=dirs('desi-directory-'+survey,base+survey+'/') or ['bright','dark','backup','other']
        except Exception:
            programs=['bright','dark','backup','other']
        for program in programs:
            for pix in pixels:
                jobs.append((survey,program,pix))
    svc='DESI DR1 direct HEALPix'
    def process(job):
        survey,program,pix=job
        prefix=f'{survey}-{program}-{pix}'
        directory=f'{base}{survey}/{program}/{pix//100}/{pix}/'
        try:
            body=fetch('desi-redrock-'+prefix,directory+'redrock-'+prefix+'.fits',timeout=35)
        except Exception as e:
            if getattr(getattr(e,'response',None),'status_code',None)==404:
                return prefix,'no released redrock file',0
            return prefix,'unavailable '+str(e)[:180],0
        try:
            with sp.unpack(body) as h:
                fm=h['FIBERMAP'].data
                red=h['REDSHIFTS'].data
                match=[(int(r['TARGETID']),sp.dist(r['TARGET_RA'],r['TARGET_DEC'])) for r in fm
                       if sp.dist(r['TARGET_RA'],r['TARGET_DEC'])<=3.]
                ztab={int(r['TARGETID']):float(r['Z']) for r in red if r['ZWARN']==0 and 0<r['Z']<1}
            if not match:
                return prefix,'file searched, no nuclear target',0
            name='coadd-'+prefix+'.fits'
            with sp.unpack(fetch('desi-pixels-'+prefix,directory+name,timeout=110)) as h:
                all_ex=h['EXP_FIBERMAP'].data
                separate=[]
                for tid,sep in match:
                    z=ztab.get(tid,ZREF)
                    ex=all_ex[all_ex['TARGETID']==tid]
                    mjd=min(ex['MJD']) if len(ex) else ''
                    end=max(ex['MJD']) if len(ex) else ''
                    store(svc,f'{name}-{tid}',*sp.desi_arrays(h,tid),z,mjd,end,sp.STD,sep,
                          f'{len(ex)} exposure records; z from DESI unflagged fit or BOSS reference')
                    for r in ex:
                        event(svc,'optical exposure in coadd',float(r['MJD']),kind='exposure',source=str(r['EXPID']),
                              note=f"{name}, target {tid}, tile {r['TILEID']}")
                    if len(set(ex['NIGHT']))>1:
                        separate.extend([(tid,z,sep,*a) for a in sorted(set((int(v['NIGHT']),int(v['TILEID']),int(v['PETAL_LOC'])) for v in ex))])
            for tid,z,sep,night,tile,petal in separate:
                n=f'coadd-{petal}-{tile}-{night}.fits'
                try:
                    with sp.unpack(fetch('desi-night-'+n,f'{root}/tiles/pernight/{tile}/{night}/{n}',timeout=65)) as h:
                        ex=h['EXP_FIBERMAP'].data
                        ex=ex[ex['TARGETID']==tid]
                        store(svc,n+f'-{tid}',*sp.desi_arrays(h,tid),z,min(ex['MJD']),max(ex['MJD']),
                              sp.STD,sep,'Per-night product overlapping HEALPix coadd.')
                except Exception as e:
                    status(svc,'unavailable',f'{n}: {e}')
            return prefix,'located and analyzed target',len(match)
        except Exception as e:
            return prefix,'unavailable '+str(e)[:180],0
    with ThreadPoolExecutor(max_workers=5) as pool:
        results=list(pool.map(process,jobs))
    for prefix,state,count in results:
        INDEX.append(dict(archive=svc,table=prefix,state=state,records=count))
    n=sum(r[2] for r in results)
    searched=sum('file searched' in r[1] for r in results)
    failures=sum(r[1].startswith('unavailable') for r in results)
    status(svc,'measured' if n else 'coverage_only' if failures or not searched else 'no_match',
           f'NSIDE 64 nested pixel(s) {pixels}; {len(results)} survey/program/pixel combinations; {searched} catalogs with zero matches; {failures} unreadable, {n} nuclear targets.')


def votables(key,url,params):
    raw=fetch(key,url,params,timeout=35)
    tree=BeautifulSoup(raw,'xml')
    errors=[i.get_text(' ',strip=True) for i in tree.find_all('INFO') if str(i.get('value','')).upper()=='ERROR']
    if errors:
        raise ValueError(' '.join(errors)[:300])
    (OUT/(key+'.xml')).write_bytes(raw) if len(raw)<180000 else None
    tables=list(parse(io.BytesIO(raw)).iter_tables())
    data=[]
    for tab in tables:
        # The service may put a display label in FIELD.name and a column name in FIELD.ID.
        for row in tab.array:
            record={}
            for i,field in enumerate(tab.fields):
                val=row[i]
                if np.ma.is_masked(val):
                    continue
                val=str(val.decode(errors='replace') if isinstance(val,bytes) else val).strip()
                for k in (field.ID,field.name,field.ucd):
                    if k:
                        record.setdefault(re.sub(r'[^a-z0-9]','',k.lower()),val)
            data.append(record)
    return tables,data


def attempt_lamost(svc,table,row):
    def v(*keys):
        return next((row[k] for k in keys if row.get(k)), '')
    rr=v('nucra','raj2000','ra','posEqRa'.lower(),'poseqrametamain')
    dd=v('nucdec','dej2000','dec','poseqdec','poseqdecmetamain')
    sep=sp.dist(rr,dd) if rr and dd else np.nan
    if np.isfinite(sep) and sep>3.:
        return
    obsid=v('obsid','obsid1','obsid2','designation','ssaaccref','accref')
    if not obsid:
        obsid='row'+str(len(META))
    urls=[v('accessurl','accref','ssaaccref','metarefurl','ref','url')]
    if obsid.isdigit():
        urls.extend([f'https://www.lamost.org/dr5/v3/spectrum/fits/{obsid}',
                     f'https://www.lamost.org/dr7/v2.0/spectrum/fits/{obsid}',
                     f'https://www.lamost.org/dr9/v2.0/spectrum/fits/{obsid}'])
    urls=[url.replace('http://','https://',1) if url.startswith('http') else urljoin('https://dc.g-vo.org/getproduct/',url)
          for url in urls if url]
    if not urls:
        status(svc,'coverage_only',f'{table} matched nuclear row lacks a downloadable FITS identifier')
        return
    for i,url in enumerate(urls):
        try:
            with sp.unpack(fetch('lamost-pixels-'+sp.keyname(table+obsid)+f'-{i}',url,timeout=25)) as h:
                w,f,iv,unit=sp.flux_array(h)
                z=sp.num(v('z','redshift','srcRedshift'.lower()))
                if not np.isfinite(z) or not 0<z<.3:
                    z=ZREF
                mjd=sp.num(v('mjd','mjdobs','mjdate','ssadatamidpoint'))
                if not 40000<mjd<70000:
                    mjd=sp.num(h[0].header.get('MJD-OBS',h[0].header.get('MJD',np.nan)))
                if not np.isfinite(mjd):
                    mjd=sp.asdate(h[0].header.get('DATE-OBS',v('dateobs','obsdate')))
                name=sp.keyname('LAMOST-'+obsid)
                if any(r['file']==name+'.npz' for r in META):
                    return
                store(svc,name,w,f,iv,z,mjd,'',unit,sep,table+'; LAMOST continuum scale can be relative')
            return
        except Exception as e:
            last=e
    status(svc,'unavailable',f'{table}, {obsid}: {last}')


def lamost_direct():
    svc='LAMOST GAVO VOTable'
    for table,ra,dec in [('lamost5.data','raj2000','dej2000'),('lamost6.ssa_lrs','ra','dec'),('lamost6.ssa_mrs','ra','dec')]:
        try:
            q=f'SELECT TOP 100 {ra} AS nuc_ra,{dec} AS nuc_dec,t.* FROM {table} AS t WHERE '+sp.box(ra,dec)
            _,hits=votables('lamost-return-'+sp.keyname(table),'https://dc.g-vo.org/tap/sync',
                            dict(REQUEST='doQuery',LANG='ADQL',FORMAT='votable',QUERY=q))
            INDEX.append(dict(archive=svc,table=table,state='coordinate query',records=len(hits)))
            status(svc,'coverage_only' if hits else 'no_match',f'{table}: {len(hits)} rows in 10 arcsec bounding box')
            for row in hits:
                attempt_lamost(svc,table,row)
        except Exception as e:
            status(svc,'unavailable',f'{table}: {e}')
    url='https://vizier.cfa.harvard.edu/viz-bin/votable'
    for src in ('V/164/dr5','VII/249'):
        label='6dF Galaxy Survey US mirror' if src.startswith('VII/') else 'LAMOST DR5 US mirror'
        try:
            ts,hits=votables('catalog-US-'+sp.keyname(src),url,{'-source':src,'-c':f'{HOST_RA} {HOST_DEC}','-c.rs':'3','-out.all':'1','-out.max':'80'})
            if not ts:
                raise ValueError('no returned catalog table')
            status(label,'coverage_only' if hits else 'no_match',f'{len(hits)} catalog rows in 3 arcsec')
            if src.startswith('V/'):
                for row in hits:
                    attempt_lamost(label,src,row)
        except Exception as e:
            status(label,'unavailable',str(e))


def template(wave,flux,ivar,z,line):
    """Fixed-shape broad component with free local continuum and tied narrow doublet."""
    hb=line=='Hbeta'
    center=4862.68 if hb else 6564.61
    lo,hi=(4780.,5075.) if hb else (6425.,6695.)
    main,weak,ratio=(5008.24,4960.30,2.98) if hb else (6585.27,6549.86,2.96)
    good=np.isfinite(wave)&np.isfinite(flux)
    if np.any(ivar>0):
        good &= np.isfinite(ivar)&(ivar>0)
    x=wave/(1+z)
    good &= (x>lo)&(x<hi)
    x,y,iv=x[good],flux[good],ivar[good]
    if len(x)<60 or not all(np.any(abs(x-t)<12.) for t in (lo+12,center,hi-12)):
        return {}
    order=np.argsort(x)
    x,y,iv=x[order],y[order],iv[order]
    if not np.any(iv>0):
        noise=np.median(abs(np.diff(y)-np.median(np.diff(y))))*1.4826/np.sqrt(2)
        if not noise>0:
            return {}
        iv=np.full(len(x),noise**-2)
    def g(mu,dv,sig):
        return np.exp(-.5*((x-mu*(1+dv/299792.458))/(mu*sig/299792.458))**2)
    weights=np.sqrt(iv)
    best=None
    # A fixed 5,000 km/s template is specified by the measured BOSS H-alpha widths.
    # The centered template gives an interpretable sensitivity even if H-beta is blended.
    bsig=5000/2.35482
    broad=g(center,0,bsig)
    for vel in (-160,-80,0,80,160):
        for width in (80,120,180,250,320):
            A=np.column_stack([np.ones_like(x),(x-center)/100,
                               g(center,vel,width),g(main,vel,width)+g(weak,vel,width)/ratio,broad])
            B=A*weights[:,None]
            coeff,_,rank,_=lstsq(B,y*weights)
            if rank!=5:
                continue
            chi=np.sum(((y-A@coeff)*weights)**2)
            if best is None or chi<best[0]:
                best=(chi,coeff,B)
    if best is None:
        return {}
    chi,coef,B=best
    err=np.sqrt(np.linalg.pinv(B.T@B)[-1,-1]*max(1,chi/(len(x)-5)))
    norm=np.sqrt(2*np.pi)*center*bsig/299792.458*(1+z)
    f=coef[-1]*norm
    e=err*norm
    return dict(template_flux=f'{f:.7g}',error=f'{e:.7g}',sensitivity_3sigma=f'{3*e:.7g}',
                upper95=f'{max(0.,f)+1.645*e:.7g}',snr=f'{f/e:.5g}' if e>0 else '',
                reduced_chi2=f'{chi/(len(x)-5):.4g}')


def fits_from_pixels():
    meta=json.loads((PIX/'index.json').read_text())
    results=[]
    for m in meta:
        with np.load(PIX/m['file']) as h:
            w,f,iv=(h[k] for k in ('wave','flux','ivar'))
            for line in ('Halpha','Hbeta'):
                out=template(w,f,iv,m['z'],line)
                results.append(dict(product=m['product'],mjd=m['mjd'],line=line,**out))
    return results


def comparison():
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(8.3,5.8),constrained_layout=True)
    for m in sorted(META,key=lambda v:v['mjd'] or 0):
        with np.load(PIX/m['file']) as h:
            w,f,iv=(h[k] for k in ('wave','flux','ivar'))
            x=w/(1+m['z'])
            for ax,(lo,hi) in zip(axes,[(4760,5080),(6420,6710)]):
                good=(x>lo)&(x<hi)&np.isfinite(f)
                if np.any(iv>0):
                    good &= iv>0
                if sum(good)<30:
                    continue
                baseline=good & ((x<lo+35)|(x>hi-35))
                if not baseline.any() or not np.median(f[baseline])>0:
                    continue
                ii=np.argsort(x[good])
                ax.plot(x[good][ii],f[good][ii]/np.median(f[baseline]),lw=.8,alpha=.8,
                        label=f'{m["archive"]} MJD {float(m["mjd"]):.0f}' if m['mjd']!='' else m['archive'])
    for ax in axes:
        ax.set_ylabel('Flux / local continuum')
        if ax.lines:
            ax.legend(frameon=False,fontsize=8)
    axes[-1].set_xlabel('Rest wavelength (Å)')
    for x,label,ax in [(4862.68,'Hβ',axes[0]),(5008.24,'[O III]',axes[0]),(6564.61,'Hα',axes[1])]:
        ax.axvline(x,color='.6',lw=.65,ls='--')
        ax.annotate(label,xy=(x,.90),xycoords=('data','axes fraction'),fontsize=8)
    fig.savefig(OUT/'balmer-comparison.png',dpi=150)
    plt.close(fig)


def write_results():
    PIX.mkdir(exist_ok=True)
    (PIX/'index.json').write_text(json.dumps(sorted(META,key=lambda m:(m['mjd'] or 0,m['product'])),indent=2)+'\n')
    (OUT/'spectral-search-pass2.json').write_text(json.dumps(sorted(INDEX,key=lambda r:r['table']),indent=2)+'\n')
    rs=fits_from_pixels()
    free=[]
    for m in META:
        with np.load(PIX/m['file']) as pix:
            w,f,iv=(pix[k] for k in ('wave','flux','ivar'))
            use=np.isfinite(w)&np.isfinite(f)
            has_ivar=np.any(iv>0)
            if has_ivar:
                use &= np.isfinite(iv)&(iv>0)
            w,f,iv=w[use],f[use],iv[use]
            order=np.argsort(w)
            w,f,iv=w[order],f[order],iv[order]
            err=1/np.sqrt(iv) if has_ivar else np.full(w.shape,np.nan)
            for line in ('Halpha','Hbeta'):
                ans=sp.fit_complex(w,f,err,m['z'],line)
                free.append(dict(product=m['product'],mjd=m['mjd'],end_mjd=m['end_mjd'],component=line,
                                 unit=m['unit'],**{k:v for k,v in ans.items() if k!='component'}))
            for c in (5100,6200):
                sub=abs(w/(1+m['z'])-c)<25.
                if sum(sub)>=8:
                    value=float(np.median(f[sub]))
                    free.append(dict(product=m['product'],mjd=m['mjd'],end_mjd=m['end_mjd'],
                                     unit=m['unit'],component=f'C{c}',continuum=f'{value:.7g}'))
                    if m['archive']!='BOSS DR17':
                        event(m['archive'],f'optical continuum {c} A rest',m['mjd'],f'{value:.7g}',
                              unit=m['unit'],source=m['product'],end_mjd=m['end_mjd'],note='50-A rest window; includes host')
    fc='product mjd end_mjd component unit continuum fit broad_flux broad_flux_err fwhm_kms delta_bic snr rest_ew_A forbidden_flux'.split()
    with (OUT/'free-profile-balmer.csv').open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=fc)
        writer.writeheader()
        writer.writerows(free)
    cols='product mjd line template_flux error sensitivity_3sigma upper95 snr reduced_chi2'.split()
    with (OUT/'fixed-width-balmer.csv').open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=cols)
        writer.writeheader()
        writer.writerows(rs)
    with (OUT/'spectra-pass2.md').open('w') as fh:
        fh.write('# Spectral comparison and sensitivity\n\n')
        fh.write('| Product | MJD | Rest-frame feature | Fixed-width flux ± formal error | 3σ sensitivity |\n|---|---:|---|---|---:|\n')
        for r in rs:
            fh.write(f"| {r['product']} | {r['mjd']} | {r['line']} | {r.get('template_flux','—')} ± {r.get('error','—')} | {r.get('sensitivity_3sigma','—')} |\n")
        fh.write('\nThe common template has FWHM 5,000 km/s, centered on each vacuum Balmer transition. It uses each native observed flux density unit × observed Å. BOSS and DESI use 10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹. LAMOST flux scales may be relative. Both BOSS spectra precede the 2014 slew. Coverage and dates for any later spectrum appear in the pixel index.\n\n')
        fh.write('The continuum is linear, with a shared narrow velocity and width and a fixed doublet ratio for [O III] or [N II]. The error comes from weighted regression at the best narrow parameters, inflated by residual scatter. The 3σ sensitivity is the formal error times three. The upper95 column adds 1.645 errors to the positive part of the fitted amplitude. Both are conditional on this profile, the continuum model and the recorded aperture. The Hβ template can absorb Fe II, stellar absorption and mismodeled narrow lines. Use the free-profile fit and spectra as well when interpreting changes. The continuum-normalized plot permits visual line-shape comparison across different flux units.\n\n')
        fh.write('Pixel arrays, masks expressed as zero inverse variance, flux units, input provenance and observation dates are deposited. Run `python spec_extra.py --replay` from the study folder with the pinned dependencies. This repeats the fixed-width measurements without querying archives. The original full FITS bytes are identified by SHA-256 in the receipts.\n')
    if META:
        comparison()


def run():
    jobs=[('BOSS pixel records',boss),('DESI direct products',desi_direct),('LAMOST and 6dF spectra',lamost_direct)]
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        for name,job in [(n,pool.submit(fn)) for n,fn in jobs]:
            try:
                job.result()
            except Exception as exc:
                status(name,'unavailable',str(exc))
    write_results()


if __name__=='__main__' and sys.argv[1:]==['--replay']:
    actual=fits_from_pixels()
    with (OUT/'fixed-width-balmer.csv').open() as fh:
        expected=list(csv.DictReader(fh))
    if len(actual)!=len(expected):
        raise SystemExit('pixel replay count differs')
    for a,b in zip(actual,expected):
        if (a['product'],a['line'])!=(b['product'],b['line']):
            raise SystemExit('pixel replay product differs')
        for k in ('template_flux','error','sensitivity_3sigma','upper95'):
            if k in a and not np.isclose(float(a[k]),float(b[k]),rtol=3e-5):
                raise SystemExit('pixel replay fit differs')
    print(f'Reproduced {len(actual)} Balmer template rows from deposited pixels.')
