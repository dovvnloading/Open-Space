"""Public spectra at the 2MASX nucleus. Run from the follow-up driver."""
import csv
import gzip
import io
import math
import re
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urljoin

import numpy as np
from astropy.io import fits
from astropy.io.votable import parse
from astropy.time import Time
from bs4 import BeautifulSoup
from scipy.optimize import least_squares

from common import fetch, event, status, OUT, HOST_RA, HOST_DEC

records = []
RADIUS = 10 / 3600
DRA = RADIUS / math.cos(math.radians(HOST_DEC))
STD = '10^-17 erg cm^-2 s^-1 Angstrom^-1'


def dist(ra, dec):
    return 3600 * math.hypot((float(ra)-HOST_RA)*math.cos(math.radians(HOST_DEC)), float(dec)-HOST_DEC)


def tap(key, base, query):
    raw = fetch(key, base, dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=query))
    return rows(raw)


def rows(raw):
    txt = raw.decode('utf-8-sig', errors='replace').strip()
    if txt.startswith('#Table'):
        txt = '\n'.join(txt.splitlines()[1:])
    if not txt or txt.startswith('<') or txt.lower().startswith(('error', 'exception', 'message:')):
        raise ValueError('unexpected archive response: '+txt[:180])
    return [{k.lower().strip(): v for k,v in r.items() if k} for r in csv.DictReader(io.StringIO(txt))]


def box(ra, dec):
    return (f'{ra} BETWEEN {HOST_RA-DRA:.9f} AND {HOST_RA+DRA:.9f} AND '
            f'{dec} BETWEEN {HOST_DEC-RADIUS:.9f} AND {HOST_DEC+RADIUS:.9f}')


def keyname(s):
    return re.sub('[^A-Za-z0-9.-]+','-',str(s))[:96]


def asdate(x):
    try:
        return float(Time(str(x).strip(),scale='utc').mjd)
    except Exception:
        return ''


def num(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return np.nan


def compact(x):
    return f'{x:.5g}' if np.isfinite(x) else ''


def fit_complex(w, flux, err, z, which):
    """Fit a shared narrow-line kinematic model, then add one broad Balmer line."""
    c = 299792.458
    hb = which == 'Hbeta'
    rest = w / (1+z)
    lo, hi, center = (4730,5100,4862.68) if hb else (6400,6710,6564.61)
    ok = (rest>=lo) & (rest<=hi)
    if ok.sum()<60 or not all(np.any(abs(rest-t)<12) for t in (lo+15, center, hi-15)):
        return {}
    x, y, s = rest[ok], flux[ok], err[ok]
    scaling = np.nanmedian(abs(y)) or 1.
    y, s = y/scaling, s/scaling
    # Bound the noise floor by local pixel-to-pixel scatter for spectra without inverse variance.
    empirical = 1.4826*np.median(abs(np.diff(y)-np.median(np.diff(y))))/np.sqrt(2)
    s = np.where(np.isfinite(s) & (s>0), s, max(empirical,1.e-5))
    companion = (5008.24,4960.30,2.98) if hb else (6585.27,6549.86,2.96)
    def profile(w0, dv, sig):
        return np.exp(-0.5*((x-w0*(1+dv/c))/(sig*w0/c))**2)
    def residual(p):
        a,b,v,width,balmer,forbidden=p[:6]
        model=a+b*(x-center)/100 + balmer*profile(center,v,width)
        model+=forbidden*(profile(companion[0],v,width)+profile(companion[1],v,width)/companion[2])
        if len(p)>6:
            model+=p[6]*profile(center,p[7],p[8])
        return (model-y)/s
    init = [np.median(y),0.,0.,150.,max(.01,np.max(y)-np.median(y)),.5]
    low=[-np.inf,-np.inf,-500.,40.,0.,0.]
    high=[np.inf,np.inf,500.,450.,np.inf,np.inf]
    null=least_squares(residual,init,bounds=(low,high),max_nfev=450)
    candidates=[least_squares(residual,[*null.x,.3,0.,width],
                              bounds=(low+[0.,-1500.,510.],high+[np.inf,1500.,6000.]),max_nfev=650)
                for width in (800.,2100.)]
    alt=min(candidates,key=lambda a:a.cost)
    amp,dv,sig=alt.x[6:]
    cov=np.linalg.pinv(alt.jac.T@alt.jac)*max(1.,2*alt.cost/max(1,len(x)-9))
    norm=np.sqrt(2*np.pi)*center/c*(1+z)*scaling
    area=norm*amp*sig
    var=norm**2*(sig**2*cov[6,6]+amp**2*cov[8,8]+2*amp*sig*cov[6,8])
    area_err=np.sqrt(max(0.,var))
    sn=area/area_err if area_err>0 else np.nan
    dbic=2*(null.cost-alt.cost)-3*np.log(len(x))
    flags='broad favored' if (alt.success and dbic>=10 and sn>=5 and 525<sig<5850
                               and abs(dv)<1450) else 'inconclusive'
    forbidden=alt.x[5]*alt.x[3]*np.sqrt(2*np.pi)*companion[0]/c*(1+z)*scaling
    return dict(component=which,fit=flags,broad_flux=compact(area),broad_flux_err=compact(area_err),
                fwhm_kms=compact(2.35482*sig),delta_bic=compact(dbic),snr=compact(sn),
                rest_ew_A=compact(area/(alt.x[0]*scaling*(1+z)) if alt.x[0]>0 else np.nan),
                forbidden_flux=compact(forbidden))


def measure(archive, ident, w, f, ivar, z, mjd, end, unit, sep, note=''):
    base=dict(archive=archive,product=ident,mjd=mjd,end_mjd=end,z=z,separation_arcsec=compact(sep),unit=unit,note=note)
    w=np.ravel(np.asarray(w,dtype=float))
    f=np.ravel(np.asarray(f,dtype=float))
    good=np.isfinite(w)&np.isfinite(f)&(w>3200)&(w<11500)
    if ivar is not None:
        ivar=np.ravel(np.asarray(ivar,dtype=float))
        good &= np.isfinite(ivar)&(ivar>0)
    order=np.argsort(w[good])
    err= np.full(sum(good),np.nan) if ivar is None else 1/np.sqrt(ivar[good][order])
    w,f=w[good][order],f[good][order]
    if len(w)<50 or not (0<z<2.):
        records.append(dict(base,component='spectrum',fit='insufficient coverage or redshift'))
        return
    measured=[]
    for c in (5100,6200):
        mask=abs(w/(1+z)-c)<25
        if mask.sum()>=8:
            val=float(np.median(f[mask]))
            records.append(dict(base,component=f'C{c}',continuum=compact(val),fit='median flux density'))
            event(archive,f'optical continuum {c} A rest',mjd,compact(val),unit=unit,source=ident,
                  note=f'50-A rest window, observed flux density; {note}',end_mjd=end)
    for line in ('Hbeta','Halpha'):
        fit=fit_complex(w,f,err,z,line)
        if not fit:
            records.append(dict(base,component=line,fit='insufficient wavelength coverage'))
            continue
        records.append(dict(base,**fit))
        measured.append(line+' '+fit['fit'])
    event(archive,'optical spectrum',mjd,kind='spectrum',source=ident,end_mjd=end,note='; '.join(measured+[note]))


def unpack(payload):
    if payload.startswith(b'\x1f\x8b'):
        payload=gzip.decompress(payload)
    return fits.open(io.BytesIO(payload),memmap=False)


def sdss():
    service='SDSS Legacy/BOSS/eBOSS'
    url='https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch'
    query=('SELECT s.specObjID,s.plate,s.mjd,s.fiberID,s.run2d,s.ra,s.dec,s.z,s.zWarning,s.survey '
           'FROM SpecObjAll s WHERE '+box('s.ra','s.dec'))
    response=fetch('spectra-sdss-dr17-all',url,dict(cmd=query,format='csv'))
    if b'specobjid' not in response.lower():
        raise ValueError('SDSS did not return the requested schema: '+response[:100].decode(errors='replace'))
    data=rows(response)
    near=[r for r in data if dist(r['ra'],r['dec'])<=RADIUS*3600]
    done=0
    for r in near:
        sep=dist(r['ra'],r['dec'])
        if sep>3.:
            status(service,'nearby',f"fiber {r['specobjid']} lies {sep:.2f} arcsec from host; excluded")
            continue
        p,m,fi,redux=int(r['plate']),int(r['mjd']),int(r['fiberid']),r['run2d'].strip()
        name=f'spec-{p:04d}-{m:05d}-{fi:04d}.fits'
        paths=(['sdss','eboss'] if redux.isdigit() else ['eboss','sdss'])
        for survey,layout in [(survey,layout) for survey in paths for layout in ('spectra/lite','spectra/full','spectra')]:
            path=f'https://data.sdss.org/sas/dr17/{survey}/spectro/redux/{redux}/{layout}/{p:04d}/{name}'
            try:
                with unpack(fetch('sdss-'+name+'-'+survey+'-'+keyname(layout),path,timeout=18)) as h:
                    t=h['COADD'].data
                    iv=t['ivar'].astype(float)
                    if 'and_mask' in [s.lower() for s in t.names]:
                        iv[t['and_mask']!=0]=0
                    if int(r['zwarning']) and not np.isfinite(reference_z):
                        raise ValueError('flagged redshift and unavailable independent redshift')
                    z=float(r['z']) if not int(r['zwarning']) else reference_z
                    measure(service,name,10**t['loglam'],t['flux'],iv,z,m,'',STD,sep,
                            f"run2d {redux}, fiber {fi}, survey {r['survey']}, zWarning {r['zwarning']}")
                done+=1
                break
            except Exception as e:
                last=e
        else:
            status(service,'unavailable',f'{name}: {last}')
    status(service,'measured' if done else 'coverage_only' if near else 'no_match',
           f'{len(near)} SpecObjAll rows in 10 arcsec, {done} analyzed within 3 arcsec')


def sdss_alternates():
    service='SDSS independent VizieR/FITS cross-check'
    sources=('J/ApJS/267/19','V/154/sdss16','V/147/sdss12')
    candidates={}
    for src in sources:
        try:
            raw=fetch('sdss-id-'+keyname(src),'https://vizier.cds.unistra.fr/viz-bin/votable',
                      {'-source':src,'-c':f'{HOST_RA} {HOST_DEC}','-c.rs':'3','-out.all':'','-out.max':'100'})
            for tab in parse(io.BytesIO(raw)).iter_tables():
                for row in tab.to_table(use_names_over_ids=True):
                    vals={re.sub(r'[^a-z0-9]','',k.lower()):str(row[k]).strip() for k in row.colnames}
                    def value(*keys):
                        return next((vals[k] for k in keys if vals.get(k) not in ('','--',None)),'')
                    p=value('plate','plateid','plateid1','plate1','nplate')
                    m=value('mjd','smjd','specmjd','mjd1')
                    f=value('fiber','fiberid','fiber1','fib')
                    if not (p.isdigit() and m.isdigit() and f.isdigit()):
                        # The DR8+ specObjID integer packs plate, fiber and (MJD-50000).
                        sid=value('spobjid','specobjid','spid','specid')
                        if sid.isdigit() and int(sid)>2**50:
                            sid=int(sid)
                            p,m,f=str(sid>>50),str(((sid>>24)&16383)+50000),str((sid>>38)&4095)
                        else:
                            triplet=re.search(r'(?<![0-9])([0-9]{3,5})[- /]([5-6][0-9]{4})[- /]([0-9]{3,4})(?![0-9])',' '.join(vals.values()))
                            if not triplet:
                                continue
                            p,m,f=triplet.groups()
                    if 0<int(p)<20000 and 48000<int(m)<70000 and 0<int(f)<=1500:
                        candidates.setdefault((int(p),int(m),int(f)),(value('zsp','zspec','z'),src))
            status(service,'queried',f'{src} searched within 3 arcsec')
        except Exception as exc:
            status(service,'unavailable',f'{src}: {exc}')
    # This path is explicitly independent of SkyServer SQL. Do not refit its successful spectra.
    for (p,m,fi),(zcat,src) in sorted(candidates.items()):
        name=f'spec-{p:04d}-{m:05d}-{fi:04d}.fits'
        if any(r['product']==name for r in records):
            continue
        versions=[('sdss','26'),('sdss','103'),('sdss','104'),('eboss','v5_13_2')] if m<55500 else [('eboss','v5_13_2'),('sdss','26')]
        success=False
        for sv,v in versions:
            for path in ('spectra/lite','spectra/full','spectra'):
                url=f'https://data.sdss.org/sas/dr17/{sv}/spectro/redux/{v}/{path}/{p:04d}/{name}'
                try:
                    with unpack(fetch(f'sdss-alt-{p}-{m}-{fi}-{v}-'+keyname(path),url,timeout=12)) as h:
                        t=h['COADD'].data
                        z=num(zcat)
                        for ext in h:
                            names={k.lower():k for k in (getattr(ext.data,'names',None) or [])}
                            if 'z' in names and len(ext.data)==1:
                                z=num(ext.data[names['z']][0])
                                break
                        if not np.isfinite(z) or z<=0:
                            z=reference_z
                        iv=t['ivar'].astype(float)
                        if 'and_mask' in [k.lower() for k in t.names]:
                            iv[t['and_mask']!=0]=0
                        measure(service,name,10**t['loglam'],t['flux'],iv,z,m,'',STD,np.nan,
                                f'VizieR catalog {src}, 3 arcsec cone, run2d {v}')
                    status(service,'measured',name)
                    success=True
                    break
                except Exception as e:
                    last=e
            if success:
                break
        if not success:
            status(service,'unavailable',f'{name}, tried standard DR17 layouts: {last}')
    if not candidates:
        status(service,'coverage_only','Catalog cone searches yielded no usable plate-MJD-fiber identifier')


def other_optical():
    for service,base in [('ESO optical spectra','https://archive.eso.org/tap_obs/sync'),
                         ('MAST optical spectra','https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync'),
                         ('CADC optical spectra','https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap/sync')]:
        try:
            q=("SELECT TOP 150 obs_id,s_ra,s_dec,t_min,t_max,access_url,access_format "
               "FROM ivoa.obscore WHERE dataproduct_type='spectrum' AND em_min<8e-7 AND em_max>4e-7 AND "+box('s_ra','s_dec'))
            hits=tap(keyname(service),base,q)
            good=[r for r in hits if dist(r['s_ra'],r['s_dec'])<=3.]
            for r in good:
                key=service.split()[0]+'-'+keyname(r['obs_id'])
                if r['access_url'].startswith('http') and 'fits' in r['access_format'].lower():
                    try:
                        measure_other(service,key,r['access_url'],sep=dist(r['s_ra'],r['s_dec']),
                                      mjd=r['t_min'],note='ObsCore associated by center position')
                    except Exception as e:
                        status(service,'unavailable',f'{key}: {e}')
                else:
                    event(service,'optical spectrum metadata',r['t_min'],kind='coverage',source=key,
                          note='No direct FITS access URL',end_mjd=r['t_max'])
            status(service,'coverage_only' if good else 'no_match',f'{len(good)} nuclear 1D optical products by ObsCore center')
        except Exception as e:
            status(service,'unavailable',str(e))


def desi_arrays(h, target):
    ids=h['FIBERMAP'].data['TARGETID']
    i=int(np.flatnonzero(ids==target)[0])
    chunks=[]
    for band in ('B','R','Z'):
        w=h[band+'_WAVELENGTH'].data
        f=h[band+'_FLUX'].data[i]
        iv=h[band+'_IVAR'].data[i].copy()
        if band+'_MASK' in h:
            iv[h[band+'_MASK'].data[i]!=0]=0
        # DESI cameras overlap. Assign each wavelength to just one camera.
        mask=((w<5750) if band=='B' else ((w>=5750)&(w<7600)) if band=='R' else (w>=7600))
        chunks.append((w[mask],f[mask],iv[mask]))
    return [np.concatenate(a) for a in zip(*chunks)]


def desi():
    service='DESI DR1'
    base='https://datalab.noirlab.edu/tap/sync'
    query='SELECT TOP 200 targetid,target_ra,target_dec,survey,program,healpix,z,zwarn FROM desi_dr1.zpix WHERE '+box('target_ra','target_dec')
    data=tap('desi-dr1-zpix',base,query)
    done=0
    for row in data:
        sep=dist(row['target_ra'],row['target_dec'])
        if sep>3.:
            status(service,'nearby',f"{row['targetid']} is {sep:.2f} arcsec from host; excluded")
            continue
        s,p,hp,tid=row['survey'],row['program'],int(row['healpix']),int(row['targetid'])
        root='https://data.desi.lbl.gov/public/dr1/spectro/redux/iron'
        name=f'coadd-{s}-{p}-{hp}.fits'
        url=f'{root}/healpix/{s}/{p}/{hp//100}/{hp}/{name}'
        with unpack(fetch('desi-'+name,url,timeout=115)) as h:
            target_exps=h['EXP_FIBERMAP'].data
            exps=target_exps[target_exps['TARGETID']==tid]
            start=min(exps['MJD']) if len(exps) else ''
            end=max(exps['MJD']) if len(exps) else ''
            z=float(row['z']) if int(row['zwarn'])==0 else reference_z
            measure(service,name+f'#{tid}',*desi_arrays(h,tid),z,start,end,STD,sep,
                    f'HEALPix coadd of {len(exps)} exposure records; survey {s}/{p}')
            done+=1
            groups=sorted(set((int(v['NIGHT']),int(v['TILEID']),int(v['PETAL_LOC'])) for v in exps))
        # Separate nights are informative only for a stack spanning multiple nights.
        if len({v[0] for v in groups})>1:
            for night,tile,petal in groups:
                name=f'coadd-{petal}-{tile}-{night}.fits'
                url=f'{root}/tiles/pernight/{tile}/{night}/{name}'
                try:
                    with unpack(fetch('desi-'+name,url,timeout=110)) as h:
                        ex=h['EXP_FIBERMAP'].data
                        ex=ex[ex['TARGETID']==tid]
                        measure(service,name+f'#{tid}',*desi_arrays(h,tid),z,min(ex['MJD']),max(ex['MJD']),STD,sep,
                                'per-night tile coadd; overlaps HEALPix coadd')
                except Exception as e:
                    status(service,'unavailable',f'Per-night product {name}: {e}')
    status(service,'measured' if done else 'coverage_only' if data else 'no_match',
           f'{len(data)} zpix entries in search box, {done} nuclear coadds analyzed')


def flux_array(h):
    """Read common optical FITS conventions. Reject unknown wavelength encodings."""
    for ext in h:
        data=ext.data
        names={x.lower():x for x in (getattr(data,'names',[]) or [])}
        wk=next((x for x in ('loglam','wavelength','wave','lambda') if x in names),None)
        if wk and 'flux' in names:
            wave=np.asarray(data[names[wk]],dtype=float).ravel()
            iunit=str(ext.columns[names[wk]].unit or '').lower()
            if wk=='loglam':
                wave=10**wave
            elif iunit in ('nm','nanometer','nanometers'):
                wave*=10.
            elif iunit in ('um','micron','microns'):
                wave*=10000.
            elif iunit not in ('','a','angstrom','angstroms','angstroem','aa'):
                continue
            f=np.asarray(data[names['flux']]).ravel()
            iv=np.asarray(data[names['ivar']]).ravel() if 'ivar' in names else None
            return wave,f,iv,str(ext.columns[names['flux']].unit or h[0].header.get('BUNIT','native units'))
    a=h[0].header
    f=np.squeeze(h[0].data)
    if 'COEFF0' in a and 'COEFF1' in a:
        # LAMOST LRS primary image planes: flux, inverse variance, wavelength, masks.
        w=10**(a['COEFF0']+a['COEFF1']*np.arange(f.shape[-1]))
        return w,(f[0] if f.ndim>1 else f),(f[1] if f.ndim>1 and len(f)>=5 else None),a.get('BUNIT','native units')
    if f.ndim!=1 or 'CRVAL1' not in a or ('CD1_1' not in a and 'CDELT1' not in a):
        raise ValueError('unrecognized optical dispersion or image format')
    w=a['CRVAL1']+(np.arange(len(f))+1-a.get('CRPIX1',1))*a.get('CD1_1',a.get('CDELT1'))
    if a.get('DC-FLAG',0)==1 or 'LOG' in str(a.get('CTYPE1','')).upper():
        w=10**w
    if 'nm'==str(a.get('CUNIT1','')).strip().lower():
        w*=10
    return w,f,None,a.get('BUNIT','native units')


def measure_other(service,key,url,z=math.nan,sep=math.nan,mjd='',note=''):
    with unpack(fetch(key,url,timeout=65)) as h:
        head=h[0].header
        w,f,iv,unit=flux_array(h)
        if np.count_nonzero((w>3500)&(w<10500))<50:
            raise ValueError('no verified optical wavelength overlap')
        if mjd=='':
            mjd=head.get('MJD-OBS',head.get('MJD',asdate(head.get('DATE-OBS',''))))
        if not np.isfinite(z):
            z=num(head.get('Z',head.get('REDSHIFT',np.nan)))
        if not np.isfinite(z) or z<=0:
            z=reference_z
        measure(service,key,w,f,iv,z,mjd,'',unit,sep,note)


def lamost():
    svc='LAMOST (GAVO public releases)'
    base='https://dc.g-vo.org/tap/sync'
    data=tap('lamost-registry-gavo',base,
             "SELECT table_name,column_name FROM TAP_SCHEMA.columns WHERE table_name LIKE 'lamost%.%'")
    tables={}
    for r in data:
        tables.setdefault(r['table_name'],set()).add(r['column_name'].lower())
    searched=matched=0
    seen=set()
    for table,columns in sorted(tables.items(),reverse=True):
        def col(options):
            return next((c for c in options if c in columns),None)
        ra=col(('ra','raj2000','s_ra','alpha','raj2000j2'))
        dec=col(('dec','dej2000','s_dec','delta','dej2000j2'))
        loc=col(('ssa_location','s_region'))
        if ra and dec:
            where=box(ra,dec)
        elif loc=='ssa_location':
            where=f"1=CONTAINS(ssa_location,CIRCLE('ICRS',{HOST_RA},{HOST_DEC},{RADIUS}))"
        else:
            continue
        if not re.fullmatch(r'[\w.]+',table):
            continue
        try:
            found=tap('lamost-'+table,base,f'SELECT TOP 100 * FROM {table} WHERE {where}')
            searched+=1
        except Exception as e:
            status(svc,'unavailable',f'{table}: {e}')
            continue
        for row in found:
            sep=dist(row[ra],row[dec]) if ra and dec else np.nan
            if np.isfinite(sep) and sep>3.:
                continue
            matched+=1
            link=row.get('accref') or row.get('access_url','')
            if link and link not in seen:
                seen.add(link)
                if not link.startswith('http'):
                    link='https://dc.g-vo.org/getproduct/'+link.lstrip('/')
                try:
                    measure_other(svc,'lamost-'+keyname(row.get('obsid',link.split('/')[-1])),link,
                                  num(row.get('z')),sep,note=table)
                except Exception as e:
                    status(svc,'unavailable',f'{table} FITS: {e}')
    status(svc,'coverage_only' if matched else 'no_match' if searched else 'unavailable',
           f'{searched} GAVO release tables coordinate-queried, {matched} possible records. Release coverage is limited to served tables.')


def vizier():
    for label,cat in [('6dF Galaxy Survey','VII/249/6dfgs'),('LAMOST VizieR DR5','V/164/dr5')]:
        try:
            raw=fetch('spectra-'+keyname(label),'https://vizier.cds.unistra.fr/viz-bin/votable',
                      {'-source':cat,'-c':f'{HOST_RA} {HOST_DEC}','-c.rs':'10','-out.all':'','-out.max':'100'})
            soup=BeautifulSoup(raw,'xml')
            if any(i.get('value','').upper()=='ERROR' or 'error' in i.get('name','').lower() for i in soup.find_all('INFO')):
                raise ValueError('VizieR reports a query error')
            tables=list(parse(io.BytesIO(raw)).iter_tables())
            if not tables:
                if 'no astronomical object found' in raw.decode(errors='replace').lower():
                    status(label,'no_match','VizieR returned no object in 10 arcsec')
                    continue
                raise ValueError('VizieR returned no recognized catalog table')
            n=sum(len(t.array) for t in tables)
            status(label,'coverage_only' if n else 'no_match',f'{n} public catalog entries in 10 arcsec ({cat})')
            if label.startswith('LAMOST'):
                for t in tables:
                    for r in t.to_table(use_names_over_ids=True):
                        row={k.lower():r[k] for k in r.colnames}
                        obsid=row.get('obsid')
                        if obsid is None:
                            continue
                        sep=dist(row.get('raj2000'),row.get('dej2000'))
                        if sep<=3.:
                            u=f'https://www.lamost.org/dr5/v3/spectrum/fits/{obsid}'
                            measure_other(label,f'lamost-dr5-{obsid}',u,num(row.get('z')),sep,note=cat)
        except Exception as e:
            status(label,'unavailable',str(e))


def ned():
    service='NED contributed spectra'
    u='https://ned.ipac.caltech.edu/cgi-bin/NEDspectra'
    raw=fetch('ned-spectrum-index',u,dict(objname='KUG 1256+241',extend='multi',detail='0',preview='0',numpp='5000'))
    soup=BeautifulSoup(raw,'html.parser')
    links=[]
    for a in soup.find_all('a',href=True):
        href=a['href']
        if re.search(r'\.fits?(\.gz)?($|\?)',href,re.I):
            link=urljoin(u,href)
            if link not in links:
                links.append(link)
    done=0
    for i,link in enumerate(links):
        try:
            measure_other(service,'ned-'+keyname(link.split('/')[-1]),link,note='NED object-name match; may repeat original archive')
            done+=1
        except Exception as e:
            status(service,'unavailable',f'FITS {i+1}: {e}')
    txt=' '.join(soup.stripped_strings)
    no_match=bool(re.search(r'no (?:published |available |spectroscopic )?spectra|0 spectra|no data (?:were )?found',txt,re.I))
    status(service,'measured' if done else 'no_match' if not links and no_match else 'unavailable',
           f'{len(links)} linked FITS, {done} analyzed. '+('' if done or no_match else 'Unrecognized or unreadable spectrum listing.'))


def run():
    global reference_z
    reference_z=np.nan
    try:
        ref=tap('simbad-spec-redshift','https://simbad.cds.unistra.fr/simbad/sim-tap/sync',
                'SELECT TOP 10 main_id,ra,dec,rvz_redshift FROM basic WHERE '+box('ra','dec'))
        close=sorted((dist(r['ra'],r['dec']),num(r['rvz_redshift'])) for r in ref if num(r['rvz_redshift'])>0)
        if close and close[0][0]<=3:
            reference_z=close[0][1]
    except Exception as e:
        status('SIMBAD redshift','unavailable',str(e))
    def sdss_combined():
        try:
            sdss()
        except Exception as e:
            status('SDSS DR17','unavailable',str(e))
        sdss_alternates()
    jobs=[('SDSS DR17',sdss_combined),('Other optical archives',other_optical),('DESI DR1',desi),('LAMOST GAVO',lamost),('6dF/LAMOST VizieR',vizier),('NED spectra',ned)]
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures=[(name,pool.submit(job)) for name,job in jobs]
        for name,job in futures:
            try:
                job.result()
            except Exception as e:
                status(name,'unavailable',str(e))
    fields=('archive','product','mjd','end_mjd','z','separation_arcsec','component','fit',
            'continuum','unit','broad_flux','broad_flux_err','fwhm_kms','delta_bic','snr','rest_ew_A','forbidden_flux','note')
    ordered=sorted(records,key=lambda r:(num(r['mjd']) if str(r['mjd']) else 1e12,r['archive'],r['component']))
    with (OUT/'spectra.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(ordered)
    with (OUT/'spectra.md').open('w') as f:
        f.write('# Nuclear optical spectroscopy\n\n')
        if not ordered:
            f.write('No optical spectrum yielded a numerical line comparison in this run. Archive statuses distinguish empty searches from failures.\n\n')
        else:
            f.write('| Product | MJD | Feature | Measurement |\n|---|---:|---|---|\n')
            for r in ordered:
                if r.get('broad_flux'):
                    detail=f"{r['fit']}, broad FWHM {r['fwhm_kms']} km/s, ΔBIC {r['delta_bic']}, flux {r['broad_flux']} ± {r['broad_flux_err']} in native flux-density units × Å"
                else:
                    detail=(r['continuum']+' '+r['unit']) if r.get('continuum') else r['fit']
                f.write(f"| {r['product']} | {r['mjd']} | {r['component']} | {detail} |\n")
        f.write('\nSearches center on the 2MASX nucleus. A 10-arcsec inventory radius and a 3-arcsec association cut apply when target coordinates are available. Name-indexed NED products can duplicate survey spectra. A DESI stack spans the exposure MJD range recorded in the timeline.\n\n')
        f.write('The fit uses a local linear continuum and a common narrow-line velocity and width, with fixed [O III] and [N II] doublet ratios. The alternate model adds one Gaussian Balmer component with FWHM between 1,201 and 14,129 km/s. “Broad favored” requires ΔBIC ≥ 10, integrated-flux signal-to-noise ≥ 5 and an interior broad solution. The formal errors use pixel weights and residual scatter. Stellar absorption, Fe II, aperture differences and flux calibration affect comparisons. An inconclusive fit gives no disappearance limit. Continuum windows are 50 Å wide in the rest frame and include host light. Integrated line fluxes in the CSV use the native flux-density scale times observed Å. Rows on the same night or repeated between NED and a survey can represent the same exposure.\n')
        f.write('\nAn independent broad-line transition requires distinct epochs with comparable line sensitivity. Catalog spectral types and multiple reductions of one epoch cannot establish one. Searches are limited to responding public holdings.\n')
