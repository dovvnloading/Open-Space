"""Check BOSS pipeline labels and line spectra for all five high-vtan DR20 tail entries."""
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).parent
CACHE=ROOT/'.cache'
CACHE.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT.parent/'spectra'))
from fit_control import vacuum, factor, velocity_fit

TAIL=[r for r in csv.DictReader((ROOT.parent/'sdss/tail.csv').open()) if r['release']=='DR20']
IDLIST=','.join(r['source_id'] for r in TAIL)
BASE='https://api.sdss.org/valis/target'


def get(url):
    key=hashlib.sha256(url.encode()).hexdigest()
    p=CACHE/key
    if not p.exists():
        with urlopen(url,timeout=90) as response:
            d=response.read()
        p.write_bytes(d)
    d=p.read_bytes()
    return d,dict(url=url,sha256=hashlib.sha256(d).hexdigest(),bytes=len(d))


def tap(base,q):
    return get(base+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q)))


def pEW(w,f,ivar,mask,air,rv,pad=7.0):
    # Measure a local flux deficit with continuum covariance. O/R mask retained.
    x=w-vacuum(air)*factor(rv)
    valid=(mask==0)&np.isfinite(f)&(ivar>0)
    line=valid&(abs(x)<=pad)
    side=valid&(abs(x)>pad+4)&(abs(x)<pad+14)
    dx=np.gradient(w)[line]
    if dx.sum()<1.8*pad or side.sum()<8: return None,None
    design=np.c_[np.ones_like(x),x]
    a=design[side]*np.sqrt(ivar[side,None])
    pars=np.linalg.lstsq(a,f[side]*np.sqrt(ivar[side]),rcond=None)[0]
    c=design[line]@pars
    if np.min(c)<=0: return None,None
    grad=(dx*f[line]/c**2)@design[line]
    var=np.sum(dx**2/(ivar[line]*c**2)) + grad@np.linalg.inv(a.T@a)@grad
    chi2=np.sum((f[side]-design[side]@pars)**2*ivar[side])
    return float(np.sum(dx*(1.-f[line]/c))),float(np.sqrt(var*max(1.,chi2/(sum(side)-2))))


def spectra(row):
    name=row['source_id'];sdss=row['sdss_id'];receipts={}
    def get_json(key,url):
        try:
            data,receipts[key]=get(url)
            return json.loads(data)
        except Exception as err:
            receipts[key]={'url':url,'error':str(err)}
            return {}
    params=get_json('boss_net',f'{BASE}/astra/boss_net/{sdss}?release=DR20')
    pipes=get_json('pipelines',f'{BASE}/pipelines/{sdss}?release=DR20')
    boss=get_json('boss',f'{BASE}/pipe/boss/{sdss}?release=DR20')
    records=pipes.get('boss',[])
    # Each distinct date-coadd is independent. Epoch and all-time products overlap.
    dates=[r for r in records if r['coadd']=='daily']
    if not dates and records:
        dates=[next((r for r in records if r.get('specprimary')),records[0])]
    if not dates: raise RuntimeError(f'No BOSS data returned for {name}')
    output=[]
    labels=[]
    for record in dates:
        url='https://dr20.sdss.org/sas/'+record['location']
        raw,receipts[record['location']]=get(url)
        with fits.open(io.BytesIO(raw)) as hdus:
            d=hdus[1].data
            w=10.**d['loglam'].astype(float);f=d['flux'].astype(float)
            iv=d['ivar'].astype(float);mask=d['or_mask']
            sp=hdus[2].data[0]
            v=float(sp['Z'])*299792.458
            # Independent Gaussian Balmer-core fit in windows shifted by the pipeline
            # Doppler factor to allow relative fitting +/-350 km/s for a fast target.
            data=np.c_[w/factor(v),f,np.zeros_like(f),iv,mask]
            valid=(mask==0)&(iv>0)&np.isfinite(data).all(axis=1)
            lines=[a for a in [4340.472,4861.333,6562.800] if
                   np.sum(valid&(abs(data[:,0]-vacuum(a))<13))>=12]
            measured=velocity_fit([data],lines) if lines else {}
            relative=measured.get('rv')
            vr=(v+relative)/(1+v*relative/299792.458**2) if relative is not None else None
            flags='|'.join(f'{k}={sp[k]}' for k in ('ZWARNING','OBJTYPE','FIRSTCARTON') if k in hdus[2].columns.names)
            labels.append(dict(source_id=name,sdss_id=sdss,mjd=record['mjd'],
                      class_sdss=str(sp['CLASS']),subclass_sdss=str(sp['SUBCLASS']),
                      sdss_v_kms=v,sdss_e_v_kms=float(sp['Z_ERR'])*299792.458,
                      boss_net_teff=params.get('teff'),boss_net_logg=params.get('logg'),
                      boss_net_mh=params.get('fe_h'),boss_net_rv=params.get('bn_v_r'),
                      boss_net_result_flags=params.get('result_flags'),
                      exact_gaia_match=str(pipes.get('astra',{}).get('source',{}).get('gaia_dr3_source_id'))==name,
                      balmer_air_A='/'.join(map(str,lines)),balmer_rv_kms=vr,
                      balmer_rv_err=measured.get('formal_sigma'),boss_flags=flags))
            for elem,air,pad in [('Ca II K',3933.6614,5.),('H gamma',4340.472,8.),('H beta',4861.333,8.),
                                 ('H alpha',6562.800,8.),('Mg I',4702.991,5.),('O I triplet',7774.17,8.),
                                 ('Si II',6347.11,5.),('He I',5875.62,5.),('Na I D doublet',5892.94,8.)]:
                ew,err=pEW(w,f,iv,mask,air,v,pad)
                output.append(dict(source_id=name,mjd=record['mjd'],feature=elem,air_A=air,window_half_A=pad,
                                    pEW_A=ew,sigma_A=err,lineframe_kms=v))
    return labels,output,receipts


def writecsv(path,rows):
    with (ROOT/path).open('w') as fp:
        w=csv.DictWriter(fp,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)


def main():
    q=('SELECT g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.bp_rp,g.l,g.b,g.has_xp_continuous,'
       'a.teff_gspphot,a.logg_gspphot,a.ag_gspphot,a.mh_gspphot,a.classprob_dsc_combmod_star,'
       'a.classprob_dsc_combmod_whitedwarf FROM gaiadr3.gaia_source AS g LEFT OUTER JOIN '
       f'gaiadr3.astrophysical_parameters AS a ON g.source_id=a.source_id WHERE g.source_id IN ({IDLIST})')
    gr,receipt=tap('https://gea.esac.esa.int/tap-server/tap/sync',q)
    gaia=list(csv.DictReader(io.StringIO(gr.decode())))
    if len(gaia)!=len(TAIL): raise RuntimeError('Incomplete Gaia query')
    by_id={r['source_id']:r for r in TAIL}
    for r in gaia:
        sn=by_id[r['source_id']]
        par=float(sn['parallax'])-float(sn['zp_mas'])
        r['M_G_corrected_inverse_parallax_no_extinction']=float(r['phot_g_mean_mag'])+5*math.log10(par)-10
    writecsv('gaia.csv',gaia)
    receipts={'gaia':receipt}
    quoted=','.join("'Gaia DR3 "+r['source_id']+"'" for r in TAIL)
    simbad='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
    fields='b.main_id,b.ra,b.dec,b.otype,b.sp_type,b.rvz_radvel,b.nbref'
    q=f"SELECT i.id,{fields} FROM ident AS i JOIN basic AS b ON i.oidref=b.oid WHERE i.id IN ({quoted})"
    raw,receipts['simbad_exact']=tap(simbad,q)
    if not raw.startswith(b'id,'): raise ValueError('Invalid SIMBAD exact-ID result')
    (ROOT/'simbad.csv').write_bytes(raw)
    cones=[f"1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{r['ra']},{r['dec']},0.00277778))" for r in gaia]
    raw,receipts['simbad_cones']=tap(simbad,f"SELECT {fields} FROM basic AS b WHERE "+' OR '.join(cones))
    if not raw.startswith(b'main_id,'): raise ValueError('Invalid SIMBAD positional result')
    (ROOT/'simbad_cones.csv').write_bytes(raw)
    labels=[];features=[]
    with ThreadPoolExecutor(5) as pool:
        for row,(lab,fts,rs) in zip(TAIL,pool.map(spectra,TAIL)):
            labels.extend(lab);features.extend(fts);receipts[row['source_id']]=rs
    writecsv('boss.csv',labels);writecsv('features.csv',features)
    (ROOT/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__': main()
