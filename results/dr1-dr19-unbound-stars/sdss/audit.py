"""Reproduce controls, bright-neighbor search, independent RV checks and spectrum audits."""
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
import io
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
from fetch_survey import HERE,WORK,TABLES,query
from fetch_tail import GAIA
from rank import read_raw,rank,selected
A=HERE/'audit';A.mkdir(exist_ok=True)
IDS=['57494348','62358985','70819612','88827329','78976052','80699036','88803696','61251755','88392741','59531645','66849877','70337258','73553012','80845910','61938705','61938721','61938747','88734819','73879498','101310339']
SPECTRA=['59531645','88392741','101310339','73879498','73553012','61938705','80845910','88734819','70337258','66849877','88827329']
BOSS_LINES={'Hbeta':4862.68,'Halpha':6564.61,'CaII_1':8500.36,'CaII_2':8544.44,'CaII_3':8664.52}
def csv_read(content):return list(csv.DictReader(content.decode().splitlines()))
def write_rows(p,rows):
    if not rows: return
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def tap(service,q):
    req=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q)
    with urlopen(service+'?'+urlencode(req),timeout=120) as response:return response.read()
def spectrum(k):
    rawdir=WORK/'spectra';rawdir.mkdir(exist_ok=True)
    url=f'https://data.sdss.org/sas/dr19/spectro/astra/0.6.0/spectra/visit/{k[-4:-2]}/{k[-2:]}/mwmVisit-0.6.0-{k}.fits'
    p=rawdir/f'{k}.fits'
    if not p.exists():
        with urlopen(url,timeout=120) as r:p.write_bytes(r.read())
    hsh=sha256(p.read_bytes()).hexdigest(); minima=[]
    with fits.open(p) as f:
        h=f['BOSS/APO']
        if len(h.data):
            row=h.data[0];y=row['flux'].astype(float);iv=row['ivar']
            good=np.isfinite(y)&(iv>0)
            y[~good]=np.nanmedian(y[good])
            w=10**(h.header['CRVAL']+np.arange(len(y))*h.header['CDELT'])
            med=median_filter(y,91)
            yn=np.divide(y,med,out=np.full_like(y,np.nan),where=med>0)
            for name,c in BOSS_LINES.items():
                dv=(w/c-1)*299792.458
                use=good&(np.abs(dv)<900)&np.isfinite(yn)
                if not any(use):continue
                i=np.argmin(np.where(use,yn,np.inf))
                minima.append(dict(sdss_id=k,mjd=row['mjd'],line=name,xcsao_kms=row['xcsao_v_rad'],
                                   vacuum_min_angstrom=round(w[i],3),rest_residual_kms=round(dv[i],1),
                                   normalized_flux=round(yn[i],4)))
    return dict(sdss_id=k,url=url,bytes=p.stat().st_size,sha256=hsh),minima
def neighbor(a):
    ra,dec=a['ra'],a['dec']
    q=f"SELECT source_id,ra,dec,phot_g_mean_mag,parallax,parallax_error,pmra,pmdec,ruwe,DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra},{dec}))*3600 AS sep_arcsec FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.001666667)) ORDER BY source_id"
    return [dict(sdss_id=a['sdss_id'],**r) for r in csv_read(tap(GAIA,q))]
def simbad(a):
    ra,dec=a['ra'],a['dec']
    q=f"SELECT TOP 10 b.main_id,b.otype,b.ra,b.dec,b.rvz_radvel,b.rvz_err,b.rvz_bibcode,b.sp_type,b.plx_value FROM basic b WHERE 1=CONTAINS(POINT('ICRS',{ra},{dec}),CIRCLE('ICRS',b.ra,b.dec,0.001))"
    return [dict(sdss_id=a['sdss_id'],**r) for r in csv_read(tap('https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',q))]
def distance_model(row):
    rng=np.random.default_rng(876793);n=100000
    mass=rng.uniform(.6,1.6,n)
    t=rng.normal(float(row['teff']),np.hypot(float(row['e_teff']),150),n)
    grav=rng.normal(float(row['logg']),np.hypot(float(row['e_logg']),.30),n)
    u=t-5772
    bc=np.polynomial.polynomial.polyval(u,[.06,6.731e-5,-6.647e-8,2.859e-11,-7.197e-15])+rng.normal(0,.15,n)
    lum=mass*10**(4.438-grav)*(t/5772)**4
    ag=np.maximum(0,rng.normal(2.74*float(row['ebv']),.10,n))
    return 10**((float(row['g_mag'])-ag-(4.74-2.5*np.log10(lum)-bc)-10)/5)
def controls():
    out=[]
    for t in TABLES:
        for lbl in ['prefilter','control']:
            rows=list(read_raw(t,lbl))
            for subset,f in [('all',lambda a:True),('rv_ge650',lambda a:abs(float(a['v_rad']))>=650),('tail',lambda a:selected(rank(a)))]:
                s=[a for a in rows if f(a)]
                good=[a for a in s if 0<float(a['gaia_e_v_rad'])<=15]
                out.append(dict(table=t,selection=lbl,subset=subset,rows=len(s),gaia_rv_error_le15=len(good),
                             discordance_gt100=sum(abs(float(a['v_rad'])-float(a['gaia_v_rad']))>100 for a in good)))
    write_rows(A/'rv_controls.csv',out)
if __name__=='__main__':
    controls()
    sel=("SELECT sdss_id,gaia_dr3_source_id,teff,e_teff,logg,e_logg,m_h_atm,e_m_h_atm,j_mag,h_mag,k_mag,ebv,g_mag,spectral_type,classification,ra,dec,plx,pmra,pmde "
         "FROM lite_all_star WHERE sdss_id IN ("+','.join(IDS)+") ORDER BY sdss_id")
    body=query(sel).split(b'\n',1)[1];(A/'astra_parameters.csv').write_bytes(body);rs=csv_read(body)
    with ThreadPoolExecutor(5) as pool:
        write_rows(A/'gaia_neighbors.csv',[r for result in pool.map(neighbor,rs) for r in result])
        write_rows(A/'simbad.csv',[r for result in pool.map(simbad,[a for a in rs if a['sdss_id'] in {'66849877','70337258','59531645','78976052','88803696'}]) for r in result])
        res=list(pool.map(spectrum,SPECTRA))
    write_rows(A/'spectra_manifest.csv',[r[0] for r in res])
    write_rows(A/'line_minima.csv',[a for r in res for a in r[1]])
    q="SELECT source_id,best_classification,p1_o,int_average_g,int_average_g_error,peak_to_peak_g,metallicity,average_rv FROM gaiadr3.vari_rrlyrae WHERE source_id=1170190610327078144"
    (A/'rrlyrae.csv').write_bytes(tap(GAIA,q))
    q="SELECT specobjid,ra,dec,class,subclass,z,zerr,zwarning,plate,mjd,fiberid FROM specObjAll WHERE ra BETWEEN 8.270 AND 8.276 AND dec BETWEEN -0.262 AND -0.256 ORDER BY mjd"
    (A/'stripe82_spectra.csv').write_bytes(query(q).split(b'\n',1)[1])
    dist=distance_model(next(r for r in rs if r['sdss_id']=='88827329'))
    (A/'gravity_distance.json').write_text(json.dumps({'sdss_id':'88827329','percentiles':[2.5,16,50,84,97.5],
          'distance_kpc':np.round(np.percentile(dist,[2.5,16,50,84,97.5]),5).tolist(),'fraction_at_least_6kpc':float((dist>=6).mean())},indent=2)+'\n')
    manifests=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha256(p.read_bytes()).hexdigest()) for p in sorted(WORK.glob('*.csv'))]
    write_rows(A/'survey_manifest.csv',manifests)
