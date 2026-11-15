"""Test the published OGLE point-lens shape against independent VIRAC2 Ks photometry."""
import argparse
import csv
import json
from pathlib import Path
import urllib.parse
import urllib.request

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
TAP = 'https://archive.eso.org/tap_cat/sync'
CANDIDATES = [
    # DN suffix, OGLE table-3 coordinates and best-fit t0, tE, u0.
    ('0022', 261.96975, -29.81086, 2455719.461, 6.253, 0.0),
    ('0067', 265.00708, -34.21981, 2456560.35, 7.203, 0.824),
    ('0083', 265.70112, -20.75825, 2456341.534, 38.291, 0.413),
    ('0135', 267.05867, -24.36825, 2456545.566, 99.097, 0.172),
    ('0481', 270.38233, -30.44986, 2456155.77, 14.314, 0.602),
    ('0513', 270.57625, -31.92611, 2455641.772, 3.144, 0.744),
    ('0787', 272.30358, -32.10417, 2456419.605, 7.991, 0.748),
]
def tap(query):
    params = urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':query})
    with urllib.request.urlopen(TAP+'?'+params, timeout=90) as response:
        return json.load(response)['data']

def fetch(suffix,ra,dec):
    q=(f'SELECT TOP 30 sourceid,ra,de,phot_ks_mean_mag,phot_ks_std_mag,phot_ks_n_epochs '
       f'FROM VVVX_VIRAC_V2_SOURCES WHERE ra BETWEEN {ra-.001} AND {ra+.001} '
       f'AND de BETWEEN {dec-.001} AND {dec+.001}')
    sources=tap(q)
    def distance(row):
        return 3600*np.hypot((row[1]-ra)*np.cos(np.radians(dec)),row[2]-dec)
    nearest=min(sources,key=distance)
    source=int(nearest[0])
    target_rows=[[int(r[0]),r[1],r[2],r[3],r[4],r[5]] for r in sources]
    (HERE/'vista'/f'{suffix}-sources.json').write_text(json.dumps(target_rows,indent=1)+'\n')
    if distance(nearest)>2: return
    phot=tap(f'SELECT mjdobs,filter,mag,emag,phot_flag,ambiguous_match FROM VVVX_VIRAC_V2_LC WHERE sourceid={source}')
    with (HERE/'vista'/f'{suffix}-phot.csv').open('w',newline='') as out:
        w=csv.writer(out);w.writerow(['mjdobs','filter','mag','emag','phot_flag','ambiguous_match'])
        w.writerows(sorted(phot,key=lambda r:r[0]))

def fit(suffix,ra,dec,t0,te,u):
    sources=json.loads((HERE/'vista'/f'{suffix}-sources.json').read_text())
    sep,source_id=min((3600*np.hypot((r[1]-ra)*np.cos(np.radians(dec)),r[2]-dec),int(r[0])) for r in sources)
    result=dict(dn=suffix,closest_source_id=source_id,offset_arcsec=float(sep))
    if sep>2:return result
    with (HERE/'vista'/f'{suffix}-phot.csv').open() as f:
        rows=list(csv.DictReader(f))
    obs=[(float(r['mjdobs'])+2400000.5,float(r['mag']),float(r['emag'])) for r in rows
         if r['filter']=='Ks' and r['phot_flag']=='0' and r['ambiguous_match']=='0' and r['mag'] and float(r['emag'])>0]
    t,m,e=np.array(obs).T
    _,index=np.unique(np.floor(t),return_inverse=True)
    weight=np.bincount(index,weights=e**-2)
    t=np.bincount(index,weights=t/e**2)/weight
    m=np.bincount(index,weights=m/e**2)/weight
    e=weight**-.5
    off_peak=abs(t-t0)>3*te
    baseline=np.median(m[off_peak]);scatter=1.4826*np.median(abs(m[off_peak]-baseline))
    flux=10**(-.4*(m-baseline))
    err=np.hypot(e,scatter)*(np.log(10)/2.5)*flux
    in_peak=abs(t-t0)<te
    result.update(ks_nights=len(t),within_te=int(sum(in_peak)),baseline_ks=float(baseline),off_peak_scatter=float(scatter))
    if not any(in_peak):return result
    def evaluate(center):
        h=((t-center)/te)**2+u*u
        A=(h+2)/np.sqrt(h*(h+4))
        X=np.stack((np.ones(len(t)),A-1),axis=1)
        B=X/err[:,None]
        coeff=np.linalg.lstsq(B,flux/err,rcond=None)[0]
        flat=np.average(flux,weights=err**-2)
        chi_improvement=np.sum(((flux-flat)/err)**2)-np.sum(((flux-X@coeff)/err)**2)
        return float(coeff[1]),float(np.linalg.inv(B.T@B)[1,1]**.5),float(chi_improvement)
    b,uncertainty,delta=evaluate(t0)
    controls=[]
    for center in np.linspace(t.min(),t.max(),1000):
        if abs(center-t0)<=3*te or sum(abs(t-center)<te)<2:continue
        amp,_,score=evaluate(center)
        controls.append(score if amp>0 else -score)
    result.update(brightest_sample_ks=float(min(m[in_peak])),relative_source_flux=b,formal_fit_error=uncertainty,delta_chi2=delta,
                  eligible_time_shifts=len(controls),stronger_time_shifts=sum(v>=delta for v in controls))
    # A flexible single intrinsic outburst also brightens in both bands. Test a
    # two-sided exponential as a competing description of the infrared time series.
    def flare_residual(par):
        center,rise,fall,constant,amplitude=par
        dt=t-t0-center
        profile=np.exp(-abs(dt)/np.where(dt<0,rise,fall))
        return (flux-constant-amplitude*profile)/err
    fits=[]
    for shift in (-10,0,10):
        for scale in (2,20,80):
            initial=[shift,scale,scale,1.,.8]
            bounds=([-te,.1,.005,.05,.0001],[te,2000,2000,5.,15.])
            fits.append(least_squares(flare_residual,initial,bounds=bounds,max_nfev=3000))
    alternate=min(fits,key=lambda v:v.cost)
    flat=np.average(flux,weights=err**-2)
    fixed_chi2=float(np.sum(((flux-flat)/err)**2)-delta)
    result.update(fixed_shape_chi2=fixed_chi2,free_flare_chi2=float(2*alternate.cost),
                  free_flare_peak_relative_days=float(alternate.x[0]),
                  free_flare_rise_days=float(alternate.x[1]),free_flare_decay_days=float(alternate.x[2]))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh',action='store_true',help='Redownload public VISTA rows prior to fitting')
    args=parser.parse_args()
    if args.refresh:
        for suffix,ra,dec,*_ in CANDIDATES:fetch(suffix,ra,dec)
    print(json.dumps([fit(*r) for r in CANDIDATES],indent=2))
