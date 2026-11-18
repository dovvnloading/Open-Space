"""Inspect OGLE-2026-BLG-0670 peak coverage, baseline, competing shapes and aliases."""
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from lxml import html

from crossmatch import coords
from fetch import OUT,SHA,write_csv
sys.path.insert(0,str(OUT.parent))
from peak_models import solve
from run_grid import read_curve

NAME='OGLE-2026-BLG-0670'
PEAK=2461193.176


def trial(args):
    factor,bins,nstarts=args
    raw=Path('.ogle-cache/raw/2026/blg-0670.dat')
    job=(NAME,raw,2461193.333,8.497)
    seed=json.loads((Path('staging/screen')/(NAME+'.json')).read_text())
    saved=OUT/'ogle0670-models.json'
    if factor==1.15 and bins==.16 and saved.exists():
        result=json.loads(saved.read_text())
    else:
        result=solve(job,seed,factor,nstarts,True,30.,bins)
    if factor==1.15 and bins==.16:saved.write_text(json.dumps(result,indent=2)+'\n')
    m=result['models']
    return dict(event=NAME,baseline_flux_factor=factor if np.isfinite(factor) else 'unconstrained',bin_d=bins,
                peak_bins=result['n_peak_window'],binary_starts=nstarts,ps_chi2=m['ps']['chi2'],fs_chi2=m['fs']['chi2'],
                fs_te=round(m['fs']['te_or_rise_d'],6),rho=round(m['fs']['rho'],6),exp_chi2=m['exp']['chi2'],bin_chi2=m['bin']['chi2'],
                binary_q=round(m['bin']['q'],5),binary_s=round(m['bin']['s'],5),binary_converged=m['bin']['converged'],
                profile_short_chi2=m['profile_short']['chi2'],profile_small_chi2=m['profile_small']['chi2'])


def compare():
    raw=Path('.ogle-cache/raw/2026/blg-0670.dat')
    t,m,e=np.loadtxt(raw,usecols=(0,1,2)).T
    use=np.isfinite(t+m+e)&(e>0)
    t,m,e=t[use],m[use],e[use]
    near=abs(t-PEAK)<25
    by_night=[]
    for jd in np.unique(np.floor(t[near])):
        ix=near&(np.floor(t)==jd)
        by_night.append(dict(utc_like_jd_floor=int(jd),hjd_min=round(t[ix].min(),6),hjd_max=round(t[ix].max(),6),
                        epochs=int(ix.sum()),median_I=round(np.median(m[ix]),4),median_error=round(np.median(e[ix]),4)))
    write_csv(OUT/'ogle0670_nightly.csv',by_night)
    ra,dec=coords('17:53:53.72','-29:43:47.3')
    rows=[]
    def nearest(tag,cat,source_url,digest):
        found=[(3600*np.hypot((a-ra)*np.cos(np.deg2rad(dec)),d-dec),name) for name,a,d in cat]
        angular,name=min(found)
        rows.append(dict(event=NAME,survey=tag,nearest_id=name,nearest_arcsec=round(angular,4),queried_rows=len(found),
                         input_url=source_url,input_sha256=digest))
    path=Path('staging/kmt/2026.html')
    if path.exists():
        records=html.fromstring(path.read_bytes()).xpath('//tr')
        cols=[e.text_content().strip() for e in records[0].xpath('./td')]
        cat=[]
        for tr in records[1:]:
            r=dict(zip(cols,[e.text_content().strip() for e in tr.xpath('./td')]))
            cat.append((r['Event'],*coords(r['RA'],r['Dec'])))
        nearest('KMT 2026 public alert roster',cat,'https://kmtnet.kasi.re.kr/ulens/event/2026/',SHA(path.read_bytes()))
    path=OUT/'index.csv'
    allalerts=list(csv.DictReader(path.open()))
    for survey in ['MOA','PRIME']:
        cat=[(r['event'],*coords(r['ra'],r['dec'])) for r in allalerts if r['event'].startswith(survey+'-2026')]
        nearest(survey+' 2026',cat,'index.csv',SHA(path.read_bytes()))
    path=OUT/'gaia.csv'
    cat=[(r['event'],float(r['ra_deg']),float(r['dec_deg'])) for r in csv.DictReader(path.open())]
    nearest('Gaia main alerts 2024–26 bulge rectangle',cat,'gaia.csv',SHA(path.read_bytes()))
    write_csv(OUT/'ogle0670_aliases.csv',rows)
    quiet=abs(t-PEAK)>40
    report={'event':NAME,'ra_deg':ra,'dec_deg':dec,'epochs':len(t),'raw_sha256':SHA(raw.read_bytes()),
            'phot_url':'https://www.astrouw.edu.pl/ogle/ogle4/ews/2026/blg-0670/phot.dat','reference_hjd':PEAK,
            'previous_hjd':max(t[t<PEAK]),'following_hjd':min(t[t>PEAK]),
            'unsampled_gap_d':round(min(t[t>PEAK])-max(t[t<PEAK]),6),
            'outside_40d_epochs':int(quiet.sum()),'outside_40d_median_I':round(float(np.median(m[quiet])),5),
            'outside_40d_I_1pct_99pct':np.quantile(m[quiet],[.01,.99]).tolist()}
    (OUT/'ogle0670_coverage.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    compare()
    with ProcessPoolExecutor(3) as ex:
        results=list(ex.map(trial,[(1.15,.16,30),(1.,.16,16),(1.05,.16,16),(np.inf,.16,16),(1.15,0.,16)]))
    write_csv(OUT/'ogle0670_profiles.csv',results)
