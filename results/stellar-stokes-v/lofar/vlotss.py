"""Match the frozen signed V-LoTSS catalogue to each sample at the radio date."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from scipy.spatial import cKDTree
from scipy.stats import chi2

URL='https://lofar-surveys.org/public/vlotss_dr2_cat_accepted.fits'
FIELDS=['sample','target','gaia_dr3_id','vlotss_name','published_name','published_gaia','date',
        'field','distance_arcsec','target_ra_deg','target_dec_deg','stokes_i_mjy','error_i_mjy',
        'stokes_v_mjy','error_v_mjy','abs_v_over_i_pct','error_abs_v_over_i_pct','rms_v_mjy_beam',
        'rms_i_mjy_beam','high_cp_abs_ge_10pct','target_position','control_pa_deg']


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--samples',nargs='+',type=Path,required=True)
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'data')
    args=ap.parse_args()
    args.out.mkdir(exist_ok=True,parents=True)
    raw=args.out/'vlotss_dr2.fits'
    if not raw.exists():raw.write_bytes(urlopen(URL,timeout=40).read())
    radio=fits.getdata(raw,ext=1)
    assert len(radio)==68
    rows,stats=[],[]
    for path in args.samples:
        targets=pd.read_csv(path,dtype={'gaia_dr3_id':str})
        cols=['ra_deg','dec_deg','pmra_masyr','pmdec_masyr','ref_epoch']
        good=targets.dropna(subset=cols)
        coord=SkyCoord(good.ra_deg.to_numpy()*u.deg,good.dec_deg.to_numpy()*u.deg,
                       pm_ra_cosdec=good.pmra_masyr.to_numpy()*u.mas/u.yr,
                       pm_dec=good.pmdec_masyr.to_numpy()*u.mas/u.yr,
                       obstime=Time(good.ref_epoch.to_numpy(),format='jyear'))
        tree=cKDTree(coord.cartesian.xyz.value.T)
        sample_rows=[]
        for source in radio:
            date=Time(source['date_obs'])
            pos=SkyCoord(source['RA_v'],source['Dec_v'],unit='deg')
            # Upper bound (all rows) on motion from target reference epoch.
            dt=np.abs(good.ref_epoch.to_numpy()-date.jyear)
            motion=np.max(dt*np.hypot(good.pmra_masyr,good.pmdec_masyr)/1000)
            cap=2*np.sin(np.deg2rad((motion+120+5.1)/3600)/2)
            idx=np.asarray(tree.query_ball_point(pos.cartesian.xyz.value,cap),dtype=int)
            if idx.size==0:continue
            at_date=coord[idx].apply_space_motion(new_obstime=date)
            for pa in [-1,0,45,90,135,180,225,270,315]:
                center=pos if pa==-1 else pos.directional_offset_by(pa*u.deg,120*u.arcsec)
                sep=center.separation(at_date).arcsec
                for n in np.flatnonzero(sep<=5):
                    t=good.iloc[idx[n]]
                    sample_rows.append(dict(zip(FIELDS,[path.stem,t.hostname,t.gaia_dr3_id,source['VLoTSS_name'],
                        source['Common_name'],source['designation_gaia'],source['date_obs'],source['field_name'],
                        sep[n],at_date.ra.deg[n],at_date.dec.deg[n],source['Total_flux_i'],source['e_Total_flux_i'],
                        source['Total_flux_v'],source['e_Total_flux_v'],source['v_i'],source['v_i_err'],
                        source['local_rms_v'],source['local_rms_i'],int(abs(source['v_i'])>=10),
                        'on' if pa==-1 else 'offset', '' if pa==-1 else pa])))
        nctrl=sum(r['target_position']=='offset' for r in sample_rows)
        stats.append(dict(sample=path.stem,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          rows=len(targets),usable_pm=len(good),epoch_or_pm_missing=len(targets)-len(good),
                          catalogue_centers=len(radio),offset_centers=8*len(radio),match_arcsec=5,
                          offset_arcsec=120,on_pairs=len(sample_rows)-nctrl,offset_pairs=nctrl,
                          mean_false_pairs_per_unshifted_search=nctrl/8,
                          poisson_upper95_false_pairs_per_search=chi2.ppf(0.95,2*(nctrl+1))/16))
        rows.extend(sample_rows)
    with (args.out/'vlotss_matches.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,FIELDS);w.writeheader();w.writerows(rows)
    report=dict(url=URL,input_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),rows=len(radio),samples=stats)
    (args.out/'vlotss_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(stats,indent=2))


if __name__=='__main__':main()
