"""Rebuild the optical control using three public catalogs, or rerun from member candidates.

python optical.py                       # local, uses included compact optical-matches.csv
python optical.py --fetch /tmp/ngc6517  # fetch the release queries, rebuild table and JSON
"""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
RA0 = 270.4605
DEC0 = -(8 + 57/60 + 31.6/3600)
PM_CENTRE = (-1.554, -4.452)  # Vasiliev & Baumgardt, mas/yr

def mast(endpoint, params):
    return 'https://catalogs.mast.stsci.edu/api/v0.1/' + endpoint + '?' + urlencode(params)

def queries():
    hsc = dict(ra='270.4605', dec='-8.9588', radius='0.02777778', pagesize='100000',
               columns='[MatchID,MatchRA,MatchDec,W3_F555W,W3_F814W,CI,CI_Sigma,DSigma,W3_F555W_N,W3_F814W_N]',
               **{'W3_F814W.gte':'16', 'W3_F814W.lte':'23', 'W3_F555W.gte':'16', 'W3_F555W.lte':'28'})
    ps = dict(ra='270.4605', dec='-8.9588', radius='0.03', pagesize='100000',
              columns='[objID,raMean,decMean,nDetections,nr,ni,nz,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,rMeanPSFMagErr,iMeanPSFMagErr,zMeanPSFMagErr]',
              **{'nDetections.gte':'3'})
    adql = ("SELECT source_id,ra,dec,pmra,pmdec,pmra_error,pmdec_error,parallax,parallax_error,"
            "phot_g_mean_mag,bp_rp,ruwe,phot_bp_rp_excess_factor FROM gaiadr3.gaia_source WHERE "
            "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',270.4605,-8.9588,0.03))")
    return {'cluster': mast('hsc/v3/summary/magaper2.csv', hsc),
            'ps1': mast('panstarrs/dr2/mean.csv', ps),
            'gaia': 'https://gaia.ari.uni-heidelberg.de/tap/sync?' +
                    urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=adql))}

def xy(ra, dec):
    return np.column_stack([(ra-RA0)*3600*np.cos(np.radians(DEC0)), (dec-DEC0)*3600])

def load_catalogs(cache):
    cache.mkdir(exist_ok=True, parents=True)
    tables, receipts = {}, {}
    for name, url in queries().items():
        path = cache / (name + '.csv')
        if not path.exists():
            with urlopen(url, timeout=150) as response:
                path.write_bytes(response.read())
        df = pd.read_csv(path)
        tables[name] = df
        receipts[name] = dict(url=url, rows=len(df), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        if len(df) >= 100000:
            raise ValueError('Query truncated: paginate before analysis')
    return tables, receipts

def match_members(tables):
    g, h, p = (tables[k] for k in ['gaia','cluster','ps1'])
    g = g.dropna(subset=['pmra','pmdec']).copy()
    g['dpm'] = np.hypot(g.pmra-PM_CENTRE[0], g.pmdec-PM_CENTRE[1])
    # Match against *all* Gaia sources, then apply member/quality cuts to the match.
    ref = tables['gaia'].set_index('source_id')
    tree = cKDTree(xy(ref.ra, ref.dec))
    g = g[(g.dpm<2) & (g.ruwe<1.4) & (g.pmra_error<.6) & (g.pmdec_error<.6)].set_index('source_id')
    catalogs = [
        ('HST', h.MatchID, h.MatchRA, h.MatchDec, h.W3_F814W, h.W3_F555W-h.W3_F814W,
         (h.CI<1.25) & (h.CI_Sigma<.1) & h.W3_F814W.between(16,20.5)),
        ('PS1', p.objID, p.raMean, p.decMean, p.iMeanPSFMag, p.rMeanPSFMag-p.iMeanPSFMag,
         (p.nr>=2) & (p.ni>=2) & (p.rMeanPSFMagErr<.1) & (p.iMeanPSFMagErr<.1)),
        ('Gaia', g.index.to_series(), g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp,
         g.bp_rp.between(1.4,3.5))]
    matches = []
    for name, ids, ra, dec, mag, color, good in catalogs:
        pos = xy(ra, dec)
        distance, idx = tree.query(pos)
        rows = pd.DataFrame(dict(catalog=name, object_id=np.asarray(ids), gaia_id=ref.index.values[idx],
                                 east=pos[:,0], north=pos[:,1], separation_as=distance,
                                 mag=np.asarray(mag), color=np.asarray(color)))
        rows = rows[np.asarray(good & mag.between(14,20.5) & color.between(.35,3.5)) &
                    (distance<.35) & rows.gaia_id.isin(g.index)]
        rows = rows.sort_values(['separation_as','object_id']).drop_duplicates('gaia_id')
        rows['dpm'] = rows.gaia_id.map(g.dpm)
        rows['radius'] = np.hypot(rows.east, rows.north)
        matches.append(rows[rows.radius.between(10,85,inclusive='neither')])
    return pd.concat(matches).sort_values(['catalog','gaia_id'])

def fit(rows):
    """Robust location of color vs. sky, magnitude/ridge curvature and crowding radius."""
    mag = rows.mag.to_numpy() - np.median(rows.mag)
    design = np.column_stack([np.ones(len(rows)), rows.east/60, rows.north/60,
                              mag, mag**2, rows.radius/60])
    return least_squares(lambda b: design@b-rows.color.to_numpy(), np.zeros(6),
                         loss='soft_l1', f_scale=.06).x

def analyze(rows):
    result = {}
    for name, stars in rows.groupby('catalog'):
        sample = stars[stars.dpm<1]
        coefficient = fit(sample)
        # Eight delete-one angular sectors; this checks sensitivity to spatial coverage.
        sector = np.floor((np.degrees(np.arctan2(sample.east,sample.north))%360)/45).astype(int)
        east = [float(fit(sample[sector!=i])[1]) for i in range(8)]
        result[name] = dict(n=len(sample), east_mag_per_arcmin=coefficient[1],
                            north_mag_per_arcmin=coefficient[2], sector_counts=np.bincount(sector,minlength=8).tolist(),
                            delete_sector_east_range=[min(east),max(east)],
                            alternative_east={})
        cuts = {'strict_pm_0.5':stars.dpm<.5, 'loose_pm_2':stars.dpm<2,
                'exclude_inner_30':(stars.dpm<1)&(stars.radius>30),
                'within_60':(stars.dpm<1)&(stars.radius<60),
                'match_0.15':(stars.dpm<1)&(stars.separation_as<.15)}
        for label, keep in cuts.items():
            sub = stars[keep]
            result[name]['alternative_east'][label] = {'n': len(sub), 'slope': fit(sub)[1]}
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', type=Path, metavar='CACHE')
    args = parser.parse_args()
    if args.fetch:
        cats, receipts = load_catalogs(args.fetch)
        match_members(cats).to_csv(HERE/'optical-matches.csv',index=False,float_format='%.9f')
        (HERE/'optical-inputs.json').write_text(json.dumps(receipts,indent=2)+'\n')
    output = analyze(pd.read_csv(HERE/'optical-matches.csv'))
    (HERE/'optical-analysis.json').write_text(json.dumps(output,indent=2)+'\n')
