"""Search cataloged sources for large changes with two independent-filter detections."""
import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from fetch_catalogs import prepare

HERE = Path(__file__).resolve().parent
RADIUS = 0.09  # arcsec, used for epoch matching and cross-filter confirmation


def read(path):
    with np.load(path) as archive:
        table = pd.DataFrame({k: archive[k] for k in archive.files})
    table.columns = [k.replace('sky_centroid.', '').replace('aper_total_', '') for k in table]
    return table.loc[np.isfinite(table[['ra','dec','flux','flux_err']]).all(axis=1) &
                     (table.flux > 0) & (table.flux_err > 0)].reset_index(drop=True)


def project(frame, origin):
    return (frame[['ra', 'dec']].to_numpy() - origin) * [3600*np.cos(np.deg2rad(origin[1])), 3600]


def mutual(a, b, radius):
    distance, ib = cKDTree(b).query(a)
    _, ia = cKDTree(a).query(b)
    ok = (distance <= radius) & (ia[ib] == np.arange(len(a)))
    return np.flatnonzero(ok), ib[ok]


def compare(row, cache):
    a, b = [read(cache/(getattr(row,k)+'.npz')) for k in ('reference','search')]
    origin = b[['ra','dec']].median().to_numpy()
    pa, pb = project(a,origin), project(b,origin)
    ia, ib = mutual(pa,pb,.35)
    offset = np.median(pb[ib]-pa[ia],axis=0)
    ia, ib = mutual(pa+offset,pb,RADIUS)
    a, b = a.iloc[ia].reset_index(drop=True), b.iloc[ib].reset_index(drop=True)
    # Isolated clean centroids. Even the finest NIRCam pixels subtend ~0.03".
    eligible = (a.nn_dist>12) & (b.nn_dist>12)
    controls = eligible & (a.flux/a.flux_err>15) & (b.flux/b.flux_err>15)
    if controls.sum()<25:
        raise ValueError(f'{row.pair} {row.band}: too few photometric controls')
    scale = np.median((b.flux/a.flux)[controls])
    fa, fb = a.flux.to_numpy()*scale, b.flux.to_numpy()
    diff = fb-fa
    meanflux = (fb+fa)/2
    error = np.hypot(a.flux_err*scale,b.flux_err)
    # In bins of mean two-epoch flux use the larger of the catalog error and the
    # observed scatter of matched non-crowded sources; bright resolved objects
    # also supply controls. Use the same cuts for both polarities.
    index = np.flatnonzero(eligible)
    groups = np.array_split(index[np.argsort(meanflux[index])], max(1,min(18,len(index)//60)))
    group_flux, scatter = [], []
    for ids in groups:
        group_flux.append(np.median(meanflux[ids]))
        scatter.append(1.4826*np.median(abs(diff[ids]-np.median(diff[ids]))))
    sigma = np.maximum.reduce([error, np.interp(meanflux,group_flux,scatter), .025*meanflux])
    sn = diff/sigma
    mag = -2.5*np.log10(fb/fa)
    valid = eligible & ((a.flux/a.flux_err)>=3) & ((b.flux/b.flux_err)>=3)
    selected = np.flatnonzero(valid & (abs(sn)>=5) & (abs(mag)>=.30))
    stats = dict(pair=row.pair,band=row.band,mutual=len(a),eligible=int(valid.sum()),controls=int(controls.sum()),
                 dx_arcsec=offset[0],dy_arcsec=offset[1],flux_scale=scale,triggers=len(selected))
    phot = pd.DataFrame(dict(pair=row.pair,band=row.band,ra=b.ra,dec=b.dec,
                    id_ref=a.label,id_search=b.label,flux_ref_njy=fa*1e9,flux_search_njy=fb*1e9,
                    err_ref_njy=a.flux_err*scale*1e9,err_search_njy=b.flux_err*1e9,delta_err_njy=sigma*1e9,
                    z_empirical=sn,delta_mag=mag,extended_ref=a.is_extended,extended_search=b.is_extended,
                    eligible=valid))
    return stats,phot,selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    args=parser.parse_args()
    prepare(args.cache)  # Verify frozen raw bytes before any comparison output is written.
    summary, triggers, chosen, curves = [],[],[],[]
    for name, pairs in pd.read_csv(HERE/'pairs.csv').groupby('pair'):
        bands={}
        for row in pairs.itertuples():
            stat, frame, subset = compare(row,args.cache)
            summary.append(stat)
            bands[row.band]=frame
            triggers.extend(frame.iloc[subset].to_dict('records'))
        subset=pd.DataFrame([r for r in triggers if r['pair']==name])
        if len(subset)<2:
            continue
        origin=subset[['ra','dec']].mean().to_numpy()
        xy=project(subset,origin)
        used=set()
        for k in np.argsort(-abs(subset.z_empirical.to_numpy())):
            if k in used:
                continue
            ix=cKDTree(xy).query_ball_point(xy[k],RADIUS)
            # Seed with the strongest remaining position; require a distinct band
            # and a common sign, and only assign a trigger once.
            ix=[i for i in ix if i not in used and subset.z_empirical[i]*subset.z_empirical[k]>0]
            group=subset.iloc[ix].sort_values('band')
            if group.band.nunique()<2:
                continue
            used.update(ix)
            center=group[['ra','dec']].median().to_numpy()
            cid=f'{name}-{sum(r["pair"]==name for r in chosen)+1:03d}'
            chosen.append(dict(candidate=cid,pair=name,ra=center[0],dec=center[1],bands=','.join(group.band),
                               peak_sigma=abs(subset.z_empirical[k]),sign=int(np.sign(subset.z_empirical[k]))))
            for band,table in bands.items():
                ds=np.linalg.norm(project(table,center),axis=1)
                j=ds.argmin()
                if ds[j]<=RADIUS:
                    curves.append(dict(candidate=cid,**table.iloc[j].to_dict()))
    for filename,rows in [('screen_counts',summary),('single_band',triggers),('coincidences',chosen),('matched_photometry',curves)]:
        pd.DataFrame(rows).to_csv(HERE/(filename+'.csv'),index=False,float_format='%.10g')
    print(len(summary),'pairs by filter,',len(triggers),'single-filter triggers,',len(chosen),'coincidences')


if __name__=='__main__':
    main()
