"""Summarize released JWST filesets and find candidate host pointings in the GCNS."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

TAP = 'https://dc.g-vo.org/tap/sync'
ADQL = 'SELECT source_id, ra, dec, pmra, pmdec, dist_50 FROM gcns.main WHERE dist_50 < 0.1'
PARAMS = {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','MAXREC':'350000','QUERY':ADQL}
MAS_TO_RAD = np.pi / (180 * 3600 * 1000)
GROUP = ['program','observtn','visit','targprop','targname','exp_type','opticalElements','is_psf','bkgdtarg']

def frame(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    c, s = np.cos(ra), np.sin(ra)
    u, v = np.cos(dec), np.sin(dec)
    return (np.column_stack([c*u, s*u, v]), np.column_stack([-s,c,np.zeros_like(ra)]),
            np.column_stack([-c*v,-s*v,u]))

def gcns_table(path, out):
    if not path.exists():
        with urlopen(TAP + '?' + urlencode(PARAMS), timeout=180) as f:
            path.write_bytes(f.read())
    tab = pd.read_csv(path, dtype={'source_id':'str'})
    if len(tab) >= int(PARAMS['MAXREC']) or tab['source_id'].duplicated().any():
        raise ValueError('Truncated catalogue or duplicate keys')
    (out / 'gcns-receipt.json').write_text(json.dumps({'endpoint':TAP,'params':PARAMS,'rows':len(tab),
        'sha256_csv':hashlib.sha256(path.read_bytes()).hexdigest(),'reference_epoch':'J2016.0',
        'catalogue':'Gaia Collaboration 2021, A&A 649 A6','dist_50_unit':'kpc'},indent=2)+'\n')
    return tab

def match_fields(visits, stars):
    stars = stars.dropna(subset=['ra','dec','pmra','pmdec']).reset_index(drop=True)
    pos, east, north = frame(stars.ra, stars.dec)
    velocity = (stars.pmra.to_numpy()[:,None]*east + stars.pmdec.to_numpy()[:,None]*north)*MAS_TO_RAD
    tree = cKDTree(pos)
    coords = frame(visits.targ_ra.to_numpy(), visits.targ_dec.to_numpy())[0]
    # The initial cone includes the catalogue's fastest measured proper motion, from J2016
    # to the furthest tested epoch, so the final tight match cannot be clipped.
    years = (pd.to_datetime(visits.date_obs_first, utc=True, format='ISO8601')-pd.Timestamp('2016-01-01T12:00:00Z')).dt.total_seconds().to_numpy() / (365.25*86400)
    widest = max(16, np.abs(years).max())*np.linalg.norm(velocity,axis=1).max() + 5e3*MAS_TO_RAD
    neighborhoods = tree.query_ball_point(coords, widest)
    matches = []
    for i, ids in enumerate(neighborhoods):
        candidates = []
        if ids:
            p = pos[ids]; v = velocity[ids]
            # MAST target coordinates can be supplied at different epochs. Record which
            # epoch matches. This is a field selection, not a motion measurement.
            for epoch, delta in [('obs',years[i]),('2016',0.),('2000',-16.)]:
                xyz = p + delta*v
                xyz /= np.linalg.norm(xyz,axis=1)[:,None]
                sep = np.linalg.norm(xyz-coords[i],axis=1) / (MAS_TO_RAD*1000)
                for j in np.flatnonzero(sep <= 5.):
                    candidates.append((float(sep[j]),ids[j],epoch))
        if candidates:
            distance, index, epoch = min(candidates)
            star = stars.iloc[index]
            matches.append((star.source_id,round(star.dist_50*1000,3),round(distance,3),epoch,len({c[1] for c in candidates})))
        else:
            matches.append(('',None,None,'',0))
    return pd.DataFrame(matches, columns=['gcns_edr3_source_id','gcns_distance_pc','gcns_separation_arcsec','match_epoch','within_5as_count'])

def summarize(out, samples):
    tables=[]
    for name in samples:
        with gzip.open(out / (name + '.csv.gz'),'rt') as f:
            t=pd.read_csv(f,dtype={'ArchiveFileID':'str'},keep_default_na=False)
        tables.append(t)
    all_rows=pd.concat(tables,ignore_index=True).drop_duplicates(subset='ArchiveFileID')
    # Exclude combined Stage 3 association filesets. Their time and optical elements
    # overlap the primary exposure filesets and would double-count integrations.
    all_rows = all_rows[all_rows.fileSetName.str.fullmatch(r'jw\d{11}_[0-9a-z]{5}_\d{5}')]
    if len(all_rows) != all_rows.fileSetName.nunique():
        raise ValueError('Multiple records per fileset; integrate only after resolving overlaps')
    all_rows=all_rows.sort_values(['date_obs','fileSetName'])
    all_rows['effexptm']=pd.to_numeric(all_rows.effexptm,errors='coerce')
    summary=(all_rows.groupby(GROUP,sort=False,as_index=False,dropna=False).agg(
        date_obs_first=('date_obs','min'),date_obs_last=('date_obs','max'),
        release_last=('publicReleaseDate','max'),release_first=('publicReleaseDate','min'),
        targ_ra=('targ_ra','first'),targ_dec=('targ_dec','first'),proposal_type=('proposal_type','first'),
        category=('targcat','first'),title=('title','first'),filesets=('fileSetName','size'),
        example_fileset=('fileSetName','first'),obslabel=('obslabel','first'),
        effective_exposure_s=('effexptm',lambda s: s.sum(min_count=1))))
    return summary

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path(__file__).parent/'archive')
    ap.add_argument('--gcns',type=Path,required=True,help='cache path for the GCNS TAP response')
    ap.add_argument('samples',nargs='+',help='one or more sample keys from fetch_inventory.py')
    args=ap.parse_args()
    summary=summarize(args.out,args.samples)
    stars=gcns_table(args.gcns,args.out)
    groups=pd.concat([summary,match_fields(summary,stars)],axis=1)
    groups=groups.sort_values(['release_last','date_obs_last','example_fileset'],ascending=[False,False,True])
    groups.to_csv(args.out/'visits.csv',index=False,lineterminator='\n')
    good=groups[(groups.within_5as_count > 0) & (~groups.proposal_type.isin(['COM','CAL'])) &
        (groups.is_psf!='t') & (groups.bkgdtarg!='t')]
    good.to_csv(args.out/'nearby-fields.csv',index=False,lineterminator='\n')
    print('all visits/filter groups',len(groups),'nearby non-PSF non-background non-COM/CAL',len(good),
          'unique GCNS possible hosts',good.gcns_edr3_source_id.nunique(),'programs',good.program.nunique())

if __name__ == '__main__':
    main()
