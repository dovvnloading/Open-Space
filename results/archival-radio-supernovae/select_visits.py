"""Rebuild the initial selection from downloaded ZTF BTS and VLASS inventories."""
import argparse
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
BTS = 'https://sites.astro.caltech.edu/ztf/bts/explorer.php?f=s&format=csv&subsample=sn&endz=0.05'
TAP = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync?'

def sex(v):
    a, b, c = map(float, v.lstrip('+-').split(':'))
    return (a + b / 60 + c / 3600) * (-1 if v.startswith('-') else 1)

def xyz(ra, dec):
    r, d = (np.radians(ra), np.radians(dec))
    return np.array([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]).T

def jobs(cache):
    sn = pd.read_csv(cache / 'bts.csv')
    t = sn.type
    near = t.str.startswith(('SN II', 'SN Ib', 'SN Ic')) & sn.redshift.le(0.015)
    special = t.isin(['SN IIn', 'SN Ibn', 'SN Icn', 'SN Ia-CSM', 'SN Ic-BL']) & sn.redshift.le(0.05)
    sn = sn[(near | special) & sn.redshift.gt(0)].copy()
    sn['ra'] = sn.RA.apply(sex) * 15
    sn['dec'] = sn.Dec.apply(sex)
    sn = sn.sort_values(['redshift', 'ZTFID'])
    out = []
    for e in [1, 2, 3, 4]:
        names = pd.read_csv(cache / f't{e}.csv').uri
        parsed = names.str.extract('\\.J(\\d\\d)(\\d\\d)(\\d\\d)([+-])(\\d\\d)(\\d\\d)(\\d\\d)\\.')
        ra = (parsed[0].astype(float) + parsed[1].astype(float) / 60 + parsed[2].astype(float) / 3600) * 15
        dec = (parsed[4].astype(float) + parsed[5].astype(float) / 60 + parsed[6].astype(float) / 3600) * np.where(parsed[3] == '+', 1, -1)
        dist, idx = cKDTree(xyz(ra, dec)).query(xyz(sn.ra, sn.dec), k=5)
        for (_, s), ids in zip(sn.iterrows(), idx):
            pool = [i for i in ids if abs(s.dec - dec[i]) < 0.57 and abs((s.ra - ra[i] + 180) % 360 - 180) * np.cos(np.radians(s.dec)) < 0.57]
            if pool:
                out.append((s, e, names[pool[0]]))
    return sorted(out, key=lambda j: (j[0].redshift, j[0].ZTFID, j[1]))

def inventories(cache):
    cache.mkdir(parents=True, exist_ok=True)
    urls = [('bts.csv', BTS)]
    for epoch in range(1, 5):
        query = f"SELECT a.uri FROM caom2.Observation o JOIN caom2.Plane p ON o.obsID=p.obsID JOIN caom2.Artifact a ON p.planeID=a.planeID WHERE o.collection='VLASS' AND o.observationID LIKE 'VLASS{epoch}%' AND p.productID LIKE '%.quicklook' AND a.uri LIKE '%tt0.subim.fits'"
        url = TAP + urlencode(dict(QUERY=query, LANG='ADQL', REQUEST='doQuery', FORMAT='csv', MAXREC=200000))
        urls.append((f't{epoch}.csv', url))
    for name, url in urls:
        path = cache / name
        if not path.exists():
            path.write_bytes(urlopen(url, timeout=180).read())
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    args = parser.parse_args()
    inventories(args.cache)
    rows = [dict(ztf=s.ZTFID, name=s.IAUID, epoch=e, uri=u, ra=s.ra, dec=s.dec, z=s.redshift, kind=s.type, peakmjd=s.peakt + 57999.5) for s, e, u in jobs(args.cache)]
    pd.DataFrame(rows).to_csv(args.cache / 'initial-selection.csv', index=False)
