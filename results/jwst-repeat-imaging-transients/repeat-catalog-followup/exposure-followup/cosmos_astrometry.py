"""Estimate the COSMOS-009 displacement against its immediate neighboring sources."""
import json
from io import StringIO
from urllib.request import urlopen
from urllib.parse import urlencode
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "cache",
    type=Path,
    help="Directory of catalog NPZ files from the parent search")
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir(exist_ok=True, parents=True)
p = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(p))
import campaign  # Load the shared catalog implementation.
from compare import read, mutual, project

origin = np.array([150.1458367, 2.283631908])
cache = args.cache
rows = []
for row in pd.read_csv(p / "pairs.csv").itertuples():
    if row.pair != "cosmos_medium":
        continue
    a, b = [read(cache / (k + ".npz")) for k in [row.reference, row.search]]
    xa, xb = [project(d, origin) for d in (a, b)]
    i, j = mutual(xa, xb, .35)
    delta = xb[j] - xa[i]
    offset = np.median(delta, axis=0)
    i, j = mutual(xa + offset, xb, .09)
    delta = xb[j] - xa[i]
    sep = np.hypot(*xb[j].T)
    s = sep.argmin()
    quality = (
        a.flux.values[i] /
        a.flux_err.values[i] > 20) & (
        b.flux.values[j] /
        b.flux_err.values[j] > 20)
    ids = np.flatnonzero(quality & (sep > 1) & (sep < 30))
    ids = ids[np.argsort(sep[ids])[:40]]
    med = np.median(delta[ids], axis=0)
    sig = 1.4826 * np.median(abs(delta[ids] - med), axis=0)
    d = (delta[s] - med) * 1000
    rows.append(dict(band=row.band, controls=len(ids), max_control_radius_arcsec=sep[ids].max(),
                     motion_e_mas=d[0], motion_n_mas=d[1], control_mad_e_mas=sig[0] * 1000, control_mad_n_mas=sig[1] * 1000,
                     controls_with_larger_2d_displacement=int((np.hypot(*(delta[ids] - med).T) > np.hypot(*(delta[s] - med))).sum())))
out = pd.DataFrame(rows)
out.to_csv(args.output / "local_motion.csv", index=False)

center = [150.1458367, 2.283631908]
query_args = {
    '-source': 'J/ApJS/214/24/3dhstall',
    '-c': '150.145837 +2.283632',
    '-c.rs': 55,
    '-out': 'ID,RAJ2000,DEJ2000,S/G,F160W,e_F160W',
    '-out.max': 900}
url = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?' + urlencode(query_args)
f = Path(__file__).resolve().parent / 'cosmos_medium-009' / 'nearby_hst.tsv'
if not f.exists():
    f.write_bytes(urlopen(url, timeout=45).read())
ls = [r for r in f.read_text().splitlines() if r.strip()
      and not r.startswith('#')]
a = pd.read_csv(StringIO('\n'.join(
    [ls[0]] + ls[3:])), sep='\t').rename(columns={'RAJ2000': 'ra', 'DEJ2000': 'dec'})
a["F160W"] = pd.to_numeric(a["F160W"], errors="coerce")
a["e_F160W"] = pd.to_numeric(a["e_F160W"], errors="coerce")
pa = project(a, center)
rows = []
for p in args.cache.glob('jw07814-o003*'):
    b = read(p)
    pb = project(b, center)
    ia, ib = mutual(pa, pb, .9)
    delta = pb[ib] - pa[ia]
    sep = np.linalg.norm(pb[ib], axis=1)
    quality = (
        sep > 1) & (
        sep < 40) & (
            b.flux.values[ib] >= 20 *
            b.flux_err.values[ib]) & (
                a['S/G'].values[ia] != 1) & (
                    a.F160W.values[ia] > 10 *
        a.e_F160W.values[ia])
    idx = np.flatnonzero(quality)
    idx = idx[np.argsort(sep[idx])[:40]]
    med = np.median(delta[idx], axis=0)
    target = np.where(a.ID.values[ia] == 10232)[0]
    j = target[0]
    diff = delta[j] - med
    rows.append(dict(product=p.stem, controls=len(idx), hst_id=10232, match_sep_arcsec=np.linalg.norm(delta[j]), east_offset_arcsec=diff[0], north_offset_arcsec=diff[1], radial_offset_arcsec=np.hypot(*diff),
                     control_mad_arcsec=(1.4826 * np.median(abs(delta[idx] - med), axis=0)).tolist(), controls_with_greater_shift=int(sum(np.linalg.norm(delta[idx] - med, axis=1) > np.hypot(*diff)))))
(args.output / 'hst_motion.json').write_text(json.dumps(dict(url=url, rows=rows), indent=2) + '\n')
