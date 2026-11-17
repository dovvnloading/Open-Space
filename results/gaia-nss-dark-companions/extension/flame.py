"""Extend the orbital mass screen to primary stars with Gaia FLAME track masses."""
import argparse
import csv
import hashlib
import json
import math
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent / 'gaia'))
from screen import BAD_ASTROMETRY, SB, mass_functions, mass_ratio_root, records

QUERY = ("SELECT DISTINCT n.source_id, a.mass_flame, a.mass_flame_lower, a.mass_flame_upper, a.flags_flame "
         "FROM gaiadr3.nss_two_body_orbit AS n JOIN gaiadr3.astrophysical_parameters AS a "
         "ON n.source_id=a.source_id WHERE a.mass_flame IS NOT NULL")
URL = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode(
    dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC='500000',QUERY=QUERY))


def main(cache,nss,output):
    if not cache.exists():
        with urllib.request.urlopen(URL,timeout=240) as r:
            cache.write_bytes(r.read())
    with cache.open() as f:
        rows=list(csv.DictReader(f))
    assert len(rows)<500000 and len({r['source_id'] for r in rows})==len(rows)
    masses={r['source_id']:r for r in rows if r['flags_flame'] in ('00','10')}
    iso={r['source_id'] for r in records(nss.parent/'binary.csv.gz') if r['m1_ref']=='IsocLum'}
    with (HERE.parent/'gaia/mass_screen.csv').open() as f:
        reference=list(csv.DictReader(f))
    original={r['source_id'] for r in reference}
    previous_bh={r['source_id'] for r in reference if float(r['m2_screen_msun'])>3}
    hits,modelled,other=[],set(),set()
    for r in records(nss):
        sid,typ=r['source_id'],r['nss_solution_type']
        if typ.startswith(('SB2','Eclipsing')):
            other.add(sid)
        if sid in iso or sid not in masses or (typ not in SB|{'AstroSpectroSB1'} and not typ.startswith('Orbital')):
            continue
        f,_,basis=mass_functions(r)
        if not math.isfinite(f) or f<=0: continue
        modelled.add(sid)
        a=masses[sid]
        m1=float(a['mass_flame'])
        if f<=1.4**3/(1.4+m1)**2: continue
        m2=mass_ratio_root(f,m1)
        lo=float(a['mass_flame_lower'])
        if a['flags_flame']=='10':
            lo=min(lo,.7*m1)
        hits.append(dict(source_id=sid,nss_solution_type=typ,period=r['period'],eccentricity=r['eccentricity'],
            mass_basis=basis,m1=m1,m1_lower=a['mass_flame_lower'],m1_upper=a['mass_flame_upper'],
            flame_flags=a['flags_flame'],f_msun=f,m2_screen_msun=m2,m2_at_lower_primary=mass_ratio_root(f,lo),
            original_screen=int(sid in original),original_bh_stratum=int(sid in previous_bh),
            significance=r['significance'],goodness_of_fit=r['goodness_of_fit'],flags=r['flags'],
            esa_withdrawn_astrometry=int(sid in BAD_ASTROMETRY and basis=='astrometry_dark')))
    for r in hits:
        r['has_sb2_or_eclipse']=int(r['source_id'] in other)
    hits.sort(key=lambda r:(int(r['source_id']),r['nss_solution_type']))
    with output.open('w') as f:
        w=csv.DictWriter(f,hits[0].keys(),lineterminator='\n');w.writeheader();w.writerows(hits)
    high={r['source_id'] for r in hits if r['m2_screen_msun']>3}
    all_ids={r['source_id'] for r in hits}
    meta=dict(query=QUERY,url=URL,bytes=cache.stat().st_size,sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
              downloaded_unique_ids=len(rows),flags_accepted=['00','10'],newly_modelled_stars=len(modelled),
              threshold_unique=len(all_ids),new_to_1p4_screen=len(all_ids-original),above_3=len(high),
              new_to_3_stratum=len(high-previous_bh),thresholds_union=[len(all_ids|original),len(high|previous_bh)],
              rows_by_solution=dict(Counter(r['nss_solution_type'] for r in hits)))
    (output.parent/'depth.json').write_text(json.dumps(meta,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path);p.add_argument('nss',type=Path)
    p.add_argument('--output',type=Path,default=HERE/'mass_screen.csv')
    a=p.parse_args();main(a.cache,a.nss,a.output)
