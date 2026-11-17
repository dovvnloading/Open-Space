"""Restore an interpolation grid from a published spacecraft-observables CSV."""
import argparse
import csv
from pathlib import Path
import numpy as np
from geometry import MOONS, PLANET


def restore(table, output, telescope):
    rows=list(csv.DictReader(table.open()))
    bodies=list(dict.fromkeys(int(r['body']) for r in rows))
    planet=next(p for p,b in PLANET.items() if b==bodies[0])
    if bodies!=[PLANET[planet],*MOONS[planet]]:raise ValueError('Unexpected catalog')
    grids=[];reference=None
    for b in bodies:
        block=[r for r in rows if int(r['body'])==b]
        dates=np.array([float(r['mjd_utc']) for r in block])
        if reference is None:reference=dates
        if not np.array_equal(dates,reference):raise ValueError('Unequal grids')
        grids.append([[float(r[k]) if r[k] else np.nan for k in ('ra_icrf_deg','dec_icrf_deg','diam_arcsec','r_au','delta_au')] for r in block])
    cuts=np.flatnonzero(np.diff(dates)>10/1440)+1
    lengths=np.diff(np.r_[0,cuts,len(dates)])
    np.savez_compressed(output,mjd=dates,ids=bodies,values=grids,planet=planet,observer=telescope,window_lengths=lengths)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('table',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--telescope',choices=['HST','JWST'],default='JWST')
    a=p.parse_args();restore(a.table,a.output,a.telescope)
