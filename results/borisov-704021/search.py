"""Replay the coarse Borisov orbit survey with two JPL element catalogues and a comet vector table.

python search.py ELEMENTS.NUMBR.gz ELEMENTS.UNNUM.gz > shortlist.csv
Optional input vector table is --vectors (the default is the released full-arc positions).
A screen of Kepler orbits is intentionally permissive and is not a completeness test.
"""
import argparse
import csv
import gzip
from pathlib import Path
import sys
import numpy as np


def chunks(path, width, number):
    with gzip.open(path, 'rt', encoding='ascii', errors='replace') as f:
        next(f)
        next(f)
        rows = []
        for line in f:
            data = line[width:].split()
            if len(data) < 7:
                continue
            orbital = list(map(float, data[:7]))
            if orbital[1] <= 0 or orbital[2] >= 1:
                continue
            name = line[:6] if number else line[:width]
            rows.append((name.strip(), orbital))
            if len(rows) == 10000:
                yield rows
                rows = []
        if rows:
            yield rows


def screen(rows, vectors, limit):
    names = np.array([r[0] for r in rows])
    epoch, a, ecc, inc, peri, node, mean = np.array([r[1] for r in rows]).T
    inc, peri, node, mean = np.deg2rad([inc, peri, node, mean])
    ci, si = np.cos(inc), np.sin(inc)
    co, so = np.cos(node), np.sin(node)
    cw, sw = np.cos(peri), np.sin(peri)
    p = np.array([cw*co-sw*ci*so, cw*so+sw*ci*co, sw*si]).T
    q = np.array([-sw*co-cw*ci*so, -sw*so+cw*ci*co, cw*si]).T
    plane = np.array([si*so, -si*co, ci]).T
    n = .01720209895 / a**1.5
    b = a*np.sqrt(1-ecc**2)
    perih, aph = a*(1-ecc), a*(1+ecc)
    mind = np.full(len(a), limit)
    when = np.zeros(len(a))
    for day, *state in vectors:
        comet_x, comet_v = np.array(state).reshape(2,3)
        r = np.linalg.norm(comet_x)
        if r >= 6:
            continue
        i = np.where((r >= perih-limit-.12) & (r <= aph+limit+.12) &
                     (np.abs(plane @ comet_x) <= limit+.12))[0]
        if not len(i):
            continue
        m = (mean[i]+(day-epoch[i])*n[i]+np.pi)%(2*np.pi)-np.pi
        E = m + np.where(ecc[i]>.65, np.copysign(.7,m), 0.)
        for _ in range(12):
            change = (E-ecc[i]*np.sin(E)-m)/(1-ecc[i]*np.cos(E))
            E -= change
            if np.max(abs(change)) < 1e-9:
                break
        sn, cs = np.sin(E), np.cos(E)
        pos = p[i]*(a[i]*(cs-ecc[i]))[:,None] + q[i]*(b[i]*sn)[:,None] - comet_x
        vel = (-p[i]*(a[i]*sn)[:,None] + q[i]*(b[i]*cs)[:,None])*(n[i]/(1-ecc[i]*cs))[:,None] - comet_v
        off = np.clip(-np.sum(pos*vel,axis=1)/np.sum(vel*vel,axis=1),-1,1)
        d = np.linalg.norm(pos+vel*off[:,None],axis=1)
        ok = d < mind[i]
        mind[i[ok]], when[i[ok]] = d[ok], day+off[ok]
    return zip(names[when!=0], mind[when!=0], when[when!=0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('numbered')
    parser.add_argument('unnumbered')
    parser.add_argument('--vectors', default=Path(__file__).parent / 'data/borisov-screen.csv')
    args = parser.parse_args()
    vectors = np.loadtxt(args.vectors, delimiter=',', skiprows=1)
    writer = csv.writer(sys.stdout)
    writer.writerow(['object', 'two_body_distance_au', 'mjd_tdb'])
    for path, width, number in [(args.numbered, 25, True), (args.unnumbered, 14, False)]:
        for rows in chunks(path, width, number):
            writer.writerows(screen(rows, vectors, .25))
