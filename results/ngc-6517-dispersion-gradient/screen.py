"""Replay the initial exploratory ATNF 2.9.0 cluster screen; use its frozen DB."""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from analyze import ROOT, sexagesimal, summary

EXPECTED_SHA256 = 'bac32471fef3e60cfd734f1db85503a9bd0a565018cb98edff54475f87eb642f'


def screen(path):
    from hashlib import sha256
    data = path.read_bytes()
    if sha256(data).hexdigest() != EXPECTED_SHA256:
        raise ValueError('Initial screen used the specific ATNF 2.9.0 database named in SOURCES.md.')
    clusters = defaultdict(list)
    for block in data.decode('latin1').split('@--------'):
        d = {parts[0]: parts[1] for line in block.splitlines()
             if len(parts := line.split()) > 1 and not line.startswith('#')}
        match = re.search(r'GC:([^,\[]+)', d.get('ASSOC', ''))
        if not match or not all(k in d for k in ('RAJ', 'DECJ','DM')):
            continue
        ra, dec = d['RAJ'], d['DECJ']
        # This original screen used printed precision as a proxy for localization.
        if (ra.count(':') != 2 or dec.count(':') != 2 or
                len(ra.split(':')[-1]) < 4 or len(dec.split(':')[-1]) < 3):
            continue
        clusters[match[1]].append([15*sexagesimal(ra), sexagesimal(dec), float(d['DM'])])
    with (ROOT/'screen.csv').open('w') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['atnf_cluster', 'pulsars', 'preliminary_r_squared'])
        for gc, group in clusters.items():
            if len(group) < 5:
                continue
            a = np.array(group)
            pos = (a[:,:2]-a[:,:2].mean(axis=0))*[3600*np.cos(np.radians(a[0,1])),3600]
            writer.writerow([gc,len(group), f"{summary(pos,a[:,2])['r2']:.8f}"])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    screen(parser.parse_args().database)
