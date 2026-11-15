"""Illustrative transverse-motion tests with the 2022 ASKAP component as the origin."""
from datetime import datetime
import sys
import numpy as np
from replay import read,chord
from collections import defaultdict

tracks=defaultdict(list)
for row in read('ephemerides.csv'):
    if row['date']=='2054-03-02' and row['ephemeris']=='JPL' and row['site'] in ('MeerKAT','GMRT'):
        tracks[row['site']].append(row)
pos=next(p for p in read('positions.csv') if p['position']=='RACS-low2')
ra,dec=float(pos['ra_deg']),float(pos['dec_deg'])
print('motion_mas_per_yr,position_angles_deg,angles_with_both_chords,angles_with_overlap,min_common_seconds,max_common_seconds')
def seconds(s):
    t=datetime.strptime(s,'%H:%M:%S.%f')
    return t.hour*3600+t.minute*60+t.second+t.microsecond*1e-6
for pm in (0,50,100,150,200,250):
    chord_count=0
    overlaps=[]
    for angle in np.deg2rad(np.arange(0,360,5)):
        radius=pm*32/1000
        r=ra+radius*np.sin(angle)/(3600*np.cos(np.deg2rad(dec)))
        d=dec+radius*np.cos(angle)/3600
        events=[chord(rows,r,d) for rows in tracks.values()]
        if all(e[3] for e in events):
            chord_count+=1
            overlaps.append(max(0.,min(seconds(e[4]) for e in events)-max(seconds(e[3]) for e in events)))
        else:
            overlaps.append(0.)
    print(f'{pm},72,{chord_count},{sum(x>0 for x in overlaps)},{min(overlaps):.1f},{max(overlaps):.1f}')
