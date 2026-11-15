"""Offline solid-disk chord check for two pulsar sightlines. Python 3, NumPy, SciPy."""
import csv
from datetime import datetime,timedelta
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parent

def read(name):
    with (ROOT/name).open() as stream:
        return list(csv.DictReader(stream))

def chord(rows, ra, dec):
    stamps = [datetime.fromisoformat(r['utc']) for r in rows]
    seconds = np.array([(t - stamps[0]).total_seconds() for t in stamps])
    pos = np.deg2rad([[float(r['ra_deg']), float(r['dec_deg'])] for r in rows])
    a, d = pos.T
    xyz = np.column_stack([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
    path = CubicSpline(seconds,xyz)
    radius = CubicSpline(seconds,[float(r['radius_arcsec']) for r in rows])
    a,d = np.deg2rad([ra,dec])
    target = np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
    def separation(t):
        v = path(t)
        return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(v,target)),v@target))*3600
    closest = minimize_scalar(separation,bounds=(seconds[0],seconds[-1]),method='bounded')
    t=closest.x
    def stamp(offset):
        return (stamps[0]+timedelta(seconds=float(offset))).strftime('%H:%M:%S.%f')[:-5]
    gap = lambda x: separation(x)-radius(x)
    ingress,egress='',''
    if gap(t)<0:
        ingress=stamp(brentq(gap,seconds[0],t))
        egress=stamp(brentq(gap,t,seconds[-1]))
    return stamp(t),round(closest.fun,3), round(float(radius(t)),3), ingress,egress

def results():
    positions=read('positions.csv')
    tracks={}
    for r in read('ephemerides.csv'):
        tracks.setdefault((r['date'],r['site'],r['ephemeris']),[]).append(r)
    out=[]
    for (date,site,model),track in tracks.items():
        for p in positions:
            if p['date']!=date:continue
            result=chord(track,float(p['ra_deg']),float(p['dec_deg']))
            alt=[float(r['alt_deg']) for r in track if r['alt_deg']]
            out.append((date,site,model,p['position'],*result,round(min(alt),1) if alt else '',round(max(alt),1) if alt else ''))
    return out

if __name__=='__main__':
    import sys
    w=csv.writer(sys.stdout,lineterminator='\n')
    w.writerow(['date','site','ephemeris','position','closest_utc','separation_arcsec','radius_arcsec','ingress_utc','egress_utc','alt_min_over_grid_deg','alt_max_over_grid_deg'])
    w.writerows(results())
