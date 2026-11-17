"""Published follow-up locations. Epochs are days since J2000 (UTC surrogate)."""
import csv
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'orbits'))
from propagate import earth, epoch, unit, norm, MU, OBLIQUITY, kepler, ray

def radec(v):
    v = np.atleast_2d(v)
    return np.degrees(np.arctan2(v[:,1], v[:,0])) % 360, np.degrees(np.arcsin(v[:,2] / norm(v)))

def act_sky(row,t):
    # Naess et al. section 1, Tables 3–4: heliocentric tangent-plane velocity.
    years = (t - (57688 - 51544.5)) / 365.25
    a = row['ra'] + row['vx'] * years / (60 * np.cos(np.radians(row['dec'])))
    d = row['dec'] + row['vy'] * years / 60
    return radec(row['r'] * unit(a,d) - earth(t))

def basis(lon, lat):
    a,d=np.radians([lon,lat]); return np.array([-np.sin(a),np.cos(a),0]),np.array([-np.cos(a)*np.sin(d),-np.sin(a)*np.sin(d),np.cos(d)])

def circular_family(name):
    """Circular bracket for unknown heading, illustrative and non-probabilistic."""
    if name=='RR':
        t0=epoch('1983-07-01'); lon,lat=13.62,71.64
        rot=np.array([[1,0,0],[0,np.cos(OBLIQUITY),-np.sin(OBLIQUITY)],[0,np.sin(OBLIQUITY),np.cos(OBLIQUITY)]])
        n=rot@unit(lon,lat); e,b=[rot@x for x in basis(lon,lat)]
        radii=np.arange(210,241,5); observed=False
    else:
        lon,lat,ds={'Chen1':(42.676,-15.021,'2007-07-28'),'Chen2':(45.297,-16.711,'2006-07-30')}[name]
        t0=epoch(ds+'T12:00:00'); n=unit(lon,lat); e,b=basis(lon,lat)
        radii=np.arange(400,901,50); observed=True
    rr,hh=np.meshgrid(radii,np.radians(np.arange(0,360,5))); rr=rr.ravel(); hh=hh.ravel()
    if observed: pos,_=ray(earth(t0),n,rr)
    else: pos=rr[:,None]*n
    heading=np.cos(hh)[:,None]*b+np.sin(hh)[:,None]*e
    v=heading-np.sum(heading*pos,axis=1)[:,None]*pos/rr[:,None]**2
    v=v/norm(v)[:,None]*np.sqrt(MU/rr)[:,None]
    return pos,v,t0

def other_sky(states,t):
    pos,vel,t0=states
    # Exact circular solution short enough for both an annulus and decades of extrapolation.
    r=norm(pos); w=np.sqrt(MU/r**3); angle=w*(t-t0)
    return radec(pos*np.cos(angle)[:,None]+vel*(np.sin(angle)/w)[:,None]-earth(t))

def act_rows():
    rows={}
    with (HERE/'act.csv').open() as f:
        for row in csv.DictReader(f):
            name=row.pop('id'); rows[name]={k:float(v) for k,v in row.items()}
    return rows

def main():
    with (HERE/'other_ephemerides.csv').open('w') as f:
        out=csv.writer(f);out.writerow(['utc','target','directions','ra_min','ra_max','dec_min','dec_max'])
        for year in range(2023,2027):
            for month in range(1,13):
                day=f'{year}-{month:02d}-01';t=epoch(day)
                for name,row in act_rows().items():
                    ra,dec=act_sky(row,t);out.writerow([day,name,1,*np.round([ra[0],ra[0],dec[0],dec[0]],6)])
                for name in ['RR','Chen1','Chen2']:
                    ra,dec=other_sky(circular_family(name),t);out.writerow([day,name,len(ra),*np.round([ra.min(),ra.max(),dec.min(),dec.max()],6)])
if __name__=='__main__':main()
