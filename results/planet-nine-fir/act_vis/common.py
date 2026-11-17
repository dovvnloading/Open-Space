"""Observer and tangent-plane geometry for the Naess et al. ACT Table 4 peak 6."""
import csv,io,json,os,sys
from pathlib import Path
import numpy as np
import erfa
from astropy.time import Time
from scipy.interpolate import CubicSpline
from metadata import HERE
CACHE=Path(os.environ.get('ACT_VIS_CACHE', str(HERE.parent/'cache_act4_06')))
CACHE.mkdir(exist_ok=True,parents=True)
def shared():
    options=[HERE.parent/'orbits'/'euclid',HERE.parent.parent/'orbits'/'euclid']
    if 'PHAN_EUCLID' in os.environ:options.insert(0,Path(os.environ['PHAN_EUCLID']))
    return next(p for p in options if (p/'vis_pixels.py').exists())
EU=shared()
sys.path[:0]=[str(EU),str(EU.parent)]
from propagate import observers,unit
rows=sorted(csv.DictReader((HERE/'euclid_dr1.response').open()),key=lambda r:float(r['t_min']))
assert len(rows)==4
T=Time([.5*(float(r['t_min'])+float(r['t_max'])) for r in rows],format='mjd',scale='utc')
def track():
    buf=json.loads((HERE/'horizons.response').read_bytes())['result'].split('$$SOE')[1].split('$$EOE')[0].strip()
    eph=np.array([[float(line[k]) for k in (0,2,3,4)] for line in csv.reader(io.StringIO(buf))])
    sc=CubicSpline(eph[:,0],eph[:,1:],extrapolate=False)(T.tdb.jd)
    earth_helio,earth_bary=observers(T.tt.jd-2451545)
    sc_heliocentric=sc-earth_bary+earth_helio
    yr=(T.tt.jd-2400000.5-57688.)/365.25
    ra=69.03+3.9*yr/60/np.cos(np.radians(-21.10)); de=-21.10+1.4*yr/60
    u=300.*unit(ra,de).T-sc_heliocentric
    return np.rad2deg(np.arctan2(u[:,1],u[:,0]))%360,np.rad2deg(np.arcsin(u[:,2]/np.linalg.norm(u,axis=1)))
RA,DE=track()
RA0,DE0=RA.mean(),DE.mean()
COS=np.cos(np.radians(DE0))
PIXSCALE=.1
STEP=PIXSCALE/3600
HALF=625
GX=np.arange(-HALF,HALF+1);GY=GX
N=len(GX)
def world(x,y):return RA0+x*STEP/COS,DE0+y*STEP
def coord(ra,de):return (ra-RA0)*COS/STEP,(de-DE0)/STEP
MOTION=np.array(coord(RA,DE)).T
if __name__=='__main__':
 print(T.isot,RA,DE,MOTION)
