"""Find pairs of consecutive eclipses visible in darkness from Chile and Australia."""
import csv
import json
import sys
import warnings
from pathlib import Path
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz,EarthLocation,SkyCoord,get_sun
from astropy.time import Time
from astropy.utils import iers
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'joint'))
from joint import geometry,initial
from orbit import EP,P,TC,Y
iers.conf.auto_download=False
iers.conf.iers_degraded_accuracy='warn'
STAR=SkyCoord(21.06275*u.deg,-73.32594444*u.deg)
SITES={'Cerro Tololo':EarthLocation.from_geodetic(-70.8065*u.deg,-30.1697*u.deg,2200*u.m),
       'Siding Spring':EarthLocation.from_geodetic(149.0685*u.deg,-31.2749*u.deg,1164*u.m)}
MODELS=('kepler','single_trough','cubic','sinusoid_offset')
fits=json.loads((HERE.parent/'joint/fits.json').read_text())['full']

def mid(which,cycle,center):
    params=np.array(fits[which]['parameters'])
    t=EP+(cycle+.5*center)*P
    for k in range(5):
        years=(t-TC)/Y
        change=geometry(which,params[2:2+len(initial(which)[0])],years)[0][center]
        t=EP+P*(cycle+.5*center+params[0]+params[1]*years/10+change)
    return t

def utc_for_hjd(hjd,site):
    t=Time(hjd,format='jd',scale='utc',location=site)
    for j in range(2):
        delay=t.light_travel_time(STAR,kind='heliocentric')
        t=Time(hjd-delay.to_value(u.day),format='jd',scale='utc',location=site)
    return t

def utc_str(jd):
    return Time(jd+0.5/1440,format='jd').strftime('%Y-%m-%d %H:%M')

def find():
    cycles=np.arange(np.ceil((Time('2027-08-01').jd-EP)/P),np.floor((Time('2027-12-31').jd-EP)/P))
    events={};epochs=[{m:mid(m,cycles,c) for m in MODELS} for c in (0,1)]
    for name,site in SITES.items():
        for c in (0,1):
            kepler=epochs[c]['kepler'];times=np.array(list(epochs[c].values()))
            # Make allowance for model spread and enough baseline on either side.
            low=times.min(axis=0)-120/1440
            high=times.max(axis=0)+120/1440
            us=utc_for_hjd(low,site);ue=utc_for_hjd(high,site)
            ut=utc_for_hjd(kepler,site)
            grid=np.linspace(us.jd,ue.jd,11)
            ts=Time(grid,format='jd',scale='utc',location=site);frame=AltAz(obstime=ts,location=site)
            obj=STAR.transform_to(frame).alt.deg.min(axis=0)
            sun=get_sun(ts).transform_to(frame).alt.deg.max(axis=0)
            events[name,c]=(us.jd,ue.jd,ut.jd,obj,sun)
    all_pairs=[]
    for i,cycle in enumerate(cycles):
        for s1 in SITES:
            for s2 in SITES:
                v,w=events[s1,0],events[s2,1]
                if min(v[3][i],w[3][i])<30 or max(v[4][i],w[4][i])>-15:continue
                all_pairs.append(dict(cycle=int(cycle),primary_site=s1,secondary_site=s2,
                        start_primary_utc=utc_str(v[0][i]),end_primary_utc=utc_str(v[1][i]),kepler_primary_utc=utc_str(v[2][i]),
                        start_secondary_utc=utc_str(w[0][i]),end_secondary_utc=utc_str(w[1][i]),kepler_secondary_utc=utc_str(w[2][i]),
                        min_target_alt_deg=round(min(v[3][i],w[3][i]),1),max_sun_alt_deg=round(max(v[4][i],w[4][i]),1),
                        kepler_interval_minutes=round((epochs[1]['kepler'][i]-epochs[0]['kepler'][i])*1440,1),
                        trough_interval_minutes=round((epochs[1]['single_trough'][i]-epochs[0]['single_trough'][i])*1440,1)))
    if all_pairs:
        with (HERE/'windows.csv').open('w') as out:
            w=csv.DictWriter(out,fieldnames=all_pairs[0]);w.writeheader();w.writerows(all_pairs)
    print(len(all_pairs))
    for row in sorted(all_pairs,key=lambda r:-r['min_target_alt_deg'])[:3]:print(row)
if __name__=='__main__':find()
