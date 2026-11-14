"""Compute visibility and conditional eclipse phase, with only archive timestamp metadata."""
import csv, datetime, json
from pathlib import Path
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.utils import iers
from astropy import units as u
import numpy as np
D=Path(__file__).resolve().parent
iers.conf.auto_download=False
field=SkyCoord(ra=21.06275*u.deg,dec=-73.32594444*u.deg)
location=EarthLocation.from_geodetic(-70.815*u.deg,-30.165*u.deg,2215*u.m)
p=np.array(json.loads((D/'model.json').read_text())['params'])
archive={r['archive_filename']:r for r in json.loads((D/'meta-final.json').read_text())[1:]}
rows=[]
for r in json.loads((D/'selection.json').read_text()):
    stamp=r['file_dateobs'].split('+')[0]
    # Public Astro Archive SIA gives file_dateobs. Use the timestamp directly, without assuming exposure duration.
    t=Time(stamp,scale='utc',location=location)
    hjd=t.utc.jd+t.light_travel_time(field,kind='heliocentric').to_value(u.day)
    year=(hjd-2456000.)/365.25
    phi=((hjd-2460000.0369)/1.1696719)%1
    ang=2*np.pi*year/p[12];shift=p[10]*np.cos(ang)+p[11]*np.sin(ang)
    c=p[3]+p[8]*year+shift
    d=(phi-c+.5)%1-.5
    b=r['archive_filename'].split('/')[-1]
    rows.append({'filename':b,'dateobs_utc':stamp,'hjd_utc_at_archive_timestamp':hjd,'filter':b.split('_')[-2],
                 'catalog_phase':phi,'model_secondary_phase':c,'offset_secondary_minutes':d*1.1696719*1440,
                 'model_sigma_minutes':p[5]*1.1696719*1440,'archive_hdu_url':r['url'],'release_date':archive.get(r['archive_filename'],{}).get('release_date',''),
                 'exposure_seconds':archive.get(r['archive_filename'],{}).get('exposure','')})
with (D/'phases.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for date in sorted(set(r['dateobs_utc'][:10] for r in rows)):
    group=sorted([r for r in rows if r['dateobs_utc'].startswith(date)],key=lambda r:r['hjd_utc_at_archive_timestamp'])
    print(date,len(group),[(j['filename'], round(j['catalog_phase'],6),round(j['offset_secondary_minutes'],1)) for j in [group[0],group[-1]]])
