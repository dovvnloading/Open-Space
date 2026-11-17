"""Audit catalogue neighbours inside GALEX's beam and a 1-arcminute local background."""
import csv
import hashlib
import io
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

HERE=Path(__file__).parent
u=json.loads((HERE/'galex.json').read_text())['data'][0]
ra,dec=u['ra'],u['dec']
q=(f"SELECT source_id,ra,dec,pmra,pmdec,ref_epoch,phot_g_mean_mag,bp_rp "
   f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
   f"CIRCLE('ICRS',{ra:.9f},{dec:.9f},0.016666667)) ORDER BY source_id")
url='https://gea.esac.esa.int/tap-server/tap/sync?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=3000,QUERY=q))
data=urllib.request.urlopen(url,timeout=50).read()
rows=list(csv.DictReader(io.StringIO(data.decode())))
assert 0<len(rows)<3000 and 'phot_g_mean_mag' in rows[0]
for r in rows:
    dt=2007-float(r['ref_epoch'] or 2016)
    lat=float(r['dec'])+dt*float(r['pmdec'] or 0)/3600000
    lon=float(r['ra'])+dt*float(r['pmra'] or 0)/3600000/math.cos(math.radians(float(r['dec'])))
    r['sep_galex_approx2007_arcsec']=round(3600*math.hypot((lon-ra)*math.cos(math.radians(dec)),lat-dec),3)
rows.sort(key=lambda r:r['sep_galex_approx2007_arcsec'])
with (HERE/'nuv_neighbours.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
(HERE/'nuv_neighbours_input.json').write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(data).hexdigest(),rows=len(rows)),indent=2)+'\n')
# A larger UV cone estimates the rate of accidental AIS matches at the observed brightness.
request=dict(service='Mast.Galex.Catalog',params=dict(ra=ra,dec=dec,radius=0.10),format='json',pagesize=5000,page=1)
url_uv='https://mast.stsci.edu/api/v0/invoke?'+urllib.parse.urlencode(dict(request=json.dumps(request)))
data_uv=urllib.request.urlopen(url_uv,timeout=90).read()
payload=json.loads(data_uv)
if payload['status']!='COMPLETE' or payload['paging']['pagesFiltered']>1:raise ValueError('incomplete local GALEX cone')
clean=[s for s in payload['data'] if s['survey']=='AIS' and 0<s['nuv_mag']<=u['nuv_mag']+u['nuv_magerr']
       and s['nuv_artifact']==0 and s['nuv_exptime']>0 and s['distance_arcmin']>8/60]
centres=[]
for s in sorted(clean,key=lambda s:s['nuv_mag']):
    xy=(s['ra']*math.cos(math.radians(dec))*3600,s['dec']*3600)
    if all(math.dist(xy,c)>3 for c in centres):centres.append(xy)
area=math.pi*(360**2-8**2)
rho=len(centres)/area
from scipy.stats import chi2
upper_rho=.5*chi2.ppf(.95,2*(len(centres)+1))/area
upper={str(r):-math.expm1(-upper_rho*math.pi*r*r) for r in [.74,6,8]}
close=[str(s['objID']) for s in payload['data'] if s['distance_arcmin']<=8/60]
(HERE/'nuv_chance.json').write_text(json.dumps(dict(url=url_uv,sha256=hashlib.sha256(data_uv).hexdigest(),
    all_rows=len(payload['data']),qualifying_detections=len(clean),unique_outside_8_arcsec=len(centres),
    area_arcsec2=area,object_ids_inside_8_arcsec=close,one_sided_95_upper_density_arcsec2=upper_rho,
    p_chance_uniform_95_upper=upper),indent=2)+'\n')
