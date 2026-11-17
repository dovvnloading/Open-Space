from pathlib import Path
import json, math
from astropy.io import fits
from matplotlib.path import Path as Poly
p=Path(__file__).resolve().parent
ids=['102158276_013663548_VIS','102158276_028106832_VIS','102158276_343226670_VIS']; point=(math.radians(273.167),math.radians(64.3687)); result=[]
for ident in ids:
 with fits.open(p/'inputs'/('layer-'+ident+'.fits')) as hd:
  d=hd[1].data; rows=[r for r in d if Poly(r['POLYGON'].reshape(4,2)).contains_point(point)]; exposures=[{'calib_image':str(r['CALIB_IMAGE']),'detector_extension':str(r['EXTNAME']),'mjd_obs_float32':float(r['MJD-OBS'])} for r in rows]
 result.append({'id':ident,'layer_rows':len(d),'point_covering_rows':len(rows),'exposures':exposures})
for i in range(len(result)):
 for j in range(i+1,len(result)):
  assert not ({e['calib_image'] for e in result[i]['exposures']} & {e['calib_image'] for e in result[j]['exposures']})
(p/'layers-summary.json').write_text(json.dumps(result,indent=2)+'\n');print([(x['id'],x['point_covering_rows']) for x in result])
