"""Check frozen source/visit lineage without reading spectral flux."""
import collections,csv,hashlib,json,math,sys
from pathlib import Path
from astropy.io import fits
root=Path(sys.argv[1]);manifest=json.loads(Path(sys.argv[2]).read_text());selection=list(csv.DictReader(open(sys.argv[3],newline='')));groups=collections.defaultdict(dict)
for r in manifest:
 assert 'error' not in r and r['kind'] not in groups[r['sdss_id']];groups[r['sdss_id']][r['kind']]=r
assert len(groups)==len(selection)==32 and set(groups)=={int(x['sdss_id']) for x in selection}
for sid,g in groups.items():
 assert set(g)=={'baseline','dev','hold'}
 for k,r in g.items():
  assert r['warning']==0 and len(r['tai_list'])>0 and len(r['tai_list'])==len(r['mjd_list'])
  assert hashlib.sha256((root/r['filename']).read_bytes()).hexdigest()==r['sha256']
  with fits.open(root/r['filename'],memmap=True) as h:
   s=h['SPALL'].data[0];assert int(s['SDSS_ID'])==sid and str(s['OBJTYPE']).strip().lower()=='science'
 for i,k in enumerate(('baseline','dev','hold')):
  for l in ('baseline','dev','hold')[i+1:]:
   a,b=g[k],g[l];assert set(a['tai_list']).isdisjoint(b['tai_list']);d=math.hypot((a['ra']-b['ra'])*math.cos(math.radians(a['dec'])),a['dec']-b['dec'])*3600;assert d<1
print('32 sources, 96 spectra: science, no redshift warning, matched coordinates and disjoint parent TAI lists')
