from pathlib import Path
import csv,gzip,hashlib,io,json,re,math
p=Path(__file__).parent/'ztf-dr24'
for sub in ('requests','ibe'):
 for r in json.loads((p/sub/'requests.json').read_text()):
  b=gzip.decompress((p/sub/(r['file']+'.gz')).read_bytes())
  assert len(b)==r['bytes'] and hashlib.sha256(b).hexdigest()==r['sha256']
def rows(name):
 return list(csv.DictReader(io.StringIO(gzip.decompress((p/'requests'/(name+'.response.gz')).read_bytes()).decode())))
g=rows('normal-g-all'); assert len(g)==90
assert {r['oid'] for r in g}=={'1795112200005698'}
assert all(int(r['catflags'])==0 for r in g)
assert all(not (59588 <= float(r['mjd']) <= 59708 or 60262 <= float(r['mjd']) <= 60382) for r in g)
assert max(math.hypot((float(r['ra'])-213.351140496)*math.cos(math.radians(53.0908107113)), float(r['dec'])-53.0908107113)*3600 for r in g)<2
print('90 clean historical g rows; range',min(float(r['mjd']) for r in g),max(float(r['mjd']) for r in g))
for n in ('normal-g-baseline','normal-g-development','normal-r-all','normal-r-baseline','normal-r-development'): assert not rows(n)
for name in ('baseline','development'):
 b=gzip.decompress((p/'ibe'/(name+'.tbl.gz')).read_bytes()).decode()
 n=int(re.search(r'RowsRetrieved\s*=\s*(\d+)',b).group(1)); assert n==({'baseline':3,'development':0}[name]); print(name,n)
print('Verified all response hashes and window results')
