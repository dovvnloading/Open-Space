"""Targeted VizieR cone requests around COSMOS-022; retain responses for exact provenance."""
import re
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
D=Path(__file__).resolve().parent
cats={'vst2015':'J/A+A/574/A112','vst2022':'J/A+A/664/A117','chandra':'J/ApJ/817/34','cosmos2015':'J/ApJS/224/24','3dhst':'J/ApJS/214/24/3dhstall','xmm':'J/ApJ/716/348','smuvs':'J/A+A/633/A159','vlaradio':'J/A+A/602/A1','mirradioagn':'J/A+A/602/A2','milliquas':'VII/294'}
def run(kv):
 name,cat=kv
 q=urlencode({'-source':cat,'-c':'150.122033 +2.323386208','-c.rs':3 if name in ('3dhst','cosmos2015') else 100,'-out.add':'_r','-sort':'_r','-out.max':8})
 for host in ('vizier.cfa.harvard.edu','vizier.cds.unistra.fr','vizier.u-strasbg.fr'):
  url='https://'+host+'/viz-bin/asu-tsv?'+q
  try:
   with urlopen(url,timeout=30) as r:b=r.read()
   t=b.decode(errors='replace')
   if '\n#Table' not in t or 'STDB_ERR' in t:continue
   (D/(name+'.tsv')).write_bytes(b)
   print(name,len(b),[s[:225] for s in t.splitlines() if s.lstrip()[:1].isdigit()][:2],flush=True)
   return dict(name=name,url=url,sha256=sha256(b).hexdigest(),bytes=len(b))
  except OSError:continue
 print('unavailable',name,flush=True);return dict(name=name,unavailable=True)
if __name__=='__main__':
 with ThreadPoolExecutor(6) as ex:r=list(ex.map(run,cats.items()))
 (D/'catalog_receipts.json').write_text(json.dumps(r,indent=1)+'\n')
