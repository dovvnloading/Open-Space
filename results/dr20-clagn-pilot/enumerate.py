"""Metadata-only candidate enumeration; no spectrum flux is read."""
from pathlib import Path
import csv,sys
from astropy.io import fits
p=Path(sys.argv[1]);out=Path(sys.argv[2]);tables={}
for file in p.glob('zbest-*.prefix'):
 m=int(file.stem.split('-')[1])
 with fits.open(file,ignore_missing_end=True,memmap=False) as h:
  tables[m]={int(r['SDSS_ID']):dict(catalogid=int(r['CATALOGID']),z=float(r['Z']),cls=str(r['CLASS']).strip(),zw=int(r['ZWARNING']),ra=float(r['FIBER_RA']),dec=float(r['FIBER_DEC'])) for r in h[1].data if int(r['SDSS_ID'])>0}
rows=[]
for sid,b in sorted(tables[59648].items()):
 if not(.1<=b['z']<=.75 and b['cls'] in ('QSO','GALAXY') and b['zw']==0 and sid in tables[60674]):continue
 early=[m for m in tables if 60156<=m<60674 and sid in tables[m]]
 if not early:continue
 dev=min(early);late=[m for m in tables if m>dev and sid in tables[m]]
 if not late:continue
 hold=max(late)
 rows.append(dict(sdss_id=sid,z=b['z'],baseline_mjd=59648,baseline_catalogid=b['catalogid'],dev_mjd=dev,dev_catalogid=tables[dev][sid]['catalogid'],hold_mjd=hold,hold_catalogid=tables[hold][sid]['catalogid'],baseline_ra=b['ra'],baseline_dec=b['dec']))
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('provisional',len(rows),'selected',min(32,len(rows)))
