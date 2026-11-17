from pathlib import Path
import csv,sys
from astropy.io import fits
p=Path(sys.argv[1]); src=list(csv.DictReader(open(sys.argv[2],newline='')))[:32];tables={}
for f in p.glob('zbest-*.prefix'):
 m=int(f.stem.split('-')[1])
 if not 60156<=m<=60498:continue
 with fits.open(f,memmap=False,ignore_missing_end=True) as h:tables[m]={int(r['SDSS_ID']):int(r['CATALOGID']) for r in h[1].data if int(r['SDSS_ID'])>0}
for r in src:
 sid=int(r['sdss_id']);good=[m for m in tables if sid in tables[m] and m>int(r['dev_mjd'])]
 r['hold_mjd']=max(good) if good else '';r['hold_catalogid']=tables[max(good)][sid] if good else ''
with open(sys.argv[3],'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=src[0].keys());w.writeheader();w.writerows(src)
print('selected',len(src),'with good later',sum(bool(r['hold_mjd']) for r in src))
