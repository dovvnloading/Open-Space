from astropy.io import fits
from pathlib import Path
import json,urllib.request,hashlib,numpy as np,csv
p=Path(r'D:\Open-Space\campaigns\dr20-clagn-e657b57-terra61\inputs');out=p.parent/'stars';out.mkdir(exist_ok=True);tables={}
for m in (59648,60322):
 with fits.open(p/f'zbest-{m}.prefix',ignore_missing_end=True,memmap=False) as h:
  tables[m]={int(r['SDSS_ID']):int(r['CATALOGID']) for r in h[1].data if int(r['SDSS_ID'])>0 and str(r['OBJTYPE']).strip()=='SPECTROPHOTO_STD'}
ids=sorted(tables[59648].keys()&tables[60322].keys())[:8];print('ids',ids,flush=True);receipts=[];measures=[]
for sid in ids:
 for m in (59648,60322):
  cid=tables[m][sid];f=f'spec-112359-{m}-{cid}.fits';u=f'https://dr20.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/spectra/daily/lite/112XXX/112359/{m}/{f}';b=urllib.request.urlopen(u,timeout=30).read();(out/f).write_bytes(b);receipts.append(dict(id=sid,mjd=m,url=u,sha256=hashlib.sha256(b).hexdigest()))
  with fits.open(out/f) as h:
   c=h['COADD'].data;s=h['SPALL'].data[0];assert int(s['SDSS_ID'])==sid
   wave=10**c['loglam'];good=(c['ivar']>0)&np.isfinite(c['flux'])&(c['and_mask']==0)
   for center in [5500,6500,7500,8500]:
    region=(wave>=center-50)&(wave<center+50);mask=region&good;coverage=mask.sum()/region.sum() if region.sum() else 0;flux=float(np.median(c['flux'][mask])) if coverage>=.8 else float('nan');measures.append(dict(sdss_id=sid,mjd=m,wavelength=center,coverage=round(float(coverage),4),median_flux=flux))
with open(r'results\dr20-clagn-pilot\calibration-stars.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=measures[0]);w.writeheader();w.writerows(measures)
Path(r'results\dr20-clagn-pilot\calibration-receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')
for center in [5500,6500,7500,8500]:
 ratios=[]
 for sid in ids:
  a=[x['median_flux'] for x in measures if x['sdss_id']==sid and x['mjd']==59648 and x['wavelength']==center][0];b=[x['median_flux'] for x in measures if x['sdss_id']==sid and x['mjd']==60322 and x['wavelength']==center][0]
  if a>0 and np.isfinite(b):ratios.append(b/a)
 med=np.median(ratios);mad=1.4826*np.median(np.abs(np.asarray(ratios)-med));print(center,len(ratios),float(med),float(mad))
