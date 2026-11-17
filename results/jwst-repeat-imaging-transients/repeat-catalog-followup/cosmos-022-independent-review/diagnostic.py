import argparse, hashlib, json, pathlib, csv, io, numpy as np
from urllib.request import Request, urlopen
from astropy.io import fits
from astropy.wcs import WCS
HERE=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--cache',type=pathlib.Path,required=True)
p.add_argument('--out',type=pathlib.Path,required=True)
p.add_argument('--fetch',action='store_true',help='Fetch missing frozen byte ranges')
a=p.parse_args();root=a.cache.resolve();out=a.out.resolve()
if root==HERE or HERE in root.parents or out==HERE or HERE in out.parents or out==root or root in out.parents:
 raise ValueError('Output must be outside the frozen package and input cache')
receipts=json.loads((HERE/'ranges.json').read_text());prior={r['exposure']:r for r in csv.DictReader((HERE/'exposures.csv').open())}
if len(receipts)!=34 or len(prior)!=34 or {r['exposure'] for r in receipts}!=set(prior): raise ValueError('Input membership mismatch')
# Verify every range before calculating or writing results. Cache cannot redefine the receipt.
for rec in receipts:
 if not rec['url'].startswith('https://stpubdata.s3.amazonaws.com/jwst/public/') or '/' in rec['exposure'] or '..' in rec['exposure']: raise ValueError('Invalid receipt')
 d=root/rec['exposure']
 for item in rec['ranges']:
  target=d/(str(item['start'])+'_'+str(item['bytes'])+'.bin')
  if target.exists(): data=target.read_bytes()
  elif a.fetch:
   end=item['start']+item['bytes']-1
   with urlopen(Request(rec['url'],headers={'Range':f"bytes={item['start']}-{end}"}),timeout=90) as response:
    if response.status!=206 or not response.headers['Content-Range'].startswith(f"bytes {item['start']}-{end}/"): raise ValueError('Range response mismatch')
    data=response.read()
  else: raise FileNotFoundError(target)
  if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']: raise ValueError(f'Range mismatch: {target}')
  if not target.exists(): d.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
# Fixed post-inspection diagnostic: 4 radii, offsets 0 and +/-0.02 arcsec in each axis;
# all flagged core pixels excluded, common .40-.58 arcsec annulus, 5x5 area sampling.
result=[]
for rec in receipts:
 r=rec['ranges']; d=root/rec['exposure']; get=lambda x:(d/(str(x['start'])+'_'+str(x['bytes'])+'.bin')).read_bytes()
 hdr=fits.Header.fromfile(io.BytesIO(get(r[1]))); w=WCS(hdr,naxis=2); xf,yf=w.world_to_pixel_values(150.122033,2.323386208)
 width=hdr['NAXIS1']; x0=max(0,int(xf)-30);x1=min(width,int(xf)+31);y0=max(0,int(yf)-30);y1=min(hdr['NAXIS2'],int(yf)+31)
 sci=np.frombuffer(get(r[2]),dtype='>f4').reshape(-1,width)[:,x0:x1];dh=fits.Header.fromfile(io.BytesIO(get(r[4]))); dq=np.frombuffer(get(r[5]),dtype='>i4').reshape(-1,width)[:,x0:x1].astype('i8');dq=(dq*dh.get('BSCALE',1)+dh.get('BZERO',0)).astype('u4')
 yy,xx=np.mgrid[y0:y1,x0:x1];scale=np.sqrt(hdr['PIXAR_SR'])*206264.806; good=np.isfinite(sci)&(dq==0)
 for east,north in [(0,0),(.02,0),(-.02,0),(0,.02),(0,-.02)]:
  ra=150.122033+east/3600/np.cos(np.deg2rad(2.323386208)); dec=2.323386208+north/3600; cx,cy=w.world_to_pixel_values(ra,dec)
  dist=np.hypot(xx-cx,yy-cy)*scale; skyvals=sci[(dist>=.40)&(dist<.58)&good]; sky=np.median(skyvals) if len(skyvals)>50 else np.nan
  for rad in (.12,.16,.24,.32):
   wt=np.zeros_like(sci,dtype=float)
   for sy in (-.4,-.2,0,.2,.4):
    for sx in (-.4,-.2,0,.2,.4): wt+=(np.hypot(xx+sx-cx,yy+sy-cy)*scale<rad)/25
   coverage=float(wt[good].sum()/wt.sum()); flux=float(((sci[good]-sky)*wt[good]).sum()*hdr['PIXAR_SR']*1e15)
   z=prior[rec['exposure']];result.append(dict(exposure=rec['exposure'],band=z['band'],mosaic=z['mosaic'],radius=rad,east=east,north=north,flux=flux,coverage=coverage,nsky=len(skyvals)))
out.mkdir(parents=True,exist_ok=True)
with (out/'diagnostics.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=result[0]);writer.writeheader();writer.writerows(result)
print('rows',len(result))
