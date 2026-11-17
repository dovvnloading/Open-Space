"""Extract one DESI coadd and plot inverse-variance weighted absorption profiles."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from astropy.io import fits
import hashlib, shutil, sys, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
OUT=Path(__file__).resolve().parent
TMP=Path(sys.argv[1]);TMP.mkdir(exist_ok=True,parents=True)
source_id=sys.argv[2] if len(sys.argv)>2 else '1717410430709663616'
import csv
with (OUT/'leads.csv').open() as f: leads=list(csv.DictReader(f))
if source_id not in {r['source_id'] for r in leads}:
    with (OUT.parent/'queries/coadds.csv').open() as f:leads.extend(csv.DictReader(f))
if source_id not in {r['source_id'] for r in leads}:
    with (OUT.parent.parent/'moderate/desi/queries/coadds.csv').open() as f:leads.extend(csv.DictReader(f))
lead=next(row for row in leads if row['source_id']==source_id)
target=int(lead['targetid']); adopted_rv=float(lead['vrad'])
survey,program,h=lead['survey'],lead['program'],lead['healpix']
prefix=h[:-2]
p=f'{survey}-{program}-{h}'
base=f'https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/{survey}/{program}/{prefix}/{h}/'
out_fit=OUT if source_id=='1717410430709663616' else ((OUT.parent/'tangential' if source_id in ['6340918677332110720','2752545725694798336'] else OUT)/source_id)
if source_id in ['4435527217525195264','4550108286326235520']:out_fit=OUT.parent/'moderate'/source_id
out_fit.mkdir(exist_ok=True)
name='coadd-'+p+'.fits'
def download(name):
    path=TMP/name
    if not path.exists():
        with urlopen(base+name,timeout=150) as f,path.with_suffix('.part').open('wb') as g: shutil.copyfileobj(f,g)
        path.with_suffix('.part').rename(path)
    return path
with ThreadPoolExecutor(2) as pool: fits_path,check_path=pool.map(download,[name,f'redux_iron_healpix_{survey}_{program}_{prefix}_{h}.sha256sum'])
with fits_path.open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
assert any(l.split()[0]==digest and l.split()[-1].lstrip('*')==name for l in check_path.read_text().splitlines())
speed=299792.458
# Vacuum rest wavelengths in Angstrom; weighted absorption profiles use a continuum slope.
# H-alpha from H I in standard air 6562.80, H-beta 4861.33, Ca II 8498.02/8542.09/8662.14, converted to vacuum.
lines=[('Hβ','B',4862.683),('Hα','R',6564.614),('Ca II 8498','Z',8500.358),('Ca II 8542','Z',8544.438),('Ca II 8662','Z',8664.520)]
fig,axs=plt.subplots(2,3,figsize=(10,5.3),constrained_layout=True)
profiles=[];results=[]
with fits.open(fits_path) as hdus:
    idx=np.flatnonzero(hdus['FIBERMAP'].data['TARGETID']==target).item()
    for (label,arm,rest),ax in zip(lines,axs.flat):
        w=hdus[arm+'_WAVELENGTH'].data; flux=hdus[arm+'_FLUX'].data[idx];iv=hdus[arm+'_IVAR'].data[idx]; mask=hdus[arm+'_MASK'].data[idx]
        v=speed*(w/rest-1)
        ok=(abs(v-adopted_rv)<850)&(mask==0)&(iv>0)&np.isfinite(flux)&np.isfinite(iv)
        v,y,e=v[ok],flux[ok],1/np.sqrt(iv[ok])
        def model(x,c,s,depth,cen,width):
            return c+s*(x-adopted_rv)/1000-depth*np.exp(-.5*((x-cen)/width)**2)
        prof,pcov=curve_fit(model,v,y,sigma=e,absolute_sigma=True,p0=[np.median(y),0,max(np.median(y)*.2,.1),adopted_rv,100.],bounds=([-np.inf,-np.inf,0,adopted_rv-500,20],[np.inf,np.inf,np.inf,adopted_rv+500,450]),maxfev=4000)
        c,s,d,cen,width=prof; err=np.sqrt(pcov[3,3]);sig=d/np.sqrt(pcov[2,2]);results.append({'line':label,'vacuum_A':rest,'rv_kms':cen,'rv_formal_error_kms':err,'line_depth_snr':sig})
        for vi,yi,ei in zip(v,y,e):profiles.append((label,rest,vi,yi,ei))
        ax.errorbar(v,y,yerr=e,color='#657481',ls='none',marker='.',ms=3,lw=.5)
        grid=np.linspace(v.min(),v.max(),500);ax.plot(grid,model(grid,*prof),color='#982f34');ax.axvline(adopted_rv,color='#177b8d',ls='--',lw=.9)
        ax.set_title(label);ax.set_xlabel('Barycentric velocity (km/s)');ax.set_ylabel('Flux')
    axs.flat[-1].axis('off')
fig.suptitle(f'DESI {target}')
fig.savefig(out_fit/f'spectrum_{source_id}.png',dpi=140)
import csv
with (out_fit/'line_fits.csv').open('w') as f:
    cw=csv.DictWriter(f,fieldnames=results[0].keys());cw.writeheader();cw.writerows(results)
with (out_fit/'line_pixels.csv').open('w') as f:
    cw=csv.writer(f);cw.writerow(['line','rest_vacuum_A','velocity_kms','flux','flux_error']);cw.writerows(profiles)
(out_fit/'spectrum_provenance.json').write_text(json.dumps({'url':base+name,'sha256':digest,'bytes':fits_path.stat().st_size,'targetid':str(target),'row':idx},indent=2)+'\n')
for a in results:print(a)
