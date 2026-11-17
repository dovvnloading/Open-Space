"""Aperture-corrected MIRI forced photometry and Sonora Bobcat conditional model comparison."""
import argparse,csv,io,json,sys,tarfile,hashlib
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from reduce import HERE,R,XX,YY,HALF,prepare,radial_reference,core,save_table

def ap(im,xc,yc,rad,skyin,skyout):
    sub=np.arange(5)/5-.4
    weight=np.mean([np.hypot(XX-xc+dx,YY-yc+dy)<=rad for dx in sub for dy in sub],axis=0)
    sep=np.hypot(XX-xc,YY-yc)
    sky=np.median(im[(sep>=skyin)&(sep<=skyout)])
    return np.sum(weight*(im-sky))
def main(cache,out):
    meta=json.loads((HERE.parent/'miri'/'calibration_inputs.json').read_text())
    for r in meta:
        p=cache/r['filename']
        if not p.exists():p.write_bytes(urlopen(r['url'],timeout=120).read())
        if hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise RuntimeError(p)
    cal=fits.getdata(cache/'jwst_miri_apcorr_0014.fits',1)
    rows=json.loads((HERE/'inputs.json').read_text())
    peaks=list(csv.DictReader((out/'peaks.csv').open()))
    records=[]; fluxes={}
    for source in peaks:
        b=source['filter']; rr=[r for r in rows if r['filters']==b]
        s=prepare(cache/rr[0]['filename']); std=[prepare(cache/r['filename']) for r in rr[1:]]
        ref=radial_reference(std); use=s['valid']&(R<=14)
        im=s['image']
        def cost(v): return im[use]-v[0]*ref[use]-v[1]
        v=least_squares(cost,[core(im).sum(),np.median(im[use])],loss='soft_l1',f_scale=2).x
        im=im-v[0]*ref-v[1]
        x,y=np.linalg.solve(s['mat'],[float(source['east_arcsec']),float(source['north_arcsec'])])
        for ee in [.7,.8]:
            row=next(c for c in cal if c['filter']==b and c['subarray']=='SUB256' and abs(c['eefraction']-ee)<.001)
            rad,si,so,corr=[float(row[k]) for k in ('radius','skyin','skyout','apcorr')]
            conversion=s['solid']*1e12*corr
            flux=ap(im,x,y,rad,si,so)*conversion
            blank=[]
            for phi in np.arange(0,360,4):
                a=np.deg2rad(phi); dx=np.hypot(x,y)*np.cos(a);dy=np.hypot(x,y)*np.sin(a)
                if np.hypot(dx-x,dy-y)<2*so or not np.all(s['valid'][np.hypot(XX-dx,YY-dy)<so+1]):continue
                blank.append(ap(im,dx,dy,rad,si,so)*conversion)
            center=np.median(blank);err=1.4826*np.median(np.abs(np.array(blank)-center))
            records.append(dict(filter=b,encircled_energy=ee,aperture_pix=rad,sky_in_pix=si,sky_out_pix=so,apcorr=corr,
                                corrected_flux_uJy=flux,blank_median_uJy=center,blank_MAD_sigma_uJy=err,blank_sites=len(blank)))
            if ee==.7:fluxes[b]=[flux,err]
    save_table(out/'source_photometry.csv',records)
    result=[]; distance=1000/57.062001208183894
    grid=[]
    with tarfile.open(cache/'bobcat.tar.gz') as tar:
        for mh in ['-0.5','+0.0','+0.5']:
            lines=tar.extractfile('photometry_tables/flux_table_JWST'+mh).read().decode().splitlines()
            i=next(i for i,s in enumerate(lines) if s.strip().startswith('Teff'))
            names=lines[i].split();bands=names[names.index('F070W'):]
            tab=np.loadtxt(io.StringIO('\n'.join(s for s in lines[i+1:] if s.strip() and s.strip()[0].isdigit() and '*' not in s)))
            optical=tar.extractfile('photometry_tables/flux_table'+mh).read().decode().splitlines()
            j=next(j for j,s in enumerate(optical) if s.strip().startswith('Teff'))
            heading=optical[j].split()
            irac_index=heading.index('[4.5]')-len(heading)  # index from the end, independent of multiword header labels
            irac={}
            for line in optical[j+1:]:
                if line.strip() and line.strip()[0].isdigit() and '*' not in line:
                    numbers=list(map(float,line.split()))
                    irac[tuple(numbers[:6])]=1000*(10/distance)**2*10**numbers[irac_index]
            for model in tab:
                ch=[1000*(10/distance)**2*10**model[6+bands.index(b)] for b in ['F444W','F2100W','F2550W']]
                residual=[(flux-fluxes[b][0])/fluxes[b][1] for flux,b in zip(ch[1:],['F2100W','F2550W'])]
                grid.append(dict(metallicity=mh,teff=int(model[0]),logg=model[1],F444W_uJy=ch[0],IRAC45_uJy=irac.get(tuple(model[:6]),np.nan),F2100W_uJy=ch[1],F2550W_uJy=ch[2],
                                 rms_sigma=np.sqrt(np.mean(np.square(residual))),maxabs_sigma=max(abs(t) for t in residual)))
    grid.sort(key=lambda t:t['rms_sigma'])
    save_table(out/'model_best_and_2sigma.csv',[g for i,g in enumerate(grid) if i<5 or g['maxabs_sigma']<=2])
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();main(a.cache,a.out)
