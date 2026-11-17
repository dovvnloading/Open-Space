"""Measure Uranus-relative point-source photometry in JWST's original integrations."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi
from link_visits import radial_psf,search,geometry
from independent_vet import save
import positive_control as pc
from photometry import diameter
ROOT=Path(__file__).resolve().parent


def run(args):
    positions=list(csv.DictReader((ROOT/'integration_check/native_centroids.csv').open()))
    originals=list(csv.DictReader((ROOT/'integration_check/integration_positions.csv').open()))
    frames=json.loads((ROOT/'inputs.json').read_text()); report=[];results={}
    with fits.open(args.flat) as hd:flat=hd['SCI'].data.copy()
    for mode,bits in [('strict_DQ',1|2|4|1024),('allow_corrected_jumps',1|2|1024)]:
        for j,frame in enumerate(frames):
            exp=frame['exposure_id'];root=exp.removesuffix('_cal')
            h=fits.getheader(args.files/(exp+'.fits'),'SCI')
            p=fits.getheader(args.files/(exp+'.fits'))
            if p['R_FLAT'].rsplit('/',1)[-1]!=args.flat.name:raise ValueError('incorrect CRDS flat')
            m=pc.metadata(args.cache,root)
            residual,x0,y0=pc.wide_images(args.cache,m,bits)
            ny,nx=residual.shape[1:];conv=1e12*h['PIXAR_SR']*h['PHOTMJSR']/flat[y0:y0+ny,x0:x0+nx]
            residual*=conv
            orig=np.array([[float(r['x']),float(r['y'])] for r in originals if r['file']==m['file'] and r['body']=='trial'])
            sub=[r for r in positions if int(r['j'])==j]
            measured=np.array([[float(r['x']),float(r['y'])] for r in sub]);index=np.array([int(r['k']) for r in sub])
            # Interpolate centroid residuals only to simulate the reference subtraction for missing integrations.
            off=measured-orig[index]
            xy=orig+np.array([np.interp(np.arange(len(orig)),index,v) for v in off.T]).T-[x0,y0]
            v=np.gradient(xy,axis=0)*m['sec']/(np.mean(np.diff(m['times']))*86400)
            scale=np.sqrt(h['PIXAR_SR'])*geometry.ARCSEC
            psf=radial_psf(args.ee,scale,'F150W2');plant=np.zeros_like(residual)
            for k,(cc,vel) in enumerate(zip(xy,v)):
                path=cc+np.linspace(-.5,.5,13)[:,None]*vel
                window,stamp=search.clipped_stamp((ny,nx),*search.trail_kernel(path[:,0],path[:,1],psf))
                plant[k][window]+=stamp
            def apertures(images):
                out=[]
                for k in index:
                    xx,yy=xy[k];ix,iy=np.rint([xx,yy]).astype(int)
                    gy,gx=np.mgrid[iy-16:iy+17,ix-16:ix+17];rr=np.hypot(gx-xx,gy-yy)
                    data=images[k,iy-16:iy+17,ix-16:ix+17]
                    sky=(rr>=10)&(rr<=15)&np.isfinite(data);ap=(rr<=7)&np.isfinite(data)
                    value=(np.sum(data[ap])-ap.sum()*np.median(data[sky])) if sky.sum()>100 and ap.sum()>45 else np.nan
                    out.append(value)
                return np.array(out)
            y=apertures(residual)
            # Symmetric source insertion probes the local throughput of each temporal median and aperture.
            response=[]
            for delta in [-.5,.5]:
                sim=residual+delta*plant
                sim-=np.nanmedian(sim,axis=0)
                response.append(apertures(sim))
            throughput=response[1]-response[0]
            ok=(throughput>.35)&(throughput<1.05)&np.isfinite(y)
            flux=y[ok]/throughput[ok]
            report.append(dict(mode=mode,exposure=exp,usable=int(ok.sum()),median_flux_uJy=round(float(np.median(flux)),4) if len(flux) else '',
                               median_response=round(float(np.median(throughput[ok])),4) if len(flux) else ''))
        med=np.array([float(r['median_flux_uJy']) for r in report if r['mode']==mode and r['usable']])
        if not len(med):
            results[mode]=dict(exposures=0,status='no integration passed aperture throughput cut')
            continue
        val=float(np.median(med));result=dict(flux_uJy=val,between_exposure_MAD_uJy=float(1.4826*np.median(abs(med-val))),
                                           exposures=len(med),AB_mag=float(23.9000656223-2.5*np.log10(val)))
        results[mode]=result
    c=search.load_npz(args.cube)
    solar=json.loads((ROOT.parent/'jwst6379_9482/f150w2_size.json').read_text())['solar_flux_at_1au_Jy']
    sun_mag=-2.5*np.log10(solar/3631)
    for result in results.values():
        if not result['exposures']:continue
        result['diameter_km_p004_p008_p016']=[float(diameter(result['AB_mag'],sun_mag,np.median(c['path_delta']),np.median(c['r_au']),a)) for a in [.04,.08,.16]]
    results['flat']={'url':'https://jwst-crds.stsci.edu/unchecked_get/references/jwst/'+args.flat.name,
                     'sha256':hashlib.sha256(args.flat.read_bytes()).hexdigest()}
    results['encircled_energy_sha256']=hashlib.sha256(args.ee.read_bytes()).hexdigest()
    save(ROOT/'integration_check/photometry.csv',report)
    (ROOT/'integration_check/photometry.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('cube','files','cache','flat','ee'):p.add_argument(k,type=Path)
    run(p.parse_args())
