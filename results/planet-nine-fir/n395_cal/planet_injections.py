"""Inject calibrated Gaia stellar profiles through the existing 2025 N395 detection filter."""
import argparse
import csv
import json
import gzip
import io
import hashlib
from urllib.request import urlopen
from scipy.constants import astronomical_unit
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as nd, signal

ROOT=Path(__file__).resolve().parent


def star_profile(im,x,y):
    ix,iy=round(x),round(y)
    patch=im[iy-26:iy+27,ix-26:ix+27].astype(float)
    yy,xx=np.indices(patch.shape)
    x,y=x-ix+26,y-iy+26
    rr=np.hypot(xx-x,yy-y)
    patch-=np.median(patch[(rr>16)&(rr<24)])
    for _ in range(3):
        ww=np.maximum(patch,0)*(np.hypot(xx-x,yy-y)<7)
        x,y=(ww*xx).sum()/ww.sum(),(ww*yy).sum()/ww.sum()
    centered=nd.shift(patch,(26-y,26-x),order=3)
    profile=np.maximum(centered[11:42,11:42],0.)
    rad=np.hypot(*np.mgrid[-15:16,-15:16])
    profile*=rad<=13
    return profile/profile.sum()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--optical-dir',type=Path,required=True)
    ap.add_argument('--orbit-dir',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True)
    args=ap.parse_args()
    sys.path[:0]=[str(args.optical_dir),str(args.orbit_dir)]
    from decam_difference import sky_subtract,psf_stars,sample,mad,register,placed
    import propagate
    stars=list(csv.DictReader((ROOT/'standards.csv').open()))
    rng=np.random.default_rng(50723)
    states=propagate.family(json.loads((args.orbit_dir/'orbit_input.json').read_text()))
    mags=np.arange(21.,24.51,.25)
    model_rng=np.random.default_rng(507231)
    model_trials=[]
    solar_file=args.cache/'sun_reference_stis_002.fits'
    if not solar_file.exists():
        with urlopen('https://ssb.stsci.edu/trds/calspec/sun_reference_stis_002.fits',timeout=90) as response:
            solar_file.write_bytes(response.read())
    with fits.open(solar_file) as sun:
        sw=np.asarray(sun[1].data['WAVELENGTH'],float)
        sf=np.asarray(sun[1].data['FLUX'],float)
    rw,trans=np.genfromtxt(ROOT/'bandpass.csv',delimiter=',',skip_header=1).T
    wave=np.linspace(rw.min()*10,rw.max()*10,2048)
    weight=np.interp(wave,rw*10,trans)*wave
    sunmag=-2.5*np.log10(np.trapezoid(np.interp(wave,sw,sf)*weight,wave)/np.trapezoid(3631e-23*2.99792458e18/wave**2*weight,wave))
    with (args.orbit_dir/'brightness.csv').open() as f:
        scenarios=sorted({(float(r['mass_Mearth']),float(r['radius_Rearth']),float(r['teff_K']),float(r['geometric_albedo_r'])) for r in csv.DictReader(f)})
    results=[]; profiles=[]; candidates=[]
    with fits.open(args.cache/'recent_n395.response') as fi, fits.open(args.cache/'recent_n395_mask.response') as fd:
        hdr=fi[0].header
        start=propagate.epoch(hdr['DATE-OBS'])
        duration=hdr['EXPTIME']/86400
        ra,dec=propagate.sky(states,start+duration/2)
        pos,vel,t0,_=states
        helobs,baryobs=propagate.observers(start+duration/2)
        position=propagate.kepler(pos,vel,start+duration/2-t0-propagate.norm(pos-baryobs)/propagate.C)
        helio=propagate.norm(position-baryobs+helobs)
        delta=propagate.norm(position-baryobs)
        a0,d0=propagate.sky(states,start)
        a1,d1=propagate.sky(states,start+duration)
        ids=rng.choice(len(ra),320,replace=False)
        for ccd in [15,21]:
            im=fi[ccd].data;dq=fd[ccd].data;h=fi[ccd].header;w=WCS(h)
            fit=[r for r in stars if r['selection']=='accepted' and int(r['ccd'])==ccd]
            zp=np.median([float(r['N395_AB'])+2.5*np.log10(float(r['ap13_adu'])) for r in fit])
            references=[star_profile(im,float(r['x']),float(r['y'])) for r in fit]
            xx,yy=w.world_to_pixel_values(ra,dec)
            left,bot=w.world_to_pixel_values(a0,d0);right,top=w.world_to_pixel_values(a1,d1)
            trailx=right-left;traily=top-bot
            ins=(xx>=40)&(yy>=40)&(xx<im.shape[1]-40)&(yy<im.shape[0]-40)
            xs=np.arange(max(0,int(min(xx[ins])-140)),min(im.shape[1],int(max(xx[ins])+140)))
            ys=np.arange(max(0,int(min(yy[ins])-140)),min(im.shape[0],int(max(yy[ins])+140)))
            xx-=xs[0];yy-=ys[0]
            bad=dq[np.ix_(ys,xs)]!=0
            im=np.nan_to_num(sky_subtract(np.where(bad,np.nan,im[np.ix_(ys,xs)])))
            tpl,th=fits.getdata(args.cache/f'dr11_n395ref_{ccd}.response',header=True)
            tb=register(tpl,WCS(th),w,xs,ys,1)
            bad|=~np.isfinite(tb)
            sm=nd.gaussian_filter(np.nan_to_num(tb),2)
            fixed=nd.maximum_filter(sm-np.median(sm)>5*mad(sm),23)
            sky_noise=mad(im[~bad])
            def sharp(im):return im-nd.median_filter(im,3)>8*sky_noise+.35*np.maximum(im,0)
            mask=nd.maximum_filter(bad,29)|nd.maximum_filter(sharp(im),29)|fixed
            mask[:32]=mask[-32:]=True;mask[:,:32]=mask[:,-32:]=True
            stars0,ps=psf_stars(im,im,bad);kernel=ps[0];norm=np.sum(kernel**2)
            fm=signal.fftconvolve(im,kernel[::-1,::-1],mode='same')/norm
            zero=np.median(fm[~mask]);noise=mad(fm[~mask])
            good=ins&(sample(mask.astype('f4'),xx,yy)==0)
            scored=(sample(fm,xx[good],yy[good])-zero)/noise
            candidates.append(dict(ccd=ccd,valid=int(good.sum()),above_5p5=int((scored>=5.5).sum()),max_snr=float(max(scored))))
            center=kernel.shape[0]//2
            rad=np.hypot(*np.mgrid[-center:center+1,-center:center+1])
            pr=np.hypot(*np.mgrid[-15:16,-15:16])
            profiles.append(dict(ccd=ccd,ZP_ap13_AB_total_ADU=float(zp),kernel_peak=float(kernel.max()),
                    stellar_peak_range=[float(min(p.max() for p in references)),float(max(p.max() for p in references))],
                    old_kernel_fraction_r2=float(kernel[rad<=2].sum()),
                    stellar_fraction_r2=[float(p[pr<=2].sum()) for p in references],
                    trail_length_px_range=[float(np.hypot(trailx[good],traily[good]).min()),float(np.hypot(trailx[good],traily[good]).max())]))
            gain=(h['GAINA']+h['GAINB'])/2
            for j in ids:
                if not good[j]:continue
                ix,iy=int(round(xx[j])),int(round(yy[j]))
                z=im[iy-31:iy+32,ix-31:ix+32];x=xx[j]-ix+31;y=yy[j]-iy+31
                k=rng.integers(len(references));stellar=references[k]
                stamp=np.zeros(z.shape)
                for alpha in np.linspace(-.5,.5,5):
                    placed(stamp,(x+alpha*trailx[j],y+alpha*traily[j]),stellar,.2)
                stamp=np.maximum(stamp,0);stamp/=stamp.sum()
                def plant(magnitude, generator):
                    flux=10**(.4*(zp-magnitude))
                    planted=z+generator.poisson(stamp*flux*gain)/gain
                    reject=sample(nd.maximum_filter(sharp(planted),29).astype('f4'),[x],[y])[0]>0
                    matched=signal.fftconvolve(planted,kernel[::-1,::-1],mode='same')/norm
                    snr=(sample(matched,[x],[y])[0]-zero)/noise
                    return dict(ccd=ccd,track=int(j),psf_gaia_id=fit[k]['id'],N395_AB=float(magnitude),
                                flux_total_ap13_ADU=flux,snr=float(snr),post_injection_sharp_reject=int(reject),recovered=int(snr>=5.5 and not reject))
                for magnitude in mags:
                    results.append(plant(magnitude,rng))
                for mass,radius,teff,albedo in scenarios:
                    magnitude=sunmag-2.5*np.log10(albedo*(radius*6371e3/astronomical_unit/(helio[j]*delta[j]))**2)
                    trial=plant(magnitude,model_rng)
                    model_trials.append(dict(mass_Mearth=mass,radius_Rearth=radius,teff_K=teff,geometric_albedo_N395=albedo,**trial))
    for name,table in [('injections.csv.gz',results),('planet_trials.csv.gz',model_trials)]:
        text=io.StringIO(newline='')
        writer=csv.DictWriter(text,table[0]);writer.writeheader();writer.writerows(table)
        (ROOT/name).write_bytes(gzip.compress(text.getvalue().encode(),mtime=0))
    model_results=[]
    for mass,radius,teff,albedo in scenarios:
        trials=[t for t in model_trials if t['radius_Rearth']==radius and t['geometric_albedo_N395']==albedo]
        mag=sunmag-2.5*np.log10(albedo*(radius*6371e3/astronomical_unit/(helio*delta))**2)
        model_results.append(dict(mass_Mearth=mass,radius_Rearth=radius,teff_K=teff,geometric_albedo_N395=albedo,
               min_grid_mag=mag.min(),max_grid_mag=mag.max(),eligible=len(trials),recovered=sum(t['recovered'] for t in trials)))
    with (ROOT/'model_recovery.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,model_results[0]);writer.writeheader();writer.writerows(model_results)
    summary={'sun_N395_AB_at_1au':float(sunmag),'solar_file_sha256':hashlib.sha256(solar_file.read_bytes()).hexdigest(),'ccds':profiles,'original_data_search':candidates,'injected_tracks':320,'by_magnitude':[]}
    for m in mags:
        rows=[r for r in results if r['N395_AB']==m]
        summary['by_magnitude'].append(dict(N395_AB=float(m),eligible=len(rows),recovered=sum(r['recovered'] for r in rows),
            ccd15_hits=sum(r['recovered'] for r in rows if r['ccd']==15),ccd21_hits=sum(r['recovered'] for r in rows if r['ccd']==21)))
    (ROOT/'injection_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':main()
