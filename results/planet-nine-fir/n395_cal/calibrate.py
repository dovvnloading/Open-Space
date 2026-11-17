"""Place the existing N395 injected-source ADU on an AB scale using Gaia XP spectra."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage
from scipy.constants import c

ROOT=Path(__file__).resolve().parent


def synthetic(wavelength,throughput,flux_lambda):
    """Photon-counted AB mean flux in Jy; Gaia F_lambda in W/(m² nm)."""
    fine=np.linspace(wavelength[0],wavelength[-1],2048)
    response=np.interp(fine,wavelength,throughput)
    sampled=np.interp(fine,np.arange(336.,1022.,2.),flux_lambda)
    return np.trapezoid(sampled*fine*response,fine)/np.trapezoid(response/fine,fine)/(c*1e9)*1e26


def photometry(im,dq,h,x,y,kernel):
    ix,iy=round(x),round(y)
    cut=im[iy-26:iy+27,ix-26:ix+27].astype(float)
    mask=dq[iy-26:iy+27,ix-26:ix+27]!=0
    yy,xx=np.indices(cut.shape)
    dx,dy=x-ix+26,y-iy+26
    rad=np.hypot(xx-dx,yy-dy)
    sky=np.median(cut[(rad>16)&(rad<24)&~mask])
    cut-=sky
    peak=np.max(cut[rad<6])+sky
    # Recenter on the image centroid, starting from Gaia proper-motion astrometry.
    for _ in range(3):
        rad=np.hypot(xx-dx,yy-dy)
        weights=(rad<7)*np.maximum(cut,0.)
        dx,dy=np.sum(weights*xx)/weights.sum(),np.sum(weights*yy)/weights.sum()
    drift=np.hypot(dx-(x-ix+26),dy-(y-iy+26))
    cy,cx=int(round(dy)),int(round(dx))
    rad=np.hypot(xx-dx,yy-dy)
    ap9=cut[rad<=9].sum()
    ap13=cut[rad<=13].sum()
    psf=ndimage.shift(kernel,[dy-cy,dx-cx],order=3,mode='constant')
    psf/=psf.sum()
    core=cut[cy-12:cy+13,cx-12:cx+13]
    fit=np.sum(core*psf)/np.sum(psf**2)
    noise=1.4826*np.median(np.abs(cut[(rad>16)&(rad<24)&~mask]))
    fit_err=np.sqrt(noise**2/np.sum(psf**2)+max(fit,0)/3.8)
    return dict(flux_adu=float(fit),stat_err_adu=float(fit_err),ap9_adu=float(ap9),ap13_adu=float(ap13),
                dq_in_12px=int(np.sum(mask&(rad<12))),center_shift_px=float(drift),peak_plus_sky_adu=float(peak),
                saturation_adu=float(h['SATURATE']))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--optical-dir',type=Path,required=True)
    ap.add_argument('--orbit-dir',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True)
    args=ap.parse_args()
    sys.path[:0]=[str(args.optical_dir),str(args.orbit_dir)]
    from decam_difference import sky_subtract, psf_stars
    import propagate
    resp=np.genfromtxt(ROOT/'bandpass.csv',delimiter=',',skip_header=1).T
    xp=Table.read(ROOT/'xp_spectra.vot.gz',format='votable')
    with gzip.open(ROOT/'gaia_field.csv.gz','rt') as f:catalog={r['source_id']:r for r in csv.DictReader(f)}
    results=[]
    fixed_fluxes=[80,150,250,400,700,1200]
    with fits.open(args.cache/'recent_n395.response') as fitsimg, fits.open(args.cache/'recent_n395_mask.response') as masks:
        hdr=fitsimg[0].header
        epoch=propagate.epoch(hdr['DATE-OBS'])+hdr['EXPTIME']/172800
        timespan=(epoch-propagate.epoch('2016-01-01'))/365.25
        states=propagate.family(json.loads((args.orbit_dir/'orbit_input.json').read_text()))
        ra,dec=propagate.sky(states,epoch)
        for ccd in [15,21]:
            im=fitsimg[ccd].data
            dq=masks[ccd].data
            h=fitsimg[ccd].header
            w=WCS(h)
            x0,y0=w.world_to_pixel_values(ra,dec)
            inside=(x0>=40)&(y0>=40)&(x0<im.shape[1]-40)&(y0<im.shape[0]-40)
            xs=np.arange(max(0,int(min(x0[inside])-140)),min(im.shape[1],int(max(x0[inside])+140)))
            ys=np.arange(max(0,int(min(y0[inside])-140)),min(im.shape[0],int(max(y0[inside])+140)))
            bad=dq[np.ix_(ys,xs)]!=0
            patch=np.nan_to_num(sky_subtract(np.where(bad,np.nan,im[np.ix_(ys,xs)])))
            stars,psf=psf_stars(patch,patch,bad)
            for spectrum in xp:
                sid=str(spectrum['datalinkID'])
                star=catalog[sid]
                de=float(star['dec'])
                a=float(star['ra'])+float(star['pmra'] or 0)*timespan/3.6e6/np.cos(np.deg2rad(de))
                de+=float(star['pmdec'] or 0)*timespan/3.6e6
                x,y=w.world_to_pixel_values(a,de);x=float(x);y=float(y)
                if not (30<x<im.shape[1]-31 and 30<y<im.shape[0]-31):continue
                flux=synthetic(*resp,spectrum['flux'])
                # Linear integration of spectral sigma is conservative under arbitrary positive wavelength correlations.
                err=synthetic(*resp,spectrum['flux_error'])
                mag=-2.5*np.log10(flux/3631.)
                tophat=synthetic(np.array([390.185,390.186,400.195,400.196]),np.array([0,1,1,0]),spectrum['flux'])
                m=photometry(im,dq,h,x,y,psf[0])
                flags=[]
                if star['phot_variable_flag']=='VARIABLE' or float(star['ruwe'])>1.4:flags.append('Gaia variability or RUWE')
                if err/flux>.22 or flux<=0:flags.append('XP SNR')
                if m['center_shift_px']>2:flags.append('centroid')
                if m['dq_in_12px']>0:flags.append('DQ')
                if m['peak_plus_sky_adu']>.80*m['saturation_adu']:flags.append('nonlinearity')
                if m['flux_adu']/m['stat_err_adu']<30:flags.append('image SNR')
                results.append(dict(ccd=ccd,id=sid,ra_deg=float(star['ra']),dec_deg=float(star['dec']),
                    x=x,y=y,gaia_G=float(star['phot_g_mean_mag']),bp_rp=float(star['bp_rp']),ruwe=float(star['ruwe']),
                    synthetic_Jy=flux,XP_err_upper_Jy=err,N395_AB=mag,
                    N395_tophat_minus_curve_mag=-2.5*np.log10(tophat/flux),**m,
                    zp_AB_per_total_ADU=mag+2.5*np.log10(m['flux_adu']),
                    selection='; '.join(flags) or 'accepted',kernel_stars=len(stars)))
    with (ROOT/'standards.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,results[0]);wr.writeheader();wr.writerows(results)
    summary={}
    for ccd in [15,21]:
        rows=[r for r in results if r['ccd']==ccd and r['selection']=='accepted']
        zp=np.array([r['zp_AB_per_total_ADU'] for r in rows])
        scatter=1.4826*np.median(abs(zp-np.median(zp)))
        s=np.array([np.hypot(1.08574*r['XP_err_upper_Jy']/r['synthetic_Jy'],
                             1.08574*r['stat_err_adu']/r['flux_adu']) for r in rows])
        median=float(np.median(zp))
        ap=np.array([r['N395_AB']+2.5*np.log10(r['ap13_adu']) for r in rows])
        apmad=1.4826*np.median(abs(ap-np.median(ap)))
        summary[str(ccd)]=dict(count=len(rows),matched_filter_zeropoint_AB_per_ADU=median,
               aperture13_zeropoint_AB_per_ADU=float(np.median(ap)),
               aperture13_rms_scatter=float(np.std(ap,ddof=1)),aperture13_mad_scatter=float(apmad),
               aperture13_formal_median_error_mag=float(1.253*max(apmad,np.median(s))/np.sqrt(len(ap))),
               rms_scatter=float(np.std(zp,ddof=1)),mad_scatter=float(scatter),
               median_single_star_upper_error_mag=float(np.median(s)),
               nominal_formal_median_error_mag=float(max(1.253*scatter/np.sqrt(len(zp)),1.253*np.median(s)/np.sqrt(len(zp)))),
               top_hat_mag_shift_median=float(np.median([r['N395_tophat_minus_curve_mag'] for r in rows])),
               AB_at_original_matched_filter_ADU={str(f):median-2.5*np.log10(f) for f in fixed_fluxes})
    (ROOT/'calibration.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':main()
