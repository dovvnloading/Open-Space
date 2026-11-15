"""Fit primary and secondary Kepler difference images with the calibrated focal-plane PRF."""
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares
from scipy.spatial import Delaunay
from catalog_screen import HERE, NEARBY, fetch, read_vizier

TPF = "kplr008110757-2013098041711_lpd-targ.fits.gz"
PRF = "kplr17.3_2011265_prf.fits"
TPF_URL = "https://stpubdata.s3.amazonaws.com/kepler/public/target_pixel_files/0081/008110757/" + TPF
PRF_URL = "https://archive.stsci.edu/missions/kepler/fpc/prf/" + PRF
PERIOD, EPOCH = .518060209, 131.9425  # DR25, days and BKJD (BJD_TDB - 2454833)
EB = "2079577157039120128"
TARGET = "2079577161336953472"
NEIGHBOR = "2079577161336952960"


def event_differences(time, flux, center, halfwidth=.05, period=PERIOD, epoch=EPOCH):
    cycles = (time-epoch)/period-center
    orbit = np.floor(cycles+.5).astype(int)
    phase = cycles-orbit
    differences, dates, baselines = [], [], []
    for number in np.unique(orbit):
        same_event = orbit == number
        inside = same_event & (abs(phase) <= halfwidth)
        left = same_event & (phase >= -.26) & (phase <= -.14)
        right = same_event & (phase >= .14) & (phase <= .26)
        if min(inside.sum(), left.sum(), right.sum()) < 2:
            continue
        fraction = (time[inside].mean()-time[left].mean())/(time[right].mean()-time[left].mean())
        baseline = (1-fraction)*flux[left].mean(axis=0)+fraction*flux[right].mean(axis=0)
        differences.append(baseline-flux[inside].mean(axis=0))
        dates.append(time[inside].mean())
        baselines.append(baseline)
    return np.array(dates), np.array(differences), np.array(baselines)


def local_prf(path, detector_xy):
    with fits.open(path) as hdus:
        samples = np.array([[h.header['CRVAL1P'], h.header['CRVAL2P']] for h in hdus[1:]])
        triangles = Delaunay(samples)
        region = triangles.find_simplex(detector_xy)
        if region < 0:
            raise ValueError("Detector location outside PRF grid")
        bary = triangles.transform[region,:2]@(detector_xy-triangles.transform[region,2])
        kernel = sum(w*hdus[i+1].data for i,w in zip(triangles.simplices[region], np.r_[bary,1-bary.sum()]))
        step = hdus[1].header['CDELT1P']
        # Even-sized grid has the PRF origin between the central samples.
        grid = (np.arange(kernel.shape[0])+.5-kernel.shape[0]/2)*step
    return RectBivariateSpline(grid,grid,kernel)


def fit_image(values, errors, x, y, kernel, initial):
    def residual(params):
        x0,y0,amplitude,sky = params
        return (amplitude*kernel.ev(y-y0,x-x0)+sky-values)/errors
    solution = least_squares(residual, [*initial, max(values)*2, 0.],
                             bounds=([-1,-1,0,-np.inf], [6,6,np.inf,np.inf]))
    return solution.x, float(solution.fun@solution.fun)


def fixed_location(values, errors, x, y, kernel, at):
    template = kernel.ev(y-at[1],x-at[0])
    design = np.column_stack((template,np.ones(len(x)))) / errors[:,None]
    amplitude, sky = np.linalg.lstsq(design,values/errors,rcond=None)[0]
    if amplitude < 0:
        amplitude, sky = 0., np.average(values,weights=1/errors**2)
    residual = (amplitude*template+sky-values)/errors
    return float(residual@residual)


def propagated_field(wcs,epoch):
    field = read_vizier(fetch("nearby.tsv", NEARBY)).rename(columns={"Source":"gaia_id"})
    for col in field.columns[1:]:
        field[col] = pd.to_numeric(field[col], errors="raise")
    dt = epoch-2016.
    ra = field.RA_ICRS + dt*field.pmRA/(3.6e6*np.cos(np.deg2rad(field.DE_ICRS)))
    dec = field.DE_ICRS + dt*field.pmDE/3.6e6
    field['x_pix'],field['y_pix'] = wcs.world_to_pixel_values(ra,dec)
    return field


def main():
    source = fetch(TPF, TPF_URL, "b5d68ed38e0a602e2ee1baca26946582c1248e4fcbab166f125e1be5389f17f9")
    model = fetch(PRF, PRF_URL, "7ab5c948bc508b9af855f9cb5d8cef39323e997aa38b541aa15c7cb2c6affac8")
    with fits.open(source) as hdus:
        table = hdus[1].data
        good = (table['QUALITY']==0) & np.isfinite(table['TIME'])
        time, stamp = table['TIME'][good], table['FLUX'][good].astype(float)
        aperture = hdus[2].data
        header = hdus[2].header
        physical = WCS(header,key='P')
        wcs = WCS(header)
    finite = np.all(np.isfinite(stamp),axis=0)
    flux = stamp[:,finite]
    optimal = ((aperture&2)>0)[finite]
    yy,xx = np.indices(aperture.shape)
    x,y = xx[finite], yy[finite]
    # PRF calibration changes over ~1000 pixels. Evaluate at stamp center.
    kernel = local_prf(model, np.asarray(physical.pixel_to_world_values(x.mean(),y.mean())))
    epoch = Time(2454833+np.mean(time),format='jd',scale='tdb').jyear
    field = propagated_field(wcs,epoch)
    field.to_csv(HERE/"gaia_field.csv",index=False)
    coords = field.set_index('gaia_id')[['x_pix','y_pix']]
    scale = 3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
    rng = np.random.default_rng(7865)
    results, maps = [], {}
    for label, center in [('primary',0.),('secondary',.5)]:
        dates,differences,baseline = event_differences(time,flux,center)
        for name,use in [('all',np.ones(len(dates),bool)),('first',dates<np.median(time)),('second',dates>=np.median(time))]:
            arr = differences[use]
            values = arr.mean(axis=0)
            errors = arr.std(axis=0,ddof=1)/np.sqrt(len(arr))
            pars,chi2 = fit_image(values,errors,x,y,kernel,[x[np.argmax(values)],y[np.argmax(values)]])
            if name=='all':
                maps[label]=values
            # Resample contiguous groups of eight events to retain short-term correlations.
            grouped = np.arange(len(arr))
            boot = []
            for _ in range(350):
                starts = rng.integers(len(arr),size=(len(arr)+7)//8)
                ix = ((starts[:,None]+np.arange(8))%len(arr)).ravel()[:len(arr)]
                pos,_ = fit_image(arr[ix].mean(axis=0),errors,x,y,kernel,pars[:2])
                boot.append(pos[:2])
            low,high = np.quantile(boot,[.025,.975],axis=0)
            raw_xy = np.array([values@x,values@y])/values.sum()
            row = dict(eclipse=label,subset=name,n_events=len(arr),prf_x=pars[0],prf_y=pars[1],
                       x_lower=low[0],x_upper=high[0],y_lower=low[1],y_upper=high[1],
                       moment_x=raw_xy[0],moment_y=raw_xy[1],chi2_free=chi2,
                       stamp_loss_e_s=values.sum(),ap_loss_e_s=values[optimal].sum(),
                       ap_depth_ppm=1e6*values[optimal].sum()/baseline[use][:,optimal].sum(axis=1).mean())
            for key,sid in [('target',TARGET),('eb',EB),('neighbor',NEIGHBOR)]:
                at=coords.loc[sid].to_numpy()
                row[f'{key}_sep_arcsec'] = np.linalg.norm(at-pars[:2])*scale
                row[f'{key}_fixed_chi2'] = fixed_location(values,errors,x,y,kernel,at)
            results.append(row)
    pd.DataFrame(results).to_csv(HERE/"eclipse_localizations.csv",index=False)
    pd.DataFrame(dict(x=x,y=y,optimal_aperture=optimal,median_flux_e_s=np.median(flux,axis=0),
                      primary_loss_e_s=maps['primary'],secondary_loss_e_s=maps['secondary'])).to_csv(HERE/"pixel_maps.csv",index=False)
    info = dict(cadences_total=len(table),cadences_kept=len(time),first_bkjd=float(min(time)),last_bkjd=float(max(time)),epoch_jyear=epoch,usable_pixels=len(x),plate_scale_arcsec_pix=scale,
                primary_secondary_map_cosine=float(np.dot(maps['primary'],maps['secondary'])/np.linalg.norm(maps['primary'])/np.linalg.norm(maps['secondary'])))
    print(json.dumps(info,indent=2))
    pd.DataFrame([info]).to_csv(HERE/'stamp_info.csv',index=False)


if __name__ == "__main__":
    main()
