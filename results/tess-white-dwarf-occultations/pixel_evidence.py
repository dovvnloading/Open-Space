"""Measure positive deficits and photocenter displacement on FITS target pixels.

Input times and ephemerides use BJD_TDB - 2457000 in days. For an isolated
event supply period_days=None. The difference fit is descriptive and uses a
single circular Gaussian with free amplitude, position, width and sky.
"""
import numpy as np
from scipy.optimize import least_squares


def profile_center(image):
    y, x = np.indices(image.shape)
    mask = np.isfinite(image)
    if mask.sum() < 9:
        return None
    x, y, z = x[mask], y[mask], image[mask]
    amplitude = np.max(z) - np.median(z)
    if amplitude <= 0:
        return None
    k = np.argmax(z)
    def resid(p):
        u,v,s,a,c = p
        return c+a*np.exp(-((x-u)**2+(y-v)**2)/(2*s*s))-z
    fit = least_squares(resid, [x[k],y[k],.9,amplitude,np.median(z)],
                        bounds=([-.5,-.5,.35,0,-np.inf],
                                [image.shape[1]-.5,image.shape[0]-.5,2.5,np.inf,np.inf]))
    return fit.x[:2]


def measure_pixels(time, cube, quality, target_xy, epoch_btjd, duration_min, period_days=None):
    t = np.asarray(time,float)
    q = np.asarray(quality)
    ix = (q==0)&np.isfinite(t)
    t, cube = t[ix], np.asarray(cube)[ix].astype(float)
    width = duration_min/1440
    if period_days is None:
        centers = [epoch_btjd]
    else:
        k = np.arange(np.ceil((t.min()-epoch_btjd)/period_days),
                      np.floor((t.max()-epoch_btjd)/period_days)+1)
        centers = epoch_btjd+k*period_days
    y,x = np.indices(cube.shape[1:])
    aperture = np.hypot(x-target_xy[0],y-target_xy[1])<=1.5
    edge = (x==0)|(y==0)|(x==x.max())|(y==y.max())
    stamps, shifts = [], []
    flux_diffs=[]
    for c in centers:
        phase = (t-c)/width
        central = abs(phase) <= .5
        before = (phase >= -2.5)&(phase <= -1)
        after = (phase >= 1)&(phase <= 2.5)
        if min(central.sum(),before.sum(),after.sum()) < 1:
            continue
        ox = .5*(np.nanmean(cube[before],axis=0)+np.nanmean(cube[after],axis=0))
        ix = np.nanmean(cube[central],axis=0)
        # Remove local background from each stamp, independently.
        ox -= np.nanmedian(ox[edge])
        ix -= np.nanmedian(ix[edge])
        delta = ox-ix
        baseline = np.nansum(ox[aperture])
        if not np.isfinite(baseline) or baseline<=0:
            continue
        stamps.append(delta)
        flux_diffs.append(np.nansum(delta[aperture])/baseline)
        center = lambda a: np.array([np.nansum(a[aperture]*v[aperture]) for v in (x,y)])/np.nansum(a[aperture])
        if np.nansum(ix[aperture])>0:
            shifts.append(center(ix)-center(ox))
    if not stamps:
        return {'status':'no supported events'}
    stamps = np.asarray(stamps)
    image = np.nanmean(stamps,axis=0)
    peak = profile_center(image)
    result = dict(events=len(stamps),target_xy=list(map(float,target_xy)),median_aperture_depth=float(np.median(flux_diffs)),
                  photocenter_shift_xy=np.nanmedian(shifts,axis=0).tolist() if shifts else None,
                  difference_image_e_s=image.tolist())
    if peak is None:
        return dict(result,status='no fit')
    dist = np.linalg.norm(peak-target_xy)
    result.update(status='measured',difference_xy=peak.tolist(),offset_pixel=float(dist))
    if len(stamps)>=3:
        rng=np.random.default_rng(79)
        samples=[profile_center(np.nanmean(stamps[rng.integers(len(stamps),size=len(stamps))],axis=0))
                 for _ in range(100)]
        boot=np.array([s for s in samples if s is not None])
        if len(boot)>1:
            result['bootstrap_xy_std_pixel']=boot.std(axis=0,ddof=1).tolist()
    return result
