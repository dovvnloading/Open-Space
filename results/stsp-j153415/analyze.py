"""Replay individual DESI line velocities, ZTF localization, and Gaia hierarchy."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
C = 299792.458
IDS = ['1375395351637562240', '1375395141182830080',
       '1375395317277825536', '1375395317276861696']
VISITS = (127032, 131968)
# Vacuum wavelengths. Na absorption references are the air lines near 8183.26 and
# 8194.82 Angstrom converted to vacuum. The red feature is a multiplet blend.
LINES = {'Halpha': ('r', [6564.61], 1), 'Hbeta': ('b', [4862.68], 1),
         'Hgamma': ('b', [4341.68], 1), 'NaI': ('z', [8185.51, 8197.08], -1)}
PERIOD = 0.252027
T0 = 2458733.686


def rows(name):
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))


def separation(a, b):
    ra, dec = np.radians([float(a['ra']), float(a['dec'])])
    rb, db = np.radians([float(b['ra']), float(b['dec'])])
    sin2 = np.sin((db-dec)/2)**2 + np.cos(dec)*np.cos(db)*np.sin((ra-rb)/2)**2
    return float(np.degrees(2*np.arcsin(np.sqrt(sin2)))*3600)


def fit_line(spectra, exposure, line, window, quadratic=False):
    band, refs, sign = LINES[line]
    key = f'{exposure}_{band}_'
    wave = spectra[key+'wave']
    ivar = spectra[key+'ivar']
    keep = ((wave > refs[0]-window) & (wave < refs[-1]+window)
            & (spectra[key+'mask'] == 0) & (ivar > 0))
    w = wave[keep]
    f = spectra[key+'flux'][keep]
    wt = np.sqrt(ivar[keep])
    x = (w-np.mean(refs))/window
    # Profile a locally sloping continuum and one Gaussian per transition. Both
    # doublet components share a velocity and width but have independent depths.
    k = 2 + int(quadratic)
    def model(a):
        v, width = a[k:k+2]
        y = a[0] + a[1]*x
        if quadratic:
            y = y + a[2]*x*x
        for ref, amplitude in zip(refs, a[k+2:]):
            y = y + sign*amplitude*np.exp(-0.5*((w-ref*(1+v/C))/width)**2)
        return y
    initial = [np.median(f), 0.] + ([0.] if quadratic else [])
    lo = [-np.inf]*k + [-500, 0.55] + [0]*len(refs)
    hi = [np.inf]*k + [500, 4.5] + [np.inf]*len(refs)
    def solve(v):
        x0 = initial + [v, 1.5] + [np.ptp(f)*0.6]*len(refs)
        return least_squares(lambda a: wt*(model(a)-f), x0, bounds=(lo,hi))
    fit = min((solve(v) for v in (-240, 0, 180)), key=lambda fit: fit.cost)
    cov = np.linalg.pinv(fit.jac.T @ fit.jac)
    redchi = 2*fit.cost/(len(f)-len(fit.x))
    err = np.sqrt(cov[k,k]*max(1,redchi))
    return {'velocity':float(fit.x[k]), 'error':float(err), 'width_A':float(fit.x[k+1]),
            'reduced_chi2':float(redchi), 'pixels':int(len(f))}


def spectrum_controls(spec, exposure):
    result = {}
    for name,(band,refs,_) in LINES.items():
        if name == 'NaI':
            continue
        k = f'{exposure}_{band}_'
        w = spec[k+'wave']
        good = (spec[k+'mask']==0) & (spec[k+'ivar']>0)
        ref = refs[0]
        line = abs(w - ref) < 8
        side = (abs(w-ref)>10) & (abs(w-ref)<19)
        height = np.max(spec[k+'flux'][line & good])-np.median(spec[k+'flux'][side & good])
        sw = spec[k+'sky_wave']
        sf = np.where(spec[k+'sky_mask']==0, spec[k+'sky_flux'], np.nan)
        sl = abs(sw - ref) < 8
        ss = (abs(sw-ref)>10) & (abs(sw-ref)<19)
        valid = np.any(np.isfinite(sf[:,sl]),axis=1) & np.any(np.isfinite(sf[:,ss]),axis=1)
        peaks = np.nanmax(sf[valid][:,sl],axis=1) - np.nanmedian(sf[valid][:,ss],axis=1)
        result[name] = {'stellar_peak_above_local':float(height), 'usable_sky':len(peaks),
                        'sky_peak_p95':float(np.percentile(peaks,95)),
                        'sky_at_least_stellar':int(np.sum(peaks>=height))}
    return result


def photometry():
    import gzip
    with gzip.open(DATA/'ztf.csv.gz','rt') as f:
        data = list(csv.DictReader(f))
    out = {}
    for sid in sorted({r['source_id'] for r in data}):
        by_band = {}
        for band in ('zg','zr','zi'):
            sample = [r for r in data if r['source_id']==sid and r['filtercode']==band
                      and int(r['catflags'])==0 and float(r['magerr'])<0.2
                      and float(r['mag'])<float(r['limitmag'])-0.3]
            t = np.array([float(r['hjd']) for r in sample])
            m = np.array([float(r['mag']) for r in sample])
            sigma = np.hypot([float(r['magerr']) for r in sample], .02)
            oids = np.array([r['oid'] for r in sample])
            phase = 2*np.pi*(t-T0)/PERIOD
            X = np.column_stack([*(oids==oid for oid in np.unique(oids)), np.sin(phase), np.cos(phase)])
            a = np.linalg.lstsq(X/sigma[:,None], m/sigma, rcond=None)[0]
            by_band[band] = {'n':len(m),'semiamplitude_mag':float(np.hypot(*a[-2:])),
                             'phase_faintest':float(np.arctan2(a[-2],a[-1])/(2*np.pi)%1)}
        out[sid] = by_band
    return out


def main():
    stars = {r['source_id']:r for r in rows('quad.csv')}
    geom = {}
    for i,j in ((0,1),(2,3),(0,2)):
        a,b = stars[IDS[i]], stars[IDS[j]]
        sep = separation(a,b)
        geom[f'{i+1}-{j+1}'] = {'arcsec':sep, 'au_at_variable_parallax':sep*1000/float(stars[IDS[0]]['parallax'])}
    velocities, differential = {}, {}
    controls = {}
    with np.load(DATA/'spectra.npz') as spectra:
        for e in VISITS:
            velocities[e] = {line:fit_line(spectra,e,line,17.) for line in LINES}
            controls[e] = spectrum_controls(spectra,e)
        for line in LINES:
            a,b = [velocities[e][line] for e in VISITS]
            trials = [fit_line(spectra,VISITS[1],line,w,q)['velocity'] -
                      fit_line(spectra,VISITS[0],line,w,q)['velocity']
                      for w in (12.,17.,22.) for q in (False,True)]
            differential[line] = dict(delta=b['velocity']-a['velocity'], formal_error=float(np.hypot(a['error'],b['error'])),
                                       min_tested=min(trials),max_tested=max(trials))
    times = [float(v['hjd_utc_mid']) for v in rows('visits.csv')]
    phases = [(t-T0)/PERIOD % 1 for t in times]
    dv = abs(differential['NaI']['delta'])
    spacing = (times[1]-times[0])/PERIOD % 1
    kmin = dv/(2*abs(np.sin(np.pi*spacing)))
    orbit = {'nominal_range_half_km_s':dv/2, 'VSX_phases':phases,
             'spacing_cycles_modulo_one':spacing, 'circular_K_min_km_s':float(kmin),
             'circular_min_mass_function_Msun':float(1.036149e-7*PERIOD*kmin**3)}
    return {'conditional_sodium_orbit':orbit, 'geometry':geom,'line_fits_km_s':velocities,'april_minus_march_km_s':differential,
            'sky_controls':controls,'fixed_period_ZTF':photometry()}


if __name__=='__main__':
    print(json.dumps(main(),indent=2))
