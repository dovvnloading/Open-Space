"""Heliocentric Ca II velocities and independent relative shifts on three public controls."""
import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
SPECTRA = {
    "D6-1": "D6-1-Gaia5805243926609660032-SALT-helio.dat",
    "D6-2": "D6-2-Gaia1798008584396457088-NOT-helio.dat",
    "D6-3": "D6-3-Gaia2156908318076164224-NOT-helio.dat",
}
BASE = "https://arxiv.org/src/1804.11163v3/anc/"
C = 299792.458
CA = np.array([3933.6614, 3968.4673])  # NIST Handbook, air Angstrom.


def doppler(v):
    return np.sqrt((1.0 + v/C)/(1.0 - v/C))


def calcium(wave, flux, instrument_fwhm, rng):
    keep = (wave >= 3912) & (wave <= 4024)
    wave, flux = wave[keep], flux[keep] / np.median(flux[keep])
    def model(params):
        velocity, width, kdepth, hdepth, level, slope = params
        k, h = CA*doppler(velocity)
        return (level + slope*(wave-3968)/60
                - kdepth*np.exp(-0.5*((wave-k)/width)**2)
                - hdepth*np.exp(-0.5*((wave-h)/width)**2))
    width_min = instrument_fwhm / 2.35482
    bounds = ([-1600, width_min, 0, 0, .1, -3], [2500, 20, 2, 2, 3, 3])
    def fit(y, start):
        return least_squares(lambda p: model(p) - y, start, bounds=bounds, max_nfev=600)
    starts = [[v, max(7, width_min), 0.3, 0.25, 1, 0] for v in range(-1500, 2001, 250)]
    best = min((fit(flux, s) for s in starts), key=lambda r: r.cost)
    # Wild residual bootstrap, one sign per resolution element; no pixel variance
    # provided. Re-use the global best basin within each bootstrap realization.
    block = ((wave-wave[0]) // instrument_fwhm).astype(int)
    pred = model(best.x)
    vel = []
    for _ in range(400):
        residual = (flux-pred)*rng.choice([-1.,1.], size=block[-1]+1)[block]
        vel.append(fit(pred+residual, best.x).x[0])
    return best.x, np.percentile(vel, [16,84])


def prepared(name, wave, flux):
    # Linear 1-A grid used only for continuum removal and template interpolation.
    grid = np.arange(3970., 6901.)
    sampled = np.interp(grid, wave, flux)
    if name == "D6-1":
        # Approximate constant resolution matching to NOT at 5400 A.
        sampled = gaussian_filter1d(sampled, np.sqrt(15**2-(5400/900)**2)/2.35482)
    return grid, sampled / gaussian_filter1d(sampled, 40.) - 1.


def relative_ccf(target, template):
    # Split the comparison across disjoint bands to expose template mismatch.
    grid, f = target
    tg, tf = template
    bands = [(4050,4370),(4370,4650),(4650,4950),(4950,5250),
             (5250,5550),(5600,5850),(5950,6200),(6400,6700)]
    velocity = np.arange(-2200., 2201., 10.)
    scores = []
    for lo, hi in bands:
        use = (grid>=lo) & (grid<hi)
        y = f[use]-f[use].mean()
        t = np.interp(grid[use][None,:] / doppler(velocity[:,None]), tg, tf)
        t -= t.mean(axis=1)[:,None]
        scores.append((t@y)/(np.linalg.norm(y)*np.linalg.norm(t,axis=1)))
    scores = np.asarray(scores)
    combined = scores.mean(axis=0)
    # Equal-weight Pearson r across intervals. Leave out one band at a time.
    leave = [velocity[np.argmax(np.delete(scores, i, axis=0).mean(axis=0))] for i in range(len(bands))]
    return velocity[np.argmax(combined)], np.max(combined), min(leave), max(leave)


def feature_indices(name, wave, flux, velocity):
    # Signed pseudo-equivalent width. Bounds are fixed in the stellar frame.
    rest = wave / doppler(velocity)
    features = [("MgII4481_region",4481.2,12), ("H_beta_region",4861.33,15),
                ("H_alpha_region",6562.81,15), ("HeI6678_region",6678.15,15),
                ("OI7774_triplet_region",7773.4,16)]
    rows = []
    for label, center, half in features:
        delta = rest - center
        if delta.min() > -(half+27) or delta.max() < half+27:
            continue
        side = (abs(delta)>half+6) & (abs(delta)<half+27)
        continuum = np.polyfit(delta[side], flux[side], 1)
        # Interpolate to exact edges to avoid pixel-phase truncation.
        use = abs(delta) < half
        pix = np.r_[-half, delta[use], half]
        norm = np.interp(pix,delta,flux) / np.polyval(continuum,pix)
        pew = np.trapezoid(1.-norm,pix)
        rows.append(dict(name=name,region=label,center_A=center,half_width_A=half,
                         pseudo_ew_A=round(pew,2),n_feature_pixels=int(use.sum())))
    return rows


def write_csv(filename, rows):
    with (HERE/filename).open("w", newline="") as out:
        d = csv.DictWriter(out, fieldnames=rows[0].keys(), lineterminator="\n")
        d.writeheader()
        d.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, help="Use previously fetched, receipt-matched input files.")
    args = parser.parse_args()
    receipts = []
    fits = []
    indices = []
    ready = {}
    prev_path = HERE / "receipts.json"
    previous = {r["name"]: r for r in json.loads(prev_path.read_text())} if prev_path.exists() else {}
    for i, (name, filename) in enumerate(SPECTRA.items()):
        if args.cache:
            content = (args.cache/filename).read_bytes()
        else:
            with urlopen(BASE+filename,timeout=45) as response:
                content = response.read()
                status = response.status
                stamp = datetime.now(timezone.utc).isoformat()
        digest = hashlib.sha256(content).hexdigest()
        if previous and digest != previous[name]["sha256"]:
            raise ValueError(f"Ancillary spectrum changed for {name}")
        assert b"air wavelengths corrected to heliocentric frame" in content[:300]
        x, y = np.loadtxt(io.BytesIO(content)).T
        if not args.cache:
            receipts.append(dict(name=name,url=BASE+filename,http_status=status,utc=stamp,
                                 bytes=len(content),sha256=digest,pixels=len(x)))
        beta, err = calcium(x, y, (3950./900 if i==0 else 15.), np.random.default_rng(2048+i))
        v, width, k, h, level, slope = beta
        fits.append(dict(name=name,ca_rv_kms=round(v,1),bootstrap_p16_kms=round(err[0],1),
                         bootstrap_p84_kms=round(err[1],1),gaussian_sigma_angstrom=round(width,2),
                         ca_k_depth=round(k,3), ca_h_depth=round(h,3)))
        ready[name] = prepared(name,x,y)
        indices.extend(feature_indices(name,x,y,v))
    if receipts:
        prev_path.write_text(json.dumps(receipts, indent=2)+"\n")
    for r in fits:
        if r["name"] == "D6-2":
            # A template's correlation with itself provides no cross-check.
            r.update(delta_v_vs_D62_kms="", ccf_mean_r="", band_jackknife_min="", band_jackknife_max="")
            continue
        dv, ccf, jlo, jhi = relative_ccf(ready[r["name"]],ready["D6-2"])
        r.update(delta_v_vs_D62_kms=int(dv), ccf_mean_r=round(ccf,3),
                 band_jackknife_min=int(jlo),band_jackknife_max=int(jhi))
    write_csv("line_measurements.csv",fits)
    write_csv("feature_indices.csv",indices)
    print(fits)
