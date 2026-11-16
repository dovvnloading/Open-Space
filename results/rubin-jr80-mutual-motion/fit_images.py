"""Fit Rubin's measured PSF to individual, WCS-calibrated alert cutouts."""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates, spline_filter
from scipy.optimize import least_squares

from cutouts import ROOT, read_cutout, rows

def fit(job):
    mid, kind, label, radius, blur = job
    (header, image), (unc, variance), (_, psf) = read_cutout(mid, kind)
    if unc.get("UTYPE") != "VarianceUncertainty":
        raise ValueError("Expected a variance plane")
    yy, xx = np.indices(image.shape)
    xx = xx - (float(header["CRPIX1"]) - 1)
    yy = yy - (float(header["CRPIX2"]) - 1)
    use = (np.hypot(xx, yy) < radius) & np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    xx, yy, image, sigma = xx[use], yy[use], image[use], np.sqrt(variance[use])
    if blur:
        psf = gaussian_filter(psf, blur)
    py, px = (np.array(psf.shape)-1)/2
    spline = spline_filter(psf, order=3)
    def profile(x, y):
        return map_coordinates(spline, [yy-y+py, xx-x+px], order=3, mode="constant", prefilter=False)
    background = np.median(image[np.hypot(xx, yy) > radius-2])
    flux_scale = max(np.sum((image - background)[np.hypot(xx,yy) < 5]), 1)
    def residual_one(params):
        flux, x, y, sky = params
        return (image - sky - flux_scale * flux * profile(x,y)) / sigma
    def residual_pair(params):
        flux, x, y, sky, dx, dy, weight = params
        # Place both PSFs relative to their common flux-weighted centroid.
        first = profile(x - weight*dx, y - weight*dy)
        second = profile(x + (1-weight)*dx, y + (1-weight)*dy)
        return (image - sky - flux_scale*flux*((1-weight)*first + weight*second))/sigma
    lower, upper = [0,-4,-4,-np.inf], [10,4,4,np.inf]
    one = least_squares(residual_one, [1,0,0,background], bounds=(lower,upper), max_nfev=100)
    starts = [(0,0)] + [(2.5*np.sin(a),2.5*np.cos(a)) for a in np.deg2rad([0,45,90,135])]
    best = None
    for dx, dy in starts:
        result = least_squares(residual_pair, [*one.x, dx, dy, .48],
                               bounds=(lower+[-6,-6,.1], upper+[6,6,.9]), max_nfev=100)
        if best is None or result.cost < best.cost:
            best = result
    if not best.success:
        best = least_squares(residual_pair, best.x, bounds=(lower+[-6,-6,.1], upper+[6,6,.9]),
                             max_nfev=1000, ftol=1e-7, xtol=1e-7)
    transform = np.array([[float(header[f"PC{i}_{j}"])*float(header[f"CDELT{i}"])
                           for j in (1,2)] for i in (1,2)]) * 3600
    # TAN local WCS coordinates are the east/north tangent plane, in arcseconds.
    en = transform @ best.x[4:6]
    if en[0] < 0:
        en = -en
    sep = np.linalg.norm(en)
    angle = np.degrees(np.arctan2(en[0], en[1])) % 180
    # Local errors conditional on the two-PSF model. Account for excess residual
    # variance; spatial correlations and seeing model errors are assessed separately.
    dof = len(image)-len(best.x)
    cov = np.linalg.pinv(best.jac.T @ best.jac) * max(1.,2*best.cost/dof)
    sky_cov = transform @ cov[4:6,4:6] @ transform.T
    g_pa = np.array([en[1],-en[0]]) / max(sep**2,1e-15)
    g_sep = en / max(sep,1e-15)
    return dict(id=mid, kind=kind, variant=label, converged=bool(best.success and one.success), n_pix=len(image), chi2_single=2*one.cost,
                chi2_pair=2*best.cost, delta_chi2=2*(one.cost-best.cost),
                east_arcsec=en[0], north_arcsec=en[1], separation_arcsec=sep, pa_deg=angle,
                sigma_east=np.sqrt(sky_cov[0,0]), sigma_north=np.sqrt(sky_cov[1,1]),
                cov_east_north=sky_cov[0,1], sigma_sep=np.sqrt(max(0,g_sep @ sky_cov @ g_sep)),
                sigma_pa_deg=np.degrees(np.sqrt(max(0,g_pa @ sky_cov @ g_pa))),
                faint_flux_fraction=min(best.x[6],1-best.x[6]), total_flux_nJy=flux_scale*best.x[0])

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    targets = {r["measurement_id"] for r in rows("target.csv")}
    jobs = []
    for r in rows("stamps.csv"):
        if r["kind"] == "Template":
            continue
        jobs.append((r["id"], r["kind"], "nominal", 8.5, 0.))
        if r["id"] in targets and r["kind"] == "Science":
            jobs.extend([(r["id"], "Science", "inner", 6., 0.),
                         (r["id"], "Science", "broadened", 8.5, .40)])
    with ProcessPoolExecutor(max_workers=args.workers) as engine:
        fits = list(engine.map(fit, jobs))
    with (ROOT / "output" / "fits.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fits[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(fits)
    print(f"Wrote {len(fits)} fits.")

if __name__ == "__main__":
    main()
