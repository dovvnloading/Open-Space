"""Verify the SDSS and DESI extracted line without downloading the original CCDs."""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent


def measure(label, wave, flux, error):
    # Spectrograph fibers see different areas, so each fit gets its own local continuum.
    mask = np.isfinite(error) & (error > 0) & np.isfinite(flux) & (np.abs(wave-6446.1)<22)
    x, y, e = wave[mask]-6446.1, flux[mask], error[mask]
    def profile(x, level, slope, area, center, width):
        return level+slope*x+area*np.exp(-.5*((x-center)/width)**2)/(np.sqrt(2*np.pi)*width)
    if label == "SDSS":
        params,cov=curve_fit(profile,x,y,p0=[np.median(y),0,65,0,1.8],sigma=e,
                            absolute_sigma=True, bounds=([-np.inf,-np.inf,0,-3,.35],
                                                          [np.inf,np.inf,np.inf,3,6]))
        level,slope,area,offset,width=params
        return float(area),float(np.sqrt(cov[2,2])),float(6446.1+offset),float(width)
    # Match the center and width from the previous epoch to avoid fitting noise for a new centroid.
    center,width=label
    g=np.exp(-.5*((x-(center-6446.1))/width)**2)/(np.sqrt(2*np.pi)*width)
    M=np.column_stack([np.ones_like(x),x,g])/e[:,None]
    cov=np.linalg.inv(M.T@M)
    coef=cov@(M.T@(y/e))
    return float(coef[2]),float(np.sqrt(cov[2,2]))


if __name__ == '__main__':
    with (ROOT/'line-pixels.csv').open() as handle:
        rows=list(csv.DictReader(handle))
    def data(epoch):
        return np.array([[float(row[k]) for k in ['wave_vac_A','flux','error']]
                         for row in rows if row['epoch']==epoch]).T
    sdss=measure('SDSS',*data('SDSS'))
    desi=measure(sdss[2:],*data('DESI'))
    print(f"SDSS integrated flux {sdss[0]:.2f} +/- {sdss[1]:.2f}, center {sdss[2]:.3f} A, sigma {sdss[3]:.3f} A")
    print(f"DESI at fixed center and width {desi[0]:.2f} +/- {desi[1]:.2f}")
