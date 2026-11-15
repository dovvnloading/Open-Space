"""Test the 15-hour component across sectors, products, and same-detector controls."""
import argparse
import csv
from pathlib import Path

from astropy.io import fits
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, lombscargle

from measure import MODES, ROOT, coefficients, design, verified_product

FREQUENCY = 1.57064


def binned(data, key):
    time = data["TIME"]
    flux = data[key]
    keep = (data["QUALITY"] == 0) & np.isfinite(time) & np.isfinite(flux)
    keep &= np.isfinite(data["PDCSAP_FLUX"])
    time = time[keep]
    flux = flux[keep].astype(float)
    if key.endswith("FLUX"):
        flux = 1e6 * (flux / np.median(flux) - 1)
    else:
        flux -= np.median(flux)
    segs = np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > .05) + 1)
    use = []
    for ids in segs:
        if len(ids) <= 600:
            continue
        flux[ids] -= gaussian_filter1d(flux[ids], 360, mode="nearest")
        use.extend(ids[300:-300])
    time, flux = time[use], flux[use]
    # Remove the rapid stellar pulsations before non-overlapping 20-minute bins.
    if key.endswith("FLUX"):
        modes = design(time - np.mean(time), MODES)
        flux -= modes @ coefficients(modes, flux)
    ix = np.floor((time - time[0]) * 72).astype(int)
    counts = np.bincount(ix)
    ok = counts >= 6
    return (np.bincount(ix, weights=time)[ok] / counts[ok],
            np.bincount(ix, weights=flux)[ok] / counts[ok])


def amplitude(time, flux, nu):
    # Absolute phase epoch is BTJD 3000 in all sectors.
    x = design(time - 3000, np.atleast_1d(nu))
    b = coefficients(x, flux - np.mean(flux))
    return float(np.hypot(*b)), float(np.arctan2(-b[1], b[0]))


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: f"{v:.7g}" if isinstance(v, (np.floating,float)) else v
                        for k,v in r.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    rows, bands, by_sector = [], [], {}
    for name in ("manifest.csv", "control_manifest.csv"):
        with (ROOT/name).open() as handle:
            entries = list(csv.DictReader(handle))
        for r in entries:
            control = name == "control_manifest.csv"
            data = fits.getdata(verified_product(r, args.cache), 1)
            columns = (["PDCSAP_FLUX"] if control else
                       ["PDCSAP_FLUX","SAP_FLUX","SAP_BKG","POS_CORR1","POS_CORR2"])
            for key in columns:
                if not np.any(np.isfinite(data[key])):
                    continue
                t,y=binned(data,key)
                if not len(t):
                    continue
                amp,phase=amplitude(t,y,FREQUENCY)
                side = np.linspace(.85,2.3,147)
                side = side[abs(side-FREQUENCY)>.12]
                local=np.array([amplitude(t,y,f)[0] for f in side])
                rows.append(dict(sector=int(r["sector"]),tic=int(r["product"].split("-")[2]),
                                 column=key,bins=len(t),amplitude=amp,phase_rad=phase,
                                 sideband_median=np.median(local),ratio=amp/np.median(local)))
                if not control and key=="PDCSAP_FLUX":
                    by_sector[int(r["sector"])]=t,y
    groups=[[4,5,6,7],[32,33,34],[87],[95],[98],[87,95,98],list(by_sector)]
    for group in groups:
        t=np.concatenate([by_sector[s][0] for s in group])
        y=np.concatenate([by_sector[s][1] for s in group])
        grid=np.arange(1.48,1.66,min(.0002,1/(8*np.ptp(t))))
        power=lombscargle(t-3000,y,2*np.pi*grid)
        peaks=find_peaks(power)[0]
        for rank,pk in enumerate(sorted(peaks,key=lambda k:power[k],reverse=True)[:4],1):
            amp,phase=amplitude(t,y,grid[pk])
            bands.append(dict(sectors=" ".join(map(str,group)),rank=rank,
                              frequency_d=grid[pk],amplitude_ppm=amp,phase_rad=phase))
    write_csv(args.out/"slow_checks.csv",rows)
    write_csv(args.out/"frequency_aliases.csv",bands)


if __name__=="__main__":
    main()
