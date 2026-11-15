"""Reproduce the independent eclipse timing, odd/even, and source checks.

Run python measure.py CACHE_DIRECTORY. Requires numpy, scipy, astropy, matplotlib.
Downloads three 7x7 MAST TESSCut archives and the IRSA ZTF source-cone CSV.
"""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
RA, DEC = 150.4282391, 51.09265768
SOURCE = "827111612020422528"
SECTORS = (21, 48, 75)
SPACING = 1.583154209  # Initial blind ZTF eclipse spacing, days
TREF = 2459000 + SPACING * 0.74  # Initial ZTF conjunction


def read_remote(path, url, compress=False):
    if not path.exists():
        with urlopen(url, timeout=120) as result:
            raw = result.read()
        path.write_bytes(gzip.compress(raw) if compress else raw)
    return path


def ztf(cache):
    args = {"POS": f"CIRCLE {RA:.8f} {DEC:.8f} 0.0010", "FORMAT": "csv"}
    name = read_remote(cache / "ztf.csv.gz", "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(args), True)
    with gzip.open(name, "rt") as stream:
        records = list(csv.DictReader(stream))
    unique = {}
    for oid in sorted({r["oid"] for r in records}):
        rows = [r for r in records if r["oid"] == oid and int(r["catflags"]) == 0 and 0.001 < float(r["magerr"]) < 0.12]
        if len(rows) < 12:
            continue
        median = np.median([float(r["mag"]) for r in rows])
        for r in rows:
            key = (r["expid"], r["filtercode"])
            if key not in unique or float(r["magerr"]) < float(unique[key]["magerr"]):
                unique[key] = dict(r, dm=float(r["mag"])-median)
    return sorted(unique.values(), key=lambda r: float(r["hjd"]))


def cutout(cache, sector):
    args = {"ra":RA, "dec":DEC, "x":7, "y":7, "units":"px", "sector":sector}
    path = read_remote(cache / f"s{sector}.zip", "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + urlencode(args))
    with ZipFile(path) as archive:
        member = next(n for n in archive.namelist() if n.endswith(".fits"))
        with fits.open(io.BytesIO(archive.read(member))) as hdus:
            data = hdus[1].data
            time = data["TIME"].copy() + 2457000
            pixels = data["FLUX"].astype(float)
            good = (data["QUALITY"] == 0) & np.isfinite(time) & np.isfinite(pixels).all(axis=(1,2))
            time, pixels = time[good], pixels[good]
            xy = WCS(hdus[2].header).world_to_pixel_values(RA, DEC)
    border = np.ones((7,7), bool)
    border[1:6,1:6] = False
    pixels -= np.median(pixels[:,border],axis=1)[:,None,None]
    aperture = pixels[:,2:5,2:5].sum(axis=(1,2))
    interval = np.median(np.diff(time))
    relative = aperture / median_filter(aperture, size=int(0.4/interval)|1, mode="nearest")
    return time, relative, pixels, xy, interval


def eclipse_times(time, relative, sector, interval):
    cycles = np.rint((time-TREF)/SPACING).astype(int)
    offsets = time-TREF-cycles*SPACING
    exposure_nodes = ((np.arange(15)+0.5)/15-0.5)*interval
    fits_out = []
    for cycle in np.unique(cycles):
        use = (cycles == cycle) & (abs(offsets) < 0.09)
        t, y = offsets[use], relative[use]
        if len(t)<5 or min(sum(t<0.004),sum(t>0.004))<2:
            continue

        def residual(pars):
            center, width, depth, base, power = pars
            phase = (t[:,None]-center+exposure_nodes)/width
            profile = np.maximum(0,1-phase**2)**power
            return base-depth*profile.mean(axis=1)-y

        fit = least_squares(residual, [0.005,0.045,0.4,1,1.5],
                            bounds=([-0.025,0.026,0.2,0.975,0.5],[0.027,0.065,0.55,1.025,3]),
                            loss="soft_l1", f_scale=0.009, ftol=1e-12, xtol=1e-12, gtol=1e-12)
        normal = fit.jac.T @ fit.jac
        if np.linalg.cond(normal)>1e11:
            continue
        error = np.sqrt(np.linalg.inv(normal)[0,0])*max(0.003,np.sqrt(np.mean(fit.fun**2)))
        center, width, depth, _, _ = fit.x
        fits_out.append(dict(sector=sector, cycle=int(cycle), bjd_tdb=TREF+cycle*SPACING+center,
                             formal_error_s=error*86400, depth=depth, duration_hr=width*48,
                             cadences=len(t)))
    return fits_out


def write_rows(name, rows):
    with (ROOT/name).open("w",newline="") as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(cache):
    cache.mkdir(parents=True,exist_ok=True)
    ground = ztf(cache)
    results, times, curves = [], [], []
    for s in SECTORS:
        t, y, pixels, xy, interval = cutout(cache,s)
        fits_s = eclipse_times(t,y,s,interval)
        times.extend(fits_s)
        cycles=np.array([r["cycle"] for r in fits_s])
        centers=np.array([r["bjd_tdb"] for r in fits_s])
        depths=np.array([r["depth"] for r in fits_s])
        design=np.column_stack([np.ones(len(cycles)),cycles-cycles.mean(),cycles%2])
        shift=(centers-TREF-cycles*SPACING)*1440
        coeff=np.linalg.lstsq(design,shift,rcond=None)[0]
        residual=shift-design@coeff
        # Difference centroid on the central 5x5 pixels, positive difference flux only.
        offset=(t-TREF+SPACING/2)%SPACING-SPACING/2
        difference=np.median(pixels[(abs(offset)>.09)&(abs(offset)<.20)],axis=0)-np.median(pixels[abs(offset-.005)<.012],axis=0)
        weight=np.maximum(difference[1:6,1:6],0)
        gy,gx=np.mgrid[1:6,1:6]
        cx,cy=(np.sum(pos*weight)/weight.sum() for pos in (gx,gy))
        results.append(dict(sector=s,cadences=len(t),fitted_events=len(fits_s),cadence_min=round(interval*1440,3),
                            even_depth=float(np.median(depths[cycles%2==0])),odd_depth=float(np.median(depths[cycles%2==1])),
                            median_duration_hr=float(np.median([r["duration_hr"] for r in fits_s])),
                            odd_minus_even_min=float(coeff[2]),timing_scatter_min=float(np.std(residual)),
                            difference_x=cx,difference_y=cy,catalog_x=float(xy[0]),catalog_y=float(xy[1])))
        curves.append((t,y))
    cycle=np.array([r["cycle"] for r in times])
    center=np.array([r["bjd_tdb"] for r in times])
    design=np.column_stack([np.ones(len(cycle)),cycle,cycle%2])
    pars=np.linalg.lstsq(design,center-TREF-cycle*SPACING,rcond=None)[0]
    epoch=TREF+pars[0]
    period=2*(SPACING+pars[1])
    residual=center-epoch-cycle*period/2-(cycle%2)*pars[2]
    formal=np.linalg.inv(design.T@design)*np.sum(residual**2)/(len(cycle)-3)
    stats=dict(epoch_even_BJD_TDB=epoch,period_days=period,formal_period_error_days=2*np.sqrt(formal[1,1]),
               odd_late_minutes=pars[2]*1440,ecos_omega_convention_even_primary=0.5*np.pi*pars[2]/period,
               ztf_points=len(ground),total_fitted_eclipses=len(times),sectors=results)
    write_rows("eclipse-times.csv",times)
    with (ROOT/"measurements.json").open("w") as handle: json.dump(stats,handle,indent=2)
    fig,ax=plt.subplots(1,2,figsize=(9,3.4),constrained_layout=True)
    for parity,color in [(0,"#21618c"),(1,"#aa482c")]:
        t,y=curves[-1]
        n=np.rint((t-epoch)/(period/2)).astype(int)
        dt=(t-epoch-n*period/2)*24
        use=(n%2==parity)&(abs(dt)<2)
        ax[0].scatter(dt[use],y[use],s=3,alpha=.24,color=color,label=("shallower","deeper")[parity])
    ax[0].set(xlabel="Hours from even half-period grid",ylabel="Normalized TESS flux",title="TESS sector 75")
    ax[0].legend(frameon=False,markerscale=2)
    for band,color in [("zg","#268a76"),("zr","#ac4635"),("zi","#62569e")]:
        rows=[r for r in ground if r["filtercode"]==band]
        phase=((np.array([float(r["hjd"]) for r in rows])-epoch)/period+.25)%1-.25
        ax[1].scatter(phase,[r["dm"] for r in rows],s=3,color=color,alpha=.46,label=band)
    ax[1].set(xlabel="Orbital phase",ylabel="Magnitude minus field median",title="Independent ZTF photometry",ylim=(.80,-.18))
    ax[1].legend(frameon=False,markerscale=2,ncol=3)
    fig.savefig(ROOT/"eclipses.png",dpi=150)
    print(json.dumps({k:v for k,v in stats.items() if k!="sectors"},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    main(parser.parse_args().cache)
