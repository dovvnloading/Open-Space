"""Blindly recover WD 1856 b in the latest archived control visit with custom target-pixel apertures."""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import ndtr

from ffi_acquire import propagated, tsv_rows
from ffi_photometry import gaia_t
from positive_controls import BASE, KNOWN, obtain
from vet_extension_pixels import sha
from wd_search import search

TIC = 267574918
GAIA = "2146576589564898688"


def profile(x, y, x0, y0, sigma):
    dx = ndtr((x+.5-x0)/sigma)-ndtr((x-.5-x0)/sigma)
    dy = ndtr((y+.5-y0)/sigma)-ndtr((y-.5-y0)/sigma)
    return dx*dy


def measure(item, row, field, cache):
    path = obtain(item["key"],item["size"],cache)
    sources = tsv_rows(field)
    with fits.open(path) as h:
        if h[0].header["TICID"]!=TIC or h[0].header["SECTOR"]!=row["sector"]:
            raise ValueError("Unexpected product")
        table = h[1].data
        cadence = round(h[1].header["TIMEDEL"]*86400)
        t = table["TIME"].astype(float)
        good = (table["QUALITY"]==0) & np.isfinite(t)
        year = 2000+(np.nanmedian(t[good])+2457000-2451545)/365.25
        wcs = WCS(h[2].header,naxis=2)
        y,x = np.indices(h[2].data.shape)
        target, = [s for s in sources if s["Source"]==GAIA]
        tx, ty = wcs.world_to_pixel_values(*propagated(target,year))
        others = [(s,wcs.world_to_pixel_values(*propagated(s,year))) for s in sources if s["Gmag"] and s["Source"]!=GAIA]

        def field_model(sigma):
            wd = profile(x,y,tx,ty,sigma)
            blend = sum(10**(-.4*(gaia_t(s)-gaia_t(target))) * profile(x,y,*xy,sigma) for s,xy in others)
            return wd, wd+blend

        signal, _ = field_model(.85)
        # Independent of event times or flux deficits. Robust pipeline pixel uncertainties
        # supply photon + sky variance; enumerate possible sums in the small target disk.
        variance = np.nanmedian(table["FLUX_ERR"][np.flatnonzero(good)[::100]].astype(float)**2,axis=0)
        disk = np.hypot(x-tx,y-ty)<=2
        pix = np.flatnonzero(disk & np.isfinite(variance) & (variance>0))
        subsets = (np.arange(1,2**len(pix))[:,None] >> np.arange(len(pix))) & 1
        snr = (subsets @ signal.ravel()[pix]) / np.sqrt(subsets @ variance.ravel()[pix])
        best = np.zeros_like(disk)
        best.ravel()[pix[subsets[np.argmax(snr)].astype(bool)]] = True
        apertures = {"radius_1p5":np.hypot(x-tx,y-ty)<=1.5, "flux_error_optimized":best}
        science = {}
        for name, mask in apertures.items():
            flux = np.sum(table["FLUX"][:,mask],axis=1,dtype=float)
            valid = good & np.isfinite(flux)
            result = search(t[valid],flux[valid],cadence)
            p0, epoch, _ = KNOWN[TIC]
            match = None
            for rank, peak in enumerate(result["periodic"],1):
                period, e = peak["period_days"], peak["epoch_btjd"]
                width = peak["duration_min"]/1440
                projected = e+round((np.median(t[valid])-e)/period)*period
                slip = abs((projected-epoch+p0/2)%p0-p0/2)
                if abs(period-p0)*np.ptp(t[valid])/period<=width and slip<=width:
                    match = dict(rank=rank, midpoint_phase_offset_minutes=1440*slip, **peak)
                    break
            fractions = []
            for sigma in (.65,.85,1.05):
                wd, blend = field_model(sigma)
                fractions.append(float(wd[mask].sum()/blend[mask].sum()))
            science[name] = dict(pixels_xy=np.stack((x[mask],y[mask]),axis=1).tolist(),
                                 gaia_model_fractions_065_085_105=fractions, usable=result["usable_cadences"],
                                 periodic_peaks=len(result["periodic"]), matching_peak=match,
                                 strongest_peak=next(iter(result["periodic"]),None))
    return dict(tic=TIC,gaia=GAIA,sector=row["sector"],target_xy=[float(tx),float(ty)],cadence_s=cadence,
                url=BASE+item["key"],bytes=item["size"],sha256=sha(path),apertures=science)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("inventory",type=Path);p.add_argument("field",type=Path)
    p.add_argument("cache",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args()
    rows=json.loads(a.inventory.read_text())["rows"]
    latest=max((r for r in rows if r["tic"]==TIC and any(q["key"].endswith("_tp.fits") for q in r["products"])),
               key=lambda r:r["sector"])
    a.cache.mkdir(parents=True,exist_ok=True)
    out=[measure(q,latest,a.field,a.cache) for q in latest["products"] if q["key"].endswith("tp.fits")]
    a.output.write_text(json.dumps(dict(inventory_sha256=sha(a.inventory),field_sha256=sha(a.field),products=out),
                                   allow_nan=False,indent=1)+"\n")
