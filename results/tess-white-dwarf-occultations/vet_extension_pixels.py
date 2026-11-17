"""Evaluate each frozen FFI ephemeris on mask-matched pixels and simultaneous SPOC SAP flux."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.special import ndtr

from wd_search import box_scores, flatten, mad, nearest


def sha(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def get_spoc(path, tic, sector, cadence, ffi_time, ffi_correction):
    with fits.open(path) as product:
        header = product[0].header
        if int(header["TICID"]) != int(tic) or int(header["SECTOR"]) != sector:
            raise ValueError(f"Unexpected target or sector: {path}")
        table = product[1].data
        interval = product[1].header["TIMEDEL"] * 86400
        if abs(interval - cadence) > 0.1:
            raise ValueError(f"Unexpected cadence: {path}")
        sc = table["TIME"].astype(float) - table["TIMECORR"]
        valid = np.isfinite(sc)
        sc = sc[valid]
        mapped_time = sc + np.interp(sc, ffi_time - ffi_correction, ffi_correction)
        flux = table["SAP_FLUX"][valid].astype(float)
        clean = (table["QUALITY"][valid] == 0) & np.isfinite(flux)
        return mapped_time, clean, flux


def pixel_flux(table, row):
    """Fixed, pixel-integrated target profile plus smooth sky and bright catalogue neighbours."""
    y, x = np.indices((9, 9))
    tx, ty = row["target_xy"]
    stamp = table["FLUX"].astype(float).reshape(-1, 81)
    def profile(x0, y0):
        dx = ndtr((x + .5 - x0)/.85) - ndtr((x - .5 - x0)/.85)
        dy = ndtr((y + .5 - y0)/.85) - ndtr((y - .5 - y0)/.85)
        return dx * dy
    u, v = (x-4)/4, (y-4)/4
    cols = [profile(tx, ty), np.ones_like(x), u, v, u*u, u*v, v*v]
    modeled = []
    positions = [(tx, ty)]
    for source, nx, ny, _, ratio in row["neighbours_top5"]:
        if (0 <= nx <= 8 and 0 <= ny <= 8 and ratio >= .001 and
                all(np.hypot(nx-xx, ny-yy) >= 1 for xx, yy in positions)):
            cols.append(profile(nx, ny))
            positions.append((nx, ny))
            modeled.append(source)
    design = np.asarray(cols).reshape(len(cols), -1).T
    weights = np.linalg.pinv(design)[0]
    forced = stamp @ weights
    ap = np.hypot(x-tx, y-ty) <= 1.5
    edge = (x==0) | (x==8) | (y==0) | (y==8)
    native = np.sum(stamp[:, ap.ravel()], axis=1) - ap.sum()*np.median(stamp[:, edge.ravel()], axis=1)
    return native, forced, modeled


class FixedSearch:
    def __init__(self, t, flux, dt):
        valid = np.isfinite(t+flux)
        t, index = np.unique(t[valid], return_index=True)
        f = flux[valid][index]
        trend = flatten(t, f)
        valid = np.isfinite(trend)
        self.time = t[valid]
        self.residual = f[valid]/trend[valid] - 1
        self.noise = mad(self.residual)
        self.dt = dt / 86400
        self.by_width = {}

    def at(self, event):
        width = event["duration_min"]/1440
        t = self.time
        if width not in self.by_width:
            self.by_width[width] = box_scores(t, self.residual, self.dt, width, self.noise)
        depth, scores = self.by_width[width]
        epoch = event["epoch_btjd"]
        period = event.get("period_days")
        if period:
            cycles = np.arange(np.ceil((t[0]-epoch)/period)-1, np.floor((t[-1]-epoch)/period)+2)
            centers = epoch + cycles*period
        else:
            centers = np.asarray([epoch])
        i = nearest(t, centers)
        support = (abs(t[i]-centers) <= .53*self.dt) & np.isfinite(scores[i])
        z = scores[i[support]]
        n = len(z)
        result = dict(covered=n, above_2p5=int(sum(z >= 2.5)), passes=False)
        if n:
            score = float(sum(z)/np.sqrt(n))
            result.update(score=score, depth=float(np.mean(depth[i[support]])))
            if period:
                other = float((sum(z)-max(z))/np.sqrt(max(1,n-1)))
                result.update(leave_deepest_out=other, passes=bool(n>=3 and sum(z>=2.5)>=3 and score>=7 and other>=4))
            else:
                result["passes"] = score>=7
        return result


def star_checks(row, rosters, cache):
    sector = row["sector"]
    tic = row["tic"]
    path = cache/f"{row['source_id']}-s{sector}.fits"
    with fits.open(path) as product:
        d = product[1].data
        for field, expected in row["arrays_sha256"].items():
            if hashlib.sha256(d[field].tobytes()).hexdigest()!=expected:
                raise ValueError(f"Original science array changed: {path.name} {field}")
        t, corr = d["TIME"].astype(float), d["TIMECORR"].astype(float)
        ok = np.isfinite(t+corr)
        native, forced, neighbors = pixel_flux(d, row)
        use = ok & (d["QUALITY"] == 0) & np.isfinite(native)
    dt = row["cadence_s"]
    products = {}
    data = {}
    for cadence in (20, 120):
        item, = [r for r in rosters if r["tic"]==tic and int(r["sector"])==sector and int(r["cadence_s"])==cadence]
        filename = item["filename"].replace("-tp.fits", "-lc.fits").replace("_tp.fits", "_lc.fits")
        file = cache / filename
        data[cadence] = get_spoc(file, tic, sector, cadence, t[ok], corr[ok])
        products[cadence] = dict(filename=filename, sha256=sha(file))
    t20, q20, _ = data[20]
    left = np.searchsorted(t20[q20],t-dt/172800)
    right = np.searchsorted(t20[q20],t+dt/172800)
    matched = use & ((right-left)>=8)
    filters = dict(ffi_original=FixedSearch(t[use],native[use],dt),
                   ffi_clean=FixedSearch(t[matched],native[matched],dt),
                   ffi_profile_clean=FixedSearch(t[matched],forced[matched],dt))
    for c, (ts, good, f) in data.items():
        filters[f"spoc_{c}"] = FixedSearch(ts[good],f[good],c)
    results = []
    for kind in ("periodic","isolated"):
        for e in row["search"][kind]:
            checks = {k:search.at(e) for k, search in filters.items()}
            measured = checks["ffi_original"]
            # The discovery search fits nearby periods/phases. Re-evaluate exactly its saved ephemeris.
            if not measured["passes"] or abs(measured["score"]-e["score"])>1e-3:
                raise ValueError(f"Original ephemeris failed replay for TIC {tic}: {e} {measured}")
            results.append(dict(kind=kind, **e, checks=checks))
    return dict(tic=tic,gaia=row["source_id"],sector=sector,ffi_file=path.name,
                ffi_sha256=sha(path),spoc_products=products,modeled_neighbors=neighbors,
                ffi_quality_zero=int(use.sum()),ffi_with_fast_support=int(matched.sum()),results=results)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("extension",type=Path);p.add_argument("inventory",type=Path)
    p.add_argument("cache",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args()
    rows=json.loads(a.extension.read_text())
    with a.inventory.open() as f:
        inventory=list(csv.DictReader(f))
    results=[star_checks(r,inventory,a.cache) for r in rows if r.get("status")=="included"
             and (r["search"]["periodic"] or r["search"]["isolated"])]
    payload=json.dumps(dict(input_sha256=sha(a.extension),products_sha256=sha(a.inventory),stars=results),
                       allow_nan=False,separators=(",",":")).encode()
    a.output.write_bytes(gzip.compress(payload,mtime=0))
