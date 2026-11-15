"""Is there still hot dust in the seven flaring nuclei? A SPHEREx colour test.

Usage: python spherex.py OUT.json [--cache DIR]

For each of the seven unreported flares this reads the 2.95-5.0 um planes of
the public SPHEREx DR1 HiPS spectral cube tile that contains it (byte ranges
only), and measures every 2MRS galaxy in the same tiles. In each channel the
nucleus value is the nearest 12.9 arcsec pixel minus the median of an annulus
5-8 pixels out. A W1-like (3.0-3.8 um) and a W2-like (4.1-5.0 um) mean give a
colour s = -2.5 log10(W1-like / W2-like).

Ordinary galaxies calibrate the colour: for 2MRS galaxies in the same tiles
with AllWISE W1-W2 below 0.5 mag, delta = s - (W1-W2 - 0.639), where 0.639 mag
is 2.5 log10 of the ratio of the W1 and W2 zero points. A flare nucleus with
residual hot dust should have s redder than its quiescent NEOWISE W1-W2
predicts: excess = s - (W1-W2_quiet - 0.639 + median delta).
"""
import argparse
import csv
import gzip
import json
import math
import pathlib
import time

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS

from spherex_tools import hdu_layout, healpix_nest, ranged, tile_url

FLARES = ["02265317+4141211", "06440345-5554310", "11514275-3116072", "21172918-8048085",
          "19265088-1147169", "03183820+3736290", "01110461-4558435"]
SOURCE_BIT = 0x200000          # "known source" flag, harmless here
RING = (5.0, 8.0)              # background ring, pixels (set by --ring)
W1_RANGE, W2_RANGE = (3.0, 3.8), (4.1, 5.0)
ZP_OFFSET = 2.5 * math.log10(309.54 / 171.787)
GATOR = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"


def fetch_planes(url, layout, k0, k1, cache):
    """IMAGE and FLAGS planes k0..k1 (inclusive) as arrays [k, y, x]."""
    out = {}
    for name, dtype in (("IMAGE", ">f4"), ("FLAGS", ">i4")):
        hdu = next(h for h in layout if h["name"] == name)
        nx, ny = hdu["shape"][0], hdu["shape"][1]
        plane = nx * ny * 4
        path = None if cache is None else cache / f"{url.rsplit('/', 1)[1]}.{name}.{k0}-{k1}.npy"
        if path is not None and path.exists():
            out[name] = np.load(path)
            continue
        raw = ranged(url, hdu["data_offset"] + k0 * plane, hdu["data_offset"] + (k1 + 1) * plane - 1)
        arr = np.frombuffer(raw, dtype=dtype).reshape(k1 - k0 + 1, ny, nx).astype(np.float32 if name == "IMAGE"
                                                                                   else np.int64)
        if path is not None:
            np.save(path, arr)
        out[name] = arr
    return out["IMAGE"], out["FLAGS"]


def channels(url, layout):
    table = layout[-1]
    raw = ranged(url, table["header_offset"], table["data_offset"] + table["data_bytes"] - 1)
    import io
    buf = io.BytesIO()
    fits.PrimaryHDU().writeto(buf)
    with fits.open(io.BytesIO(buf.getvalue() + raw)) as h:
        return np.array(h[1].data["WAVELENGTH"], float)


def measure(img, flg, wav, x, y):
    """Nucleus-minus-annulus values per channel (MJy/sr), NaN where flagged or missing."""
    ny, nx = img.shape[1:]
    ix, iy = int(round(x)), int(round(y))
    if not (0 <= ix < nx and 0 <= iy < ny):
        return None
    y0, y1, x0, x1 = max(iy - 8, 0), min(iy + 9, ny), max(ix - 8, 0), min(ix + 9, nx)   # ring may be cut at a tile edge
    yy, xx = np.mgrid[y0:y1, x0:x1]
    r = np.hypot(xx - x, yy - y)
    ring = (r >= RING[0]) & (r <= RING[1])
    vals = np.full(len(wav), np.nan)
    for k in range(len(wav)):
        c = img[k, iy, ix]
        if not np.isfinite(c) or (flg[k, iy, ix] & ~SOURCE_BIT):
            continue
        stamp, fstamp = img[k, y0:y1, x0:x1], flg[k, y0:y1, x0:x1]
        good = ring & np.isfinite(stamp) & ((fstamp & ~SOURCE_BIT) == 0)
        if good.sum() >= 20:
            vals[k] = c - np.median(stamp[good])
    return vals


def colour(vals, wav):
    m1 = (wav >= W1_RANGE[0]) & (wav <= W1_RANGE[1]) & np.isfinite(vals)
    m2 = (wav >= W2_RANGE[0]) & (wav <= W2_RANGE[1]) & np.isfinite(vals)
    if m1.sum() < 3 or m2.sum() < 3:
        return None, None, None
    f1, f2 = vals[m1].mean(), vals[m2].mean()
    if f1 <= 0 or f2 <= 0:
        return f1, f2, None
    return f1, f2, -2.5 * math.log10(f1 / f2)


def allwise(rows):
    upload = "|   ra         |   dec        |\n|   double     |   double     |\n"
    upload += "".join(f" {r['ra']:13.7f}  {r['dec']:13.7f} \n" for r in rows)
    for attempt in range(5):
        resp = requests.post(GATOR, data={"catalog": "allwise_p3as_psd", "spatial": "Upload", "uradius": "3",
                                          "uradunits": "arcsec", "outfmt": "1",
                                          "selcols": "w1mpro,w1sigmpro,w2mpro,w2sigmpro,ext_flg,cc_flags"},
                             files={"filename": ("positions.tbl", upload)}, timeout=600)
        if resp.status_code == 200 and "\\RowsRetrieved" in resp.text:
            break
        time.sleep(20 * (attempt + 1))
    else:
        raise IOError("Gator failed: " + resp.text[:200])
    lines = [ln for ln in resp.text.splitlines() if ln and not ln.startswith("\\")]
    header = [h.strip() for h in lines[0].strip("|").split("|")]
    best = {}
    for ln in lines[4:]:
        cells = ln.split()
        rec = dict(zip(header, cells[-len(header):] if len(cells) >= len(header) else cells))
        try:
            idx = int(rec["cntr_01"]) - 1
            dist = float(rec["dist_x"])
            w1, w2 = float(rec["w1mpro"]), float(rec["w2mpro"])
        except (KeyError, ValueError):
            continue
        if idx not in best or dist < best[idx][0]:
            best[idx] = (dist, w1, w2, rec.get("cc_flags", ""))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out")
    ap.add_argument("--cache", type=pathlib.Path, default=None)
    ap.add_argument("--ring", type=float, nargs=2, default=(5.0, 8.0),
                    help="inner and outer radius of the background ring, pixels")
    args = ap.parse_args()
    global RING
    RING = tuple(args.ring)
    if args.cache:
        args.cache.mkdir(parents=True, exist_ok=True)
    flares = {r["id"]: r for r in csv.DictReader(open("data/flares.csv", encoding="utf-8"))}
    targets = [flares[i] for i in FLARES]
    galaxies = [r for r in csv.DictReader(gzip.open("data/2mrs-positions.csv.gz", "rt", encoding="utf-8"))]
    tiles = sorted({healpix_nest(5, float(t["ra"]), float(t["dec"])) for t in targets})
    rows = []
    for tile in tiles:
        url = tile_url(tile)
        layout = hdu_layout(url)
        wav = channels(url, layout)
        ks = np.where((wav >= 2.95) & (wav <= 5.0))[0]
        k0, k1 = int(ks.min()), int(ks.max())
        img, flg = fetch_planes(url, layout, k0, k1, args.cache)
        w = wav[k0:k1 + 1]
        image_hdu = next(h for h in layout if h["name"] == "IMAGE")
        header = fits.Header.fromstring(image_hdu["header_text"])
        for key in list(header.keys()):           # drop the tabulated wavelength axis; keep the sky axes
            if key.endswith("3") and key[:2] in ("CT", "CU", "CR", "CD", "CN") or key.startswith(("PC3_", "PV3_",
                                                                                                   "PS3_")) \
                    or key in ("PC1_3", "PC2_3"):
                del header[key]
        header["NAXIS"], header["WCSAXES"] = 2, 2
        wcs = WCS(header)
        members = [g for g in galaxies if healpix_nest(5, float(g["ra"]), float(g["dec"])) == tile]
        for g in members:
            x, y = wcs.all_world2pix([[float(g["ra"]), float(g["dec"])]], 0)[0]
            vals = measure(img, flg, w, x, y)
            if vals is None:
                continue
            f1, f2, s = colour(vals, w)
            rows.append({"id": g["id"], "ra": float(g["ra"]), "dec": float(g["dec"]), "tile": tile,
                         "pixel": [round(float(x), 2), round(float(y), 2)],
                         "is_target": g["id"] in FLARES, "is_flare": g["id"] in flares,
                         "w1like_mjy_sr": None if f1 is None else float(f1),
                         "w2like_mjy_sr": None if f2 is None else float(f2),
                         "s": s, "n_channels": int(np.isfinite(vals).sum()),
                         "spectrum": [None if not np.isfinite(v) else round(float(v), 5) for v in vals]
                         if g["id"] in FLARES else None})
        print(f"tile {tile}: {len(members)} 2MRS galaxies, wavelengths {w.min():.3f}-{w.max():.3f} um", flush=True)
    wise = allwise(rows)
    for i, r in enumerate(rows):
        if i in wise:
            r["allwise_w1_w2"] = wise[i][1] - wise[i][2]
            r["allwise_cc_flags"] = wise[i][3]
    controls = [r for r in rows if not r["is_flare"] and r["s"] is not None and r.get("allwise_w1_w2") is not None
                and r["allwise_w1_w2"] < 0.5 and r.get("allwise_cc_flags", "0000") == "0000"]
    deltas = np.array([r["s"] - (r["allwise_w1_w2"] - ZP_OFFSET) for r in controls])
    offset = float(np.median(deltas))
    scatter = float(1.4826 * np.median(np.abs(deltas - offset)))
    results = []
    for t in targets:
        r = next((x for x in rows if x["id"] == t["id"]), None)
        q = float(t["quiescent_w1_w2"])
        res = {"id": t["id"], "host": t["host"], "quiescent_w1_w2": q, "tile": healpix_nest(5, float(t["ra"]),
                                                                                         float(t["dec"]))}
        if r is None or r["s"] is None:
            res["status"] = "not measured"
        else:
            pred = q - ZP_OFFSET + offset
            res.update({"s_observed": r["s"], "s_predicted": pred, "excess_mag": r["s"] - pred,
                        "excess_sigma": (r["s"] - pred) / scatter, "w1like_mjy_sr": r["w1like_mjy_sr"],
                        "w2like_mjy_sr": r["w2like_mjy_sr"], "n_channels": r["n_channels"],
                        "pixel": r["pixel"], "status": "measured"})
            if r.get("allwise_w1_w2") is not None:     # same catalogue as the controls
                pa = r["allwise_w1_w2"] - ZP_OFFSET + offset
                ex = r["s"] - pa
                res.update({"allwise_w1_w2": r["allwise_w1_w2"], "excess_mag_allwise": ex,
                            "excess_sigma_allwise": ex / scatter,
                            "controls_at_least_as_red": int(np.sum(deltas - offset >= ex))})
        results.append(res)
        print(res, flush=True)
    out = {"product": "SPHEREx DR1 HiPS cube, hips_20261026, order 5 (12.9 arcsec pixels)",
           "tiles": tiles, "w1_like_um": W1_RANGE, "w2_like_um": W2_RANGE, "background_ring_pixels": RING,
           "controls": len(controls), "control_offset_mag": offset, "control_scatter_mag": scatter,
           "targets": results, "galaxies": rows}
    with open(args.out, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1)
    print(f"controls {len(controls)}, offset {offset:.3f}, scatter {scatter:.3f}")


if __name__ == "__main__":
    main()
