"""Three-epoch forced PSF flux stack at each explicit circular Chen orbit.

The two single-day detections do not supply a heading. Each track here is one discrete,
conditional trajectory, not an area-integrated circular-annulus limit.
"""

from pathlib import Path
import csv, json
import numpy as np
from scipy import ndimage as nd
from scipy.spatial import cKDTree
from astropy.io import fits
from astropy.wcs import WCS
from ztf_meta import c, HERE, CACHE

hh = 17
yy, xx = np.mgrid[-hh : hh + 1, -hh : hh + 1]
rad = np.hypot(xx, yy)
ann = rad >= 13
core = rad < 9


def mad(v):
    return 1.48260222 * np.median(np.abs(v - np.median(v)))


def stars_on_ccd(im, dq, w, table):
    ny, nx = im.shape
    ra, dc = w.pixel_to_world_values(nx / 2, ny / 2)
    sel = (np.abs(table[:, 0] - ra) < 0.29) & (np.abs(table[:, 1] - dc) < 0.25)
    sky = table[sel]
    x, y = w.world_to_pixel_values(sky[:, 0], sky[:, 1])
    inb = (x > 20) & (x < nx - 20) & (y > 20) & (y < ny - 20)
    stamps = []
    shifts = []
    zps = []
    for a, b, r in zip(x[inb], y[inb], sky[inb, 2]):
        ix, iy = int(np.rint(a)), int(np.rint(b))
        sl = np.s_[iy - hh : iy + hh + 1, ix - hh : ix + hh + 1]
        if dq[sl][core].any():
            continue
        q = im[sl]
        q = q - np.median(q[ann])
        wt = (rad < 6) * np.clip(q, 0, None)
        if not wt.sum():
            continue
        dx = (wt * xx).sum() / wt.sum()
        dy = (wt * yy).sum() / wt.sum()
        flux = q[rad < 10].sum()
        if np.hypot(ix + dx - a, iy + dy - b) > 2 or flux <= 0:
            continue
        centered = nd.shift(q, (-dy, -dx), order=3)
        stamps.append(centered / flux)
        shifts.append([ix + dx - a, iy + dy - b])
        zps.append(22.5 - 2.5 * np.log10(r) + 2.5 * np.log10(flux))
    if len(stamps) < 3:
        return None
    zps = np.array(zps)
    keep = np.abs(zps - np.median(zps)) < max(0.10, mad(zps) * 3)
    if keep.sum() < 3:
        return None
    ps = np.clip(np.median(np.array(stamps)[keep], axis=0), 0, None) * (rad < 11)
    ps /= ps.sum()
    return (
        ps,
        np.median(np.array(shifts)[keep], axis=0),
        np.median(zps[keep]),
        mad(zps[keep]),
        int(keep.sum()),
    )


def stamp_measure(im, dq, x, y, ps, plant=None):
    ix, iy = int(np.rint(x)), int(np.rint(y))
    sl = np.s_[iy - hh : iy + hh + 1, ix - hh : ix + hh + 1]
    q = im[sl]
    mask = dq[sl] != 0
    if q.shape != rad.shape or np.any(mask[core]) or np.sum(~mask & ann) < 100:
        return None
    k = nd.shift(ps, (y - iy, x - ix), order=1, mode="constant")
    den = np.sum(k * k)
    sky = np.median(q[ann & ~mask])
    return (np.sum((q - sky) * k) / den, k, q - sky)


def measure(name):
    rng = np.random.default_rng(47093)
    catalog = np.loadtxt(CACHE / ("ls_" + name + ".raw"), delimiter=",", skiprows=1)
    stat = np.loadtxt(
        CACHE / ("stationary_" + name + ".raw"), delimiter=",", skiprows=1
    )
    fac = np.cos(np.deg2rad(np.mean(stat[:, 1])))
    tree = cKDTree(stat * [fac, 1])
    s = c.circular_family(name)
    n = len(s[0])
    exp = json.loads((HERE / "selected.json").read_text())[name]
    records = [{} for _ in exp]
    ccdrows = []
    mags = [19.0, 20.0, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0]
    ids = rng.choice(n, size=280, replace=False)
    idset = set(ids)
    coords = []
    for k, row in enumerate(exp):
        fid = row["md5sum"]
        science = fits.open(CACHE / (fid + ".raw"))
        masks = fits.open(CACHE / (fid + "_mask.raw"))
        hdr = science[0].header
        t = c.epoch(hdr["DATE-OBS"]) + hdr["EXPTIME"] / 172800
        ra, dec = c.other_sky(s, t)
        coords.append([ra, dec])
        sep = tree.query(np.c_[ra * fac, dec])[0] * 3600
        # The DR10 catalogs predate the modern visits. For pre-2023 exposures, shared
        # images could bias a counterpart veto, so only use masks derived independently.
        use_static_veto = hdr["DATE-OBS"] >= "2023"
        for ext in range(1, len(science)):
            h = science[ext].header
            w = WCS(h)
            nx, ny = h["NAXIS1"], h["NAXIS2"]
            x, y = w.world_to_pixel_values(ra, dec)
            inside = (x > 40) & (x < nx - 40) & (y > 40) & (y < ny - 40)
            if not inside.any():
                continue
            im = science[ext].data
            dq = masks[ext].data
            assert masks[ext].header["CCDNUM"] == h["CCDNUM"]
            cal = stars_on_ccd(im, dq, w, catalog)
            if cal is None:
                continue
            ps, shift, zp, sc, ncal = cal
            x += shift[0]
            y += shift[1]
            control = []
            for ix, iy in zip(
                rng.uniform(32, nx - 32, 155), rng.uniform(32, ny - 32, 155)
            ):
                q = stamp_measure(im, dq, ix, iy, ps)
                if q is not None:
                    control.append(q[0])
            if len(control) < 60:
                continue
            base = np.median(control)
            noise = mad(control)
            factor = 10 ** (0.4 * (22.5 - zp))
            ngood = 0
            for j in np.flatnonzero(inside):
                if use_static_veto and sep[j] < 2.2:
                    continue
                a = stamp_measure(im, dq, x[j], y[j], ps)
                if a is None:
                    continue
                flux, kernel, local = a
                fake = []
                if j in idset:
                    for mag in mags:
                        fs = 10 ** (0.4 * (zp - mag))
                        planted = local + rng.poisson(fs * kernel * 3.8) / 3.8
                        f = (
                            np.sum((planted - np.median(planted[ann])) * kernel)
                            / np.sum(kernel**2)
                            - base
                        ) * factor
                        fake.append(f)
                records[k][j] = dict(
                    f=(flux - base) * factor,
                    e=noise * factor,
                    fake=fake,
                    sep=float(sep[j]),
                )
                ngood += 1
            ccdrows.append(
                dict(
                    target=name,
                    epoch=k,
                    ccd=h["CCDNUM"],
                    ncal=ncal,
                    zp_ab=round(zp, 5),
                    zp_scatter=round(sc, 5),
                    usable=ngood,
                    noise_nmgy=round(noise * factor, 6),
                )
            )
        science.close()
        masks.close()
        print(name, k, len(records[k]), flush=True)

    def test(fs, es):
        if len(fs) < 2:
            return False, np.nan
        fs = np.array(fs)
        es = np.array(es)
        joint = np.sum(fs / es**2) / np.sqrt(np.sum(1 / es**2))
        return joint >= 5.5 and np.sum(fs / es >= 2.5) >= 2, joint

    ff = []
    eligible = 0
    counts = np.zeros(len(mags), int)
    nplant = 0
    max_sc = -1e20
    for j in range(n):
        r = [r[j] for r in records if j in r]
        if len(r) < 2:
            continue
        eligible += 1
        es = [o["e"] for o in r]
        fire, sn = test([o["f"] for o in r], es)
        max_sc = max(sn, max_sc)
        if fire:
            ff.append(
                dict(
                    track=j,
                    snr=sn,
                    nobs=len(r),
                    ra_dec=[(rr[0][j], rr[1][j]) for rr in coords],
                    per_epoch_snr=[
                        round(rec[j]["f"] / rec[j]["e"], 3) if j in rec else None
                        for rec in records
                    ],
                    nearest_dr10_arcsec=min(q["sep"] for q in r),
                )
            )
        if j not in idset:
            continue
        nplant += 1
        for l, _ in enumerate(mags):
            counts[l] += test([q["fake"][l] for q in r], es)[0]
    result = dict(
        target=name,
        circular_tracks=n,
        at_least_two_eligible=eligible,
        nflag=len(ff),
        max_joint=float(max_sc),
        individual=[len(x) for x in records],
        n_plant_draws=len(ids),
        n_plant_eligible=nplant,
        injections=[dict(r_AB=m, recovered=int(v)) for m, v in zip(mags, counts)],
    )
    (HERE / (name + "_result.json")).write_text(json.dumps(result, indent=2) + "\n")
    (HERE / (name + "_flags.json")).write_text(json.dumps(ff, indent=2) + "\n")
    with (HERE / (name + "_calibration.csv")).open("w", newline="") as f:
        wr = csv.DictWriter(f, list(ccdrows[0]))
        wr.writeheader()
        wr.writerows(ccdrows)
    print(result, flush=True)


if __name__ == "__main__":
    import sys

    for nm in sys.argv[1:] or ["Chen1", "Chen2"]:
        measure(nm)
