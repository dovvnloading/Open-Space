"""Force PSF photometry along conditional tracks on downloaded ZTF subtraction images."""

from ztf_meta import c, HERE, CACHE
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage as nd
from pathlib import Path
import numpy as np, csv, json, hashlib


def mad(x):
    return 1.48260222 * np.median(np.abs(x - np.median(x)))


def checked(tag):
    rec = CACHE / (tag + ".json")
    path = CACHE / (tag + ".raw")
    if not (rec.exists() and path.exists()):
        return None
    row = json.loads(rec.read_text())
    return (
        path
        if row.get("sha256") == hashlib.sha256(path.read_bytes()).hexdigest()
        else None
    )


def sample(img, mask, x, y, ps):
    ix, iy = int(round(x)), int(round(y))
    h = ps.shape[0] // 2
    sl = np.s_[iy - h : iy + h + 1, ix - h : ix + h + 1]
    a = img[sl]
    dq = mask[sl]
    if a.shape != ps.shape:
        return None
    yy, xx = np.mgrid[-h : h + 1, -h : h + 1]
    rr = np.hypot(xx, yy)
    if np.any(dq[rr < 4] != 0) or not np.isfinite(a).all():
        return None
    k = nd.shift(ps, (y - iy, x - ix), order=1, mode="constant")
    k /= k.sum()
    sky = np.median(a[rr > 8])
    return a - sky, k


def run(name):
    rng = np.random.default_rng(38671)
    s = c.circular_family(name)
    n = len(s[0])
    records = {}
    visits = []
    mags = [18.0, 19.0, 20.0, 20.5, 21.0, 21.5, 22.0]
    sel = json.loads((HERE / "ztf_selection.json").read_text())[name]
    pinned = HERE / "inputs.json"
    labels = (
        {
            r["label"]
            for r in json.loads(pinned.read_text())
            if r["required_for_reduction"]
        }
        if pinned.exists()
        else None
    )
    for row in sel:
        if labels is not None and "ztf_" + row["pid"] + "_diff" not in labels:
            continue
        tag = "ztf_" + row["pid"] + "_"
        paths = [checked(tag + k) for k in ["diff", "mask", "psf"]]
        if any(q is None for q in paths):
            continue
        im, h = fits.getdata(paths[0], header=True)
        dq = fits.getdata(paths[1])
        pk = fits.getdata(paths[2]).astype(float)
        pk /= pk.sum()
        assert (
            im.shape == dq.shape and pk.shape[0] == pk.shape[1] and pk.shape[0] % 2 == 1
        )
        t = (
            c.epoch(Time(float(h["OBSJD"]), format="jd", scale="utc").isot)
            + h["EXPTIME"] / 172800
        )
        tobs = Time(float(h["OBSJD"]), format="jd", scale="utc").isot
        ra, de = c.other_sky(s, t)
        w = WCS(h)
        x, y = w.world_to_pixel_values(ra, de)
        ny, nx = im.shape
        good = (x > 32) & (y > 32) & (x < nx - 32) & (y < ny - 32)
        noise = []
        for xx, yy in zip(
            rng.uniform(32, nx - 32, 1000), rng.uniform(32, ny - 32, 1000)
        ):
            val = sample(im, dq, xx, yy, pk)
            if val is not None:
                a, k = val
                noise.append(np.sum(a * k) / np.sum(k * k))
        sigma = mad(noise)
        center = np.median(noise)
        zp = h["MAGZP"]
        factor = 10 ** (0.4 * (22.5 - zp))
        used = 0
        for j in np.flatnonzero(good):
            val = sample(im, dq, x[j], y[j], pk)
            if val is None:
                continue
            a, k = val
            den = np.sum(k * k)
            signal = (np.sum(a * k) / den - center) * factor
            fake = []
            for mag in mags:
                flux = 10 ** (0.4 * (zp - mag))
                planted = a + flux * k
                fake.append((np.sum(planted * k) / den - center) * factor)
            records.setdefault(int(j), []).append((signal, sigma * factor, fake))
            used += 1
        visits.append(
            dict(
                pid=row["pid"],
                start=tobs,
                filter=row["filtercode"],
                covered=int(good.sum()),
                usable=used,
                zero_point=float(zp),
                blank_sigma_dn=float(sigma),
                controls=len(noise),
            )
        )
    eligible = 0
    flags = []
    recovered = np.zeros(len(mags), int)
    peak = None

    def decision(fs, es):
        score = np.sum(fs / es**2) / np.sqrt(np.sum(1 / es**2))
        return score >= 5.5 and np.count_nonzero(fs / es >= 2.5) >= 2, score

    for j, r in records.items():
        if len(r) < 2:
            continue
        eligible += 1
        fs, es = np.array([x[0] for x in r]), np.array([x[1] for x in r])
        hit, snr = decision(fs, es)
        peak = float(snr) if peak is None else max(peak, float(snr))
        if hit:
            flags.append(dict(track=j, snr=float(snr)))
        for q in range(len(mags)):
            recovered[q] += decision(np.array([x[2][q] for x in r]), es)[0]
    result = dict(
        target=name,
        grid_nodes=n,
        at_least_two_valid=eligible,
        n_visits_used=len(visits),
        visits=visits,
        max_joint_snr=peak,
        flags=flags,
        injections=[
            dict(ztf_r0_ab=mag, draws=n, eligible=eligible, recovered=int(num))
            for mag, num in zip(mags, recovered)
        ],
    )
    print(name, eligible, len(visits), flags, recovered, flush=True)
    return result


if __name__ == "__main__":
    result = {nm: run(nm) for nm in ["RR", "Chen1", "Chen2"]}
    (HERE / "ztf_results.json").write_text(json.dumps(result, indent=2) + "\n")
