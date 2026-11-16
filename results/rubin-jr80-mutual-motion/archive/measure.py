"""Field-star-calibrated two-component fits to CADC/NOIRLab archival image cutouts."""

import glob, json, warnings, sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings("ignore", module="astropy")
W = 8
yy, xx = np.mgrid[-W : W + 1, -W : W + 1]


def profile(x, y, s, t, b, theta):
    a = np.cos(theta) * x + np.sin(theta) * y
    c = -np.sin(theta) * x + np.cos(theta) * y
    return (1 + (a / s) ** 2 + (c / t) ** 2) ** (-b)


def fitshape(stamp):
    bg = np.median(stamp)
    peak = np.unravel_index(np.argmax(stamp), stamp.shape)
    init = [bg, stamp.max() - bg, peak[1] - W, peak[0] - W, 2.7, 2.4, 2.4, 0]

    def fun(p):
        bg, A, x, y, s, t, b, angle = p
        return (bg + A * profile(xx - x, yy - y, s, t, b, angle) - stamp).ravel()

    p = least_squares(
        fun,
        init,
        bounds=(
            [-np.inf, 0, -W, -W, 0.65, 0.65, 1.05, -7],
            [np.inf, np.inf, W, W, 15, 15, 10, 7],
        ),
        max_nfev=150,
    )
    bg, A, x, y, s, t, b, angle = p.x
    return p.x, np.sum(p.fun**2)


def fieldshape(stars):
    sh = []
    for st in stars:
        v, r = fitshape(st["d"])
        if abs(v[2]) > 2 or abs(v[3]) > 2:
            continue
        f = 2 * np.sqrt(2 ** (1 / v[6]) - 1)
        s, t = v[4:6] * f
        if not 1.5 < min(s, t) < 8 or max(s, t) / min(s, t) > 1.45:
            continue
        a = v[7]
        g = np.array([np.cos(a), np.sin(a)])
        mat = np.outer(g, g) * (s * s - t * t) + np.eye(2) * t * t
        sh.append((mat, v[6], st))
    if len(sh) < 2:
        return None
    area = np.array([np.sqrt(np.linalg.det(i[0])) for i in sh])
    med = np.median(area)
    sh = [i for i, a in zip(sh, area) if 0.68 * med < a < 1.32 * med]
    mat = np.median([i[0] for i in sh], axis=0)
    b = np.median([i[1] for i in sh])
    eig, ve = np.linalg.eigh(mat)
    theta = np.arctan2(*ve[:, 1][::-1])
    return (*np.sqrt(eig[::-1]) / (2 * np.sqrt(2 ** (1 / b) - 1)), b, theta), [
        i[2] for i in sh
    ]


def fit_source(stamp, shape, scale, noise):
    bg = np.median(stamp)
    amp = stamp.max() - bg
    y, x = np.unravel_index(np.argmax(stamp), stamp.shape)
    x -= W
    y -= W

    def basic(p, dual):
        base, A, x, y = p[:4]
        dx, dy = p[4:6] / 2 if dual else [0, 0]
        q = p[6] if dual else 1
        # A is the total scale; q is the secondary/primary flux ratio.
        mod = (
            profile(xx - x - dx, yy - y - dy, *shape)
            + q * profile(xx - x + dx, yy - y + dy, *shape)
        ) / (1 + q)
        # diagonal sky noise + Poisson source noise: use positive model-independent data estimates
        return (
            (base + A * mod - stamp) / np.sqrt(noise**2 + np.maximum(stamp - bg, 0))
        ).ravel()

    lo = [-np.inf, 0.01, -W, -W]
    hi = [np.inf, np.inf, W, W]
    p0 = least_squares(
        lambda p: basic(p, False), [bg, max(0.02, amp), x, y], bounds=(lo, hi)
    )
    bound = 1.45 / scale
    best = None
    for a in (0, np.pi / 4, np.pi / 2, np.pi * 3 / 4):
        p = least_squares(
            lambda p: basic(p, True),
            [
                *np.clip(p0.x, lo, hi),
                0.6 / scale * np.cos(a),
                0.6 / scale * np.sin(a),
                0.88,
            ],
            bounds=(lo + [-bound, -bound, 0.65], hi + [bound, bound, 1]),
        )
        if best is None or p.cost < best.cost:
            best = p
    return best.x, 2 * (p0.cost - best.cost), 2 * best.cost / (stamp.size - 7)


def measure(file):
    file = Path(file)
    meta = json.loads(file.with_suffix(".json").read_text())
    ra, dec = map(float, (meta["Object_RA"], meta["Object_Dec"]))
    with fits.open(file) as hd:
        for h in hd:
            if h.data is None or h.data.ndim != 2:
                continue
            img = h.data.astype(float)
            ny, nx = img.shape
            w = WCS(h.header)
            cx, cy = w.all_world2pix([[ra, dec]], 0)[0]
            if not (10 < cx < nx - 10 and 10 < cy < ny - 10):
                continue
            break
        else:
            return {"id": file.stem, "status": "no_coverage"}
    valid = img[np.isfinite(img)]
    bg = np.median(valid)
    noise = 1.4826 * np.median(abs(valid - bg))
    im = np.where(np.isfinite(img), img, bg)
    sm = ndimage.gaussian_filter(im, 1.1)
    mask = (sm == ndimage.maximum_filter(sm, 17)) & (sm > bg + 6 * noise) & (sm < 50000)
    ys, xs = np.where(mask)
    stamps = []
    for x, y in zip(xs, ys):
        if min(x, y, nx - 1 - x, ny - 1 - y) <= W or np.hypot(x - cx, y - cy) < 12:
            continue
        st = im[y - W : y + W + 1, x - W : x + W + 1]
        pk = st.max() - bg
        if pk < 12 * noise or pk > 45000:
            continue
        stamps.append(dict(x=int(x), y=int(y), peak=pk, d=st))
    stamps.sort(key=lambda q: np.hypot(q["x"] - cx, q["y"] - cy))
    f = fieldshape(stamps[:28])
    if f is None:
        return {"id": file.stem, "status": "insufficient_psf"}
    shape, stars = f
    ix, iy = int(round(cx)), int(round(cy))
    st = im[iy - W : iy + W + 1, ix - W : ix + W + 1]
    cp = np.array(w.all_pix2world([[cx, cy], [cx + 1, cy], [cx, cy + 1]], 0))
    cd = (cp[1:] - cp[0]).T * 3600
    cd[0] *= np.cos(np.deg2rad(dec))
    scale = abs(np.linalg.det(cd)) ** 0.5

    def fit(st):
        p, improve, red = fit_source(st, shape, scale, noise)
        east, north = cd @ p[4:6]
        return dict(
            sep=float(np.hypot(east, north)),
            pa=float(np.rad2deg(np.arctan2(east, north)) % 180),
            gain=float(improve),
            reduced=float(red),
            ratio=float(p[6]),
            peak=float((st.max() - bg) / noise),
            east=float(east),
            north=float(north),
        )

    # Similar-intensity held-in comparison stars bound empirical spurious two-source fit gains.
    stars.sort(key=lambda z: abs(np.log(z["peak"] / (st.max() - bg))))
    out = fit(st)
    controls = [fit(v["d"]) for v in stars[:6]]
    return dict(
        id=file.stem,
        epoch=float(meta["MJD"]),
        exptime=int(meta["Exptime"]),
        pixel_scale=scale,
        filter=meta["Filter"],
        psf_stars=len(stars),
        **out,
        controls=controls,
    )


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor

    d = Path(__file__).parent
    fs = sys.argv[1:] or sorted(map(str, d.glob("*p.fits")))
    with ProcessPoolExecutor(max_workers=5) as pool:
        for r in pool.map(measure, fs):
            (d / (r["id"] + ".fit.json")).write_text(json.dumps(r, indent=2) + "\n")
            print(
                r["id"],
                [
                    round(r.get(k, 0), 3)
                    for k in ("epoch", "sep", "pa", "gain", "reduced", "peak")
                ],
                r.get("status"),
                flush=True,
            )
