from pathlib import Path
import pandas as pd, numpy as np
from scipy.spatial import cKDTree
from catalog_screen import RA, DEC

bands = {}
out = []
rec = []
for b in ["277", "356"]:
    root = Path("blind-o1-f" + b + "-west")
    a = pd.read_csv(root / "detections.csv")
    c = pd.read_csv(root / "coverage.csv")
    z = pd.read_csv(root / "injections.csv")
    bands[b] = a.query("screen_pass==1").copy()
    st = np.load(root / "model.npz")
    for (mag, sign), t in z.groupby(["mag_ab", "sign"]):
        n = t.trials.sum()
        k = t.recovered.sum()
        p = k / n
        v = 1.96
        den = 1 + v * v / n
        err = v * np.sqrt(p * (1 - p) / n + v * v / (4 * n * n))
        mid = p + v * v / (2 * n)
        rec.append(
            dict(
                filter="F" + b + "W",
                mag_ab=mag,
                sign=sign,
                trials=n,
                recovered=k,
                efficiency=p,
                wilson95_low=(mid - err) / den,
                wilson95_high=(mid + err) / den,
            )
        )
    out.append(
        dict(
            filter="F" + b + "W",
            pixels=int(c.n_pix.sum()),
            area_arcmin2=c.n_pix.sum() * 9.30162342221895e-14 * (180 * 60 / np.pi) ** 2,
            peaks_7sigma=len(a),
            morphology_pass=len(bands[b]),
            median_tile_fluxerr_njy=c.median_fluxerr_njy.median(),
        )
    )
pd.DataFrame(rec).to_csv("blind_recovery.csv", index=False)
pd.DataFrame(out).to_csv("blind_area.csv", index=False)
a, b = bands["277"], bands["356"]
fn = lambda t: np.c_[(t.ra_deg - RA) * np.cos(np.deg2rad(DEC)), t.dec_deg - DEC] * 3600
dt, j = cKDTree(fn(a)).query(fn(b))
_, k = cKDTree(fn(b)).query(fn(a))
g = (dt < 0.12) & (a.snr.values[j] * b.snr.values > 0) & (k[j] == np.arange(len(b)))
c = b[g].copy()
z = a.iloc[j[g]]
d = pd.DataFrame(
    dict(
        site=np.arange(sum(g)) + 2000,
        ra=(c.ra_deg.values + z.ra_deg.values) / 2,
        dec=(c.dec_deg.values + z.dec_deg.values) / 2,
        sep_arcsec=dt[g],
        z277=z.snr.values,
        z356=c.snr.values,
        flux277_njy=z.dflux_njy.values,
        flux356_njy=c.dflux_njy.values,
    )
)
d.to_csv("blind_pairs.csv", index=False)
print(pd.DataFrame(out).round(3).to_string(index=False), "pairs", len(d))
