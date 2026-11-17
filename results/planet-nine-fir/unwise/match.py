"""Measured geometrical coverage and chance associations in one Phan-region unTimely tile.

Compares pre-specified midpoint associations against translated trajectories. It is a
catalog-level test and deliberately does not infer image completeness or exclusions.
"""

import argparse, csv, gzip, hashlib, json, sys, warnings
from pathlib import Path
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.spatial import cKDTree
import numpy as np
from catalogue import HERE, TILE, INDEX, family, sky, epoch

RADIUS = 3.0


def xyz(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return (
        np.c_[np.cos(a) * np.cos(d), np.sin(a) * np.cos(d), np.sin(d)]
        * 206264.80624709636
    )


def time(mjd):
    return epoch(Time(mjd, format="mjd", scale="utc").isot)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def run(cache):
    receipt = json.loads((HERE / "receipts.json").read_text())
    assert sha(cache / INDEX) == receipt["index_sha256"]
    idx = fits.getdata(cache / INDEX)
    idx = idx[idx["COADD_ID"] == TILE]
    states = family(json.loads((HERE.parent / "orbits/orbit_input.json").read_text()))
    n = len(states[0])
    catalogs = {}
    meta = []
    for info in receipt["catalogues"]:
        fn = cache / info["file"]
        assert sha(fn) == info["sha256"]
        h = fits.getheader(fn)
        c = fits.getdata(fn)
        k = int(c["band"][0])
        e = int(c["EPOCH"][0])
        r = idx[(idx["BAND"] == k) & (idx["EPOCH"] == e)][0]
        assert len(c) == r["N_DET"]
        w = WCS(h, naxis=2)
        xy = xyz(c["ra"], c["dec"])
        f = c["flux"]
        df = c["dflux"]
        sn = np.divide(f, df, out=np.zeros_like(f), where=df > 0)
        quality = (
            (sn >= 5)
            & (c["qf"] >= 0.8)
            & (c["fracflux"] >= 0.5)
            & (c["flags_unwise"] == 0)
            & (c["nm"] >= 3)
        )
        support = (sn >= 3) & (c["qf"] >= 0.6) & (c["nm"] >= 2)
        # keep edge detections: all coordinates tested against the true tile rectangle.
        keep = np.flatnonzero(quality)
        catalogs[k, e] = dict(
            w=w,
            xy=xy[keep],
            ra=c["ra"][keep].copy(),
            dec=c["dec"][keep].copy(),
            sn=sn[keep].copy(),
            flux=f[keep].copy(),
            dflux=df[keep].copy(),
            qf=c["qf"][keep].copy(),
            rchi=c["rchi2"][keep].copy(),
            fracflux=c["fracflux"][keep].copy(),
            flags=c["flags_info"][keep].copy(),
            detid=c["unwise_detid"][keep].copy(),
            qt=cKDTree(xy[keep]),
            support=cKDTree(xy[support]),
            row=r,
            ndet=len(c),
        )
    # Count spatial recurrence in independent visits in the SAME band, regardless of
    # scan direction. A match in >=2 other epochs marks a persistent counterpart.
    for (k, e), cat in catalogs.items():
        repeated = np.zeros(len(cat["xy"]), int)
        for (bk, be), other in catalogs.items():
            if bk == k and be != e:
                repeated += (
                    other["support"].query(cat["xy"], distance_upper_bound=RADIUS)[0]
                    < RADIUS
                )
        cat["repeat"] = repeated
        cat["nonstatic"] = np.flatnonzero(repeated < 2)
        cat["nst"] = cKDTree(cat["xy"][cat["nonstatic"]])
    # The controls preserve each orbit and cadence with a constant offset on sky.
    ang = np.arange(16) * 2 * np.pi / 16
    offsets = [(0, 0)] + [
        (r * np.cos(a), r * np.sin(a)) for r in [60.0, 120.0] for a in ang
    ]
    N = len(offsets)
    epochs = sorted(set(e for _, e in catalogs))
    E = len(epochs)
    hits = np.zeros((2, N, E, n), bool)
    cover = hits.copy()
    qhit = hits.copy()
    details = []
    swept = np.zeros((2, E, n), bool)
    out_rows = []
    for (k, e), cat in sorted(catalogs.items()):
        it = epochs.index(e)
        r = cat["row"]
        mid = float(r["MJDMEAN"])
        lo = float(r["MJDMIN"])
        hi = float(r["MJDMAX"])
        ra, dc = sky(states, time(mid))
        u = xyz(ra, dc)
        start = xyz(*sky(states, time(lo)))
        end = xyz(*sky(states, time(hi)))
        sweep = np.linalg.norm(end - start, axis=1)
        w = cat["w"]
        for j, (east, north) in enumerate(offsets):
            ra1 = ra + east / (3600 * np.cos(np.deg2rad(dc)))
            dc1 = dc + north / 3600
            q = xyz(ra1, dc1)
            x, y = w.world_to_pixel_values(ra1, dc1)
            ins = (x >= -0.5) & (x < 2047.5) & (y >= -0.5) & (y < 2047.5)
            cover[k - 1, j, it] = ins
            dclean, ids = cat["qt"].query(q, distance_upper_bound=RADIUS)
            qhit[k - 1, j, it] = ins & (dclean < RADIUS)
            dn, ind = cat["nst"].query(q, distance_upper_bound=RADIUS)
            good = ins & (dn < RADIUS)
            hits[k - 1, j, it] = good
            if j == 0:
                for tr in np.flatnonzero(good):
                    i = cat["nonstatic"][ind[tr]]
                    details.append(
                        [
                            tr,
                            k,
                            e,
                            mid,
                            cat["detid"][i],
                            cat["ra"][i],
                            cat["dec"][i],
                            dn[tr],
                            cat["flux"][i],
                            cat["dflux"][i],
                            cat["sn"][i],
                            cat["qf"][i],
                            cat["rchi"][i],
                            cat["fracflux"][i],
                            cat["flags"][i],
                            cat["repeat"][i],
                        ]
                    )
        # Relaxed diagnostic for paths traversed during a coadd. A 1-arcsec spacing
        # at maximum chord length permits at most 0.5 arcsec sampling discrepancy.
        for t in np.linspace(lo, hi, 2 + int(np.ceil(sweep.max() / 1.0))):
            at, dt = sky(states, time(t))
            q = xyz(at, dt)
            x, y = w.world_to_pixel_values(at, dt)
            ins = (x >= -0.5) & (x < 2047.5) & (y >= -0.5) & (y < 2047.5)
            swept[k - 1, it] |= ins & (
                cat["nst"].query(q, distance_upper_bound=RADIUS)[0] < RADIUS
            )
        out_rows.append(
            dict(
                band=k,
                epoch=e,
                mjd_min=lo,
                mjd_mid=mid,
                mjd_max=hi,
                mid_utc=Time(mid, format="mjd").isot,
                tile_frac_zero=float(r["FRAC_PIX_COV0"]),
                tile_min_N=int(r["COVMIN"]),
                catalog_rows=cat["ndet"],
                quality_rows=len(cat["xy"]),
                nonpersistent_quality_rows=len(cat["nonstatic"]),
                in_tile=int(cover[k - 1, 0, it].sum()),
                any_quality_match=int(qhit[k - 1, 0, it].sum()),
                nonpersistent_match=int(hits[k - 1, 0, it].sum()),
                swept_nonpersistent_match=int(swept[k - 1, it].sum()),
                coadd_sweep_arcsec_min=round(sweep.min(), 3),
                coadd_sweep_arcsec_max=round(sweep.max(), 3),
            )
        )
    summary = dict(
        tile=TILE,
        total_tracks=n,
        n_band_epochs=len(catalogs),
        distinct_epoch_labels=epochs,
        n_control_offsets=N - 1,
        midpoint_radius_arcsec=RADIUS,
        coadd_sweep_arcsec_min=min(x["coadd_sweep_arcsec_min"] for x in out_rows),
        coadd_sweep_arcsec_max=max(x["coadd_sweep_arcsec_max"] for x in out_rows),
    )
    perband = {}
    for label, k in [("W1", 0), ("W2", 1), ("either", 2)]:
        collapse = lambda z: z[k] if k < 2 else np.any(z, axis=0)
        qc = collapse(qhit).sum(1)
        nh = collapse(hits).sum(1)
        nc = collapse(cover).sum(1)
        sw = (swept[k] if k < 2 else np.any(swept, axis=0)).sum(0)
        perband[label] = dict(
            epoch_coverage_per_track=[int(nc[0].min()), int(nc[0].max())],
            track_epoch_pairs_in_tile=int(nc[0].sum()),
            pairs_near_quality=int(qc[0].sum()),
            pairs_nonpersistent=int(nh[0].sum()),
            max_distinct_epochs=int(nh[0].max()),
            max_control_distinct_epochs=int(nh[1:].max()),
            n_tracks_ge3=int((nh[0] >= 3).sum()),
            n_tracks_ge4=int((nh[0] >= 4).sum()),
            control_tracks_ge3_quantiles=np.quantile(
                (nh[1:] >= 3).sum(1), [0, 0.5, 1]
            ).tolist(),
            swept_pairs=int(sw.sum()),
            swept_tracks_ge3=int((sw >= 3).sum()),
            control_nonpersistent_pair_fraction_quantiles=np.quantile(
                nh[1:].sum(1) / nc[1:].sum(1), [0, 0.5, 1]
            ).tolist(),
        )
    summary["band_counts"] = perband
    with (HERE / "matches.csv.gz").open("wb") as raw:
        with gzip.open(raw, "wt", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(
                [
                    "track_id",
                    "band",
                    "epoch",
                    "mjd_mid",
                    "detid",
                    "ra",
                    "dec",
                    "match_arcsec",
                    "flux_vega_nmgy",
                    "dflux",
                    "snr",
                    "qf",
                    "rchi2",
                    "fracflux",
                    "flags_info",
                    "other_epoch_repeats",
                ]
            )
            wr.writerows(details)
    with (HERE / "epochs.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    (HERE / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cache", type=Path)
    run(p.parse_args().cache)
