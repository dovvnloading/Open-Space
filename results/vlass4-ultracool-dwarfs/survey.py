"""Measure VLASS quick-look Stokes I at UltracoolSheet positions."""

import argparse
import hashlib
import io
import math
import warnings
import time
from urllib.error import HTTPError, URLError
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning

CADC = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca"
REF_EPOCH = 60859.0  # 2025 July 3 UTC
warnings.filterwarnings("ignore", category=AstropyWarning)


def xyz(ra, dec):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack([np.cos(r) * np.cos(d), np.sin(r) * np.cos(d), np.sin(d)])


def target_list(root, exclude_uris=()):
    objects = pd.read_csv(root / "ucd.csv", low_memory=False)
    tiles = pd.read_csv(root / "tiles.csv").merge(
        pd.read_csv(root / "files.csv").rename(columns={"observationID": "obs_id"}),
        on="obs_id",
    )
    tiles = (
        tiles[tiles.uri.str.endswith("tt0.subim.fits")]
        .dropna(subset=["s_ra", "s_dec", "s_region"])
        .reset_index(drop=True)
    )
    tiles = tiles[~tiles.uri.isin(exclude_uris)].reset_index(drop=True)
    good = objects.sptnum_formula.ge(6) & (
        objects.spt_ir.notna() | objects.spt_opt.notna()
    )
    good &= ~objects.literature_flag.fillna("").str.contains("photSpT", case=False)
    cols = [
        "ra_20250703_formula",
        "dec_20250703_formula",
        "pmra_formula",
        "pmdec_formula",
    ]
    selected = objects[good].dropna(subset=cols)
    tree = cKDTree(xyz(tiles.s_ra, tiles.s_dec))
    ra, dec = selected[cols[0]], selected[cols[1]]
    _, near = tree.query(xyz(ra, dec), k=6)
    rows = []
    for (idx, obj), indices in zip(selected.iterrows(), near):
        r, d = obj[cols[0]], obj[cols[1]]
        # Work in a tangent plane about the target. The cutout is tested
        # against the actual data pixels after nominal footprint selection.
        for tile in indices:
            t = tiles.iloc[tile]
            coords = np.array(t.s_region.split()[2:], float).reshape(-1, 2)
            x = (coords[:, 0] - r + 180) % 360 - 180
            y = coords[:, 1] - d
            inside = False
            for j in range(len(x)):
                k = (j + 1) % len(x)
                if (y[j] > 0) != (y[k] > 0) and x[j] - y[j] * (x[k] - x[j]) / (
                    y[k] - y[j]
                ) > 0:
                    inside = not inside
            if inside:
                rows.append(
                    dict(
                        id=idx,
                        name=obj["name"],
                        ra=r,
                        dec=d,
                        pmra=obj[cols[2]],
                        pmdec=obj[cols[3]],
                        sp=obj.sptnum_formula,
                        dist=obj.dist_formula,
                        uri=t.uri,
                    )
                )
                break
    return pd.DataFrame(rows).sort_values(["dist", "name"])


def cutout_url(uri, ra, dec):
    return (
        CADC
        + "/caom2ops/sync?"
        + urlencode({"ID": uri, "CIRCLE": f"{ra:.7f} {dec:.7f} 0.0125000"})
    )


def measure(data, row):
    with fits.open(io.BytesIO(data), memmap=False) as hdus:
        header = hdus[0].header
        wcs = WCS(header).celestial
        flux = np.squeeze(hdus[0].data) * 1000  # Jy/beam to mJy/beam
    mjd = Time(header["DATE-OBS"], scale="utc").mjd
    dt = (mjd - REF_EPOCH) / 365.25
    ra = row.ra + dt * row.pmra / 3600000 / math.cos(math.radians(row.dec))
    dec = row.dec + dt * row.pmdec / 3600000
    x, y = wcs.world_to_pixel_values(ra, dec)
    yy, xx = np.indices(flux.shape)
    rr = np.hypot(xx - x, yy - y)  # VLASS quick-look pixels are 1 arcsec
    if (
        abs(abs(header["CDELT1"]) * 3600 - 1) > 0.01
        or abs(abs(header["CDELT2"]) * 3600 - 1) > 0.01
    ):
        raise ValueError("unexpected pixel scale")
    bg = flux[(rr >= 8) & (rr <= 35) & np.isfinite(flux)]
    if len(bg) < 2000 or np.isfinite(flux[rr < 3]).sum() < 20:
        raise ValueError("incomplete cutout")
    med = np.median(bg)
    noise = 1.482602 * np.median(np.abs(bg - med))

    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("invalid background noise")

    def peak(x0, y0):
        pixels = flux[np.hypot(xx - x0, yy - y0) <= 1.5]
        if not np.isfinite(pixels).any():
            return np.nan
        return float(np.nanmax(pixels) - med)

    peakflux = peak(x, y)
    controls = [
        peak(x + 20 * math.cos(i * math.pi / 4), y + 20 * math.sin(i * math.pi / 4))
        / noise
        for i in range(8)
    ]
    return dict(
        date=header["DATE-OBS"],
        ra_obs=ra,
        dec_obs=dec,
        sha256=hashlib.sha256(data).hexdigest(),
        peak_mjy=round(peakflux, 5),
        rms_mjy=round(float(noise), 5),
        snr=round(float(peakflux / noise), 3),
        **{f"control{i}": round(float(v), 3) for i, v in enumerate(controls)},
    )


def fetch_measure(row, cache):
    path = cache / f"{row.id}.fits"
    try:
        if not path.exists():
            for attempt in range(2):
                try:
                    with urlopen(
                        cutout_url(row.uri, row.ra, row.dec), timeout=60
                    ) as response:
                        data = response.read()
                    break
                except (HTTPError, URLError, TimeoutError) as err:
                    if attempt == 1 or (
                        isinstance(err, HTTPError)
                        and err.code not in {429, 500, 502, 503, 504}
                    ):
                        from direct import extract

                        data = extract(row.uri, row.ra, row.dec)
                        break
                    time.sleep(2**attempt)
            if not data.startswith(b"SIMPLE"):
                raise ValueError("archive did not return FITS")
            path.write_bytes(data)
        return {"id": row.id, **measure(path.read_bytes(), row)}
    except Exception as err:
        return {"id": row.id, "error": str(err)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    cache = args.cache
    (cache / "fits").mkdir(parents=True, exist_ok=True)
    targets = target_list(cache)
    # Resume a successful overlap choice made in a previous run.
    receipt = cache / "overlap-tiles.csv"
    if receipt.exists():
        replacements = pd.read_csv(receipt).set_index("id").uri
        targets["uri"] = targets.id.map(replacements).fillna(targets.uri)
    targets.to_csv(cache / "targets.csv", index=False)
    print(f"{len(targets)} targets", flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(fetch_measure, r, cache / "fits") for r in targets.itertuples()
        ]
        result = []
        for i, fut in enumerate(futures, 1):
            result.append(fut.result())
            if i % 50 == 0 or i == len(futures):
                frame = targets.merge(pd.DataFrame(result), on="id")
                frame.to_csv(cache / "screen.csv", index=False)
                print(
                    i,
                    "rows",
                    frame.get("snr", pd.Series(dtype=float)).ge(6).sum(),
                    "above 6",
                    flush=True,
                )

    # A catalog footprint can extend beyond its finite science pixels.
    # Try a second covering image only for targets with insufficient coverage.
    replacements = []
    for failed in frame[
        frame.get("error", pd.Series("", index=frame.index)).eq("incomplete cutout")
    ].itertuples():
        alternatives = target_list(cache, exclude_uris=[failed.uri])
        for alt in alternatives[alternatives.id.eq(failed.id)].itertuples():
            retry_dir = cache / "overlap"
            retry_dir.mkdir(exist_ok=True)
            measured = fetch_measure(alt, retry_dir)
            if "snr" in measured:
                (retry_dir / f"{alt.id}.fits").replace(
                    cache / "fits" / f"{alt.id}.fits"
                )
                targets.loc[targets.id.eq(alt.id), "uri"] = alt.uri
                for k, v in {**measured, "uri": alt.uri, "error": ""}.items():
                    frame.loc[frame.id.eq(alt.id), k] = v
                replacements.append({"id": alt.id, "uri": alt.uri})
    if replacements:
        receipt_rows = pd.DataFrame(replacements)
        if receipt.exists():
            receipt_rows = pd.concat(
                [pd.read_csv(receipt), receipt_rows]
            ).drop_duplicates("id", keep="last")
        receipt_rows.to_csv(receipt, index=False)
        targets.to_csv(cache / "targets.csv", index=False)
        frame.to_csv(cache / "screen.csv", index=False)


if __name__ == "__main__":
    main()
