"""Count Swift/XRT and Chandra/ACIS events in bands at the two optical positions."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
FIELDS = ("instrument", "obsid", "target", "band_keV", "on", "off", "scaled_off",
          "time_s", "gross_upper95_counts_per_ks")


def read(filename):
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(HERE / filename, "rt") as stream:
        return list(csv.DictReader(stream))


def count_row(instrument, obsid, name, band, on, off, alpha, exposure):
    # Exact one-sided Poisson upper bound on ALL aperture counts. Background
    # can only increase the mean, so the source-only bound is conservative.
    upper = 0.5 * chi2.ppf(0.95, 2 * (on + 1))
    return (instrument, obsid, name, band, on, off, f"{alpha * off:.3f}",
            f"{exposure:.3f}", f"{1000 * upper / exposure:.4f}")


def acis(targets):
    obs = {o["obsid"]: float(o["exposure_s"]) for o in read("observations.csv")
           if "/acisf" in o["url"] and o["date_obs"].startswith("2024")}
    events = read("photons.csv.gz")
    all_ids = "+".join(sorted(obs))
    for target in targets:
        name = target["name"]
        counts = [(np.hypot(float(e["x_offset_arcsec"]), float(e["y_offset_arcsec"])),
                   (int(e["pi"]) - 0.5) * 0.0146)
                  for e in events if e["target"] == name and e["obsid"] in obs]
        for low, high in ((0.3, 1), (1, 2), (0.3, 2), (2, 7), (0.3, 7)):
            on = sum(r <= 2 for r, energy in counts if low <= energy < high)
            off = sum(4 < r <= 8 for r, energy in counts if low <= energy < high)
            yield count_row("ACIS", all_ids, name, f"{low:g}-{high:g}", on, off, 1/12, sum(obs.values()))


def product(cache, url, expected):
    path = cache / url.rsplit("/", 1)[1]
    if not path.exists():
        with urlopen(url, timeout=90) as response:
            contents = response.read()
        if hashlib.sha256(contents).hexdigest() != expected:
            raise ValueError(f"Archive digest differs for {path.name}")
        path.write_bytes(contents)
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Cached digest differs for {path.name}")
    return path


def swift(cache, targets):
    from astropy.io import fits
    from astropy.wcs import WCS

    cache.mkdir(parents=True, exist_ok=True)
    for obs in read("swift-observations.csv"):
        event_file = product(cache, obs["event_url"], obs["event_sha256"])
        map_file = product(cache, obs["map_url"], obs["map_sha256"])
        with fits.open(event_file) as hdus, fits.open(map_file) as maps:
            data = hdus["EVENTS"].data
            header = hdus["EVENTS"].header
            image = maps[0].data
            mh = maps[0].header
            # Cleaned event X/Y are 1-based FITS sky pixels. Reject mismatched maps.
            for i, colname in enumerate(("X", "Y"), 1):
                col = hdus["EVENTS"].columns.names.index(colname) + 1
                for event_key, map_key in (("TCRPX", "CRPIX"), ("TCRVL", "CRVAL"), ("TCDLT", "CDELT")):
                    if not np.isclose(header[f"{event_key}{col}"], mh[f"{map_key}{i}"]):
                        raise ValueError("Event and exposure map coordinates differ")
            sky = WCS(mh)
            rows, cols = np.ogrid[1:image.shape[0]+1, 1:image.shape[1]+1]
            accepted = (data["GRADE"] >= 0) & (data["GRADE"] <= 12)
            for target in targets:
                x, y = sky.all_world2pix([[float(target["ra_deg"]), float(target["dec_deg"])]], 1)[0]
                radius = np.hypot(cols - x, rows - y)
                on_map = image[radius <= 6]
                off_map = image[(radius > 12) & (radius <= 24)]
                if on_map.mean() <= 0 or off_map.sum() <= 0:
                    raise ValueError("Insufficient exposure")
                alpha = on_map.sum(dtype=float) / off_map.sum(dtype=float)
                r = np.hypot(data["X"] - x, data["Y"] - y)
                for low, high in ((30, 100), (100, 200), (200, 1001), (30, 1001)):
                    sel = accepted & (data["PI"] >= low) & (data["PI"] < high)
                    on = np.count_nonzero(sel & (r <= 6))
                    off = np.count_nonzero(sel & (r > 12) & (r <= 24))
                    band = f"{low/100:g}-{min(high, 1000)/100:g}"
                    yield count_row("XRT-PC", obs["obsid"], target["name"], band,
                                    on, off, alpha, on_map.mean(dtype=float))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--swift-cache", type=Path, help="Include Swift, downloading missing FITS products")
    args = p.parse_args()
    targets = [t for t in read("targets.csv") if t["name"] in {"M31N 2017-06c", "M31N 2021-06b"}]
    import sys
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(FIELDS)
    writer.writerows(acis(targets))
    if args.swift_cache:
        writer.writerows(swift(args.swift_cache, targets))
