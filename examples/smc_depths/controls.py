"""Repeat the eclipse-depth comparison on nearby binaries and older data."""

import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from followup import depth_measure, fit_profile, phase


def fetch_ogle_two(field, number, path):
    base = "https://ogledb.astrouw.edu.pl/~ogle/photdb/"
    fields = {
        "db_target": "smc",
        "dbtyp": "dia2",
        "use_field": "on",
        "val_field": field,
        "use_starid": "on",
        "val_starid": number,
        "disp_ra": "on",
        "disp_decl": "on",
        "disp_imean": "on",
        "sort": "field",
        "sorting": "ASC",
        "pagelen": "50",
    }
    with urlopen(
        base + "query.php?qtype=phot&first=1", urlencode(fields).encode(), timeout=60
    ) as f:
        listing = f.read().decode()
    detail = re.search(r"href='(getobj.php\?[^']+)'", listing).group(1)
    with urlopen(base + detail, timeout=60) as f:
        link = re.search(
            r"""href=['"](data/[^'"]+\.dat)['"]""", f.read().decode()
        ).group(1)
    with urlopen(base + link, timeout=60) as f:
        path.write_bytes(f.read())


def historical(raw, cache, catalog):
    results = {}
    for number, field, source in [
        ("07730", "SMC_SC11", "43067"),
        ("06703", "SMC_SC5", "27872"),
        ("04349", None, "4349"),
    ]:
        star = "OGLE-SMC-ECL-" + number
        if field:
            path = cache / (field + "-" + source + ".dat")
            if not path.exists():
                fetch_ogle_two(field, source, path)
        else:
            path = cache / "OGLE-SMC-ECL-4349.dat"
            url = (
                "https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/phot/I/"
                + path.name
            )
            if not path.exists():
                with urlopen(url, timeout=60) as f:
                    path.write_bytes(f.read())
        t, m, e = np.loadtxt(path, usecols=(0, 1, 2)).T
        if not field:
            t += 2450000.0
        good = (e > 0) & (e <= 0.1)
        t, m, e = t[good], m[good], e[good]
        per, ep = map(float, catalog[star][3:5])
        p = phase(t, ep, per)
        rng = np.random.default_rng(4)
        results[star] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "file": path.name,
            "usable": len(t),
            "hjd_range": [min(t), max(t)],
            "primary": depth_measure(t, p, m, 0, rng),
            "secondary": depth_measure(t, p, m, 0.5, rng),
            "profile_parameters": fit_profile(p, m, e).tolist(),
        }
    return results


def local_controls(raw, catalog):
    source = "OGLE-SMC-ECL-07730"
    ident = {line[:18]: line for line in (raw / "ident.dat").read_text().splitlines()}
    target_i, _, target_p = map(float, catalog[source][1:4])
    chip = ".".join(ident[source][49:65].split(".")[:2])
    rows = []
    for star, row in catalog.items():
        if (
            star == source
            or ident[star][20:23].strip() != "NC"
            or "-" in (row[1], row[3], row[4])
        ):
            continue
        if ".".join(ident[star][49:65].split(".")[:2]) != chip:
            continue
        i, per, ep = float(row[1]), float(row[3]), float(row[4])
        if abs(i - target_i) > 0.75 or abs(np.log2(per / target_p)) > 1:
            continue
        path = raw / "phot/I" / (star + ".dat")
        if not path.exists():
            continue
        t, m, e = np.loadtxt(path).T
        p = phase(t, ep, per)
        core = np.minimum(abs(p), 0.5 - abs(p)) < 0.035
        quad = (abs(p) > 0.16) & (abs(p) < 0.34)
        depths = []
        for era in [t < 2457000, t >= 2459700]:
            good = era & (e > 0) & (e <= 0.1)
            if min((good & core).sum(), (good & quad).sum()) < 12:
                break
            depths.append(np.median(m[good & core]) - np.median(m[good & quad]))
        if len(depths) == 2:
            rows.append(
                {
                    "id": star,
                    "early_mag": depths[0],
                    "late_mag": depths[1],
                    "change": depths[1] - depths[0],
                }
            )
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    a = parser.parse_args()
    a.cache.mkdir(parents=True, exist_ok=True)
    catalog = {
        r.split()[0]: r.split() for r in (a.raw / "ecl.dat").read_text().splitlines()
    }
    a.output.write_text(
        json.dumps(
            {
                "historical": historical(a.raw, a.cache, catalog),
                "same_chip_controls": local_controls(a.raw, catalog),
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
