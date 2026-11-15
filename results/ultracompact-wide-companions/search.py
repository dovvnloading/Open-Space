"""Replay the Gaia search from the pinned Zenodo catalog. Python, astropy, numpy, pandas."""
import argparse
import hashlib
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
CATALOG = "https://zenodo.org/records/18842123/files/amcvn_catalogue.fits"
SERVERS = ["https://gea.esac.esa.int/tap-server/tap/sync", "https://gaia.ari.uni-heidelberg.de/tap/sync"]
COLUMNS = "source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ruwe,phot_g_mean_mag,bp_rp,radial_velocity,radial_velocity_error,duplicated_source,ipd_frac_multi_peak,visibility_periods_used,pmra_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr"
NEIGHBOR_LIMIT_AU = 100_000


def download(url, path):
    with urlopen(url, timeout=120) as response:
        data = response.read()
    tmp = path.with_suffix(".part")
    tmp.write_bytes(data)
    tmp.replace(path)


def query(adql, path):
    params = urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=adql, MAXREC=150000))
    for attempt in range(4):
        try:
            if not path.exists():
                download(SERVERS[attempt % len(SERVERS)] + "?" + params, path)
            rows = pd.read_csv(path, dtype={"source_id": str})
            if "source_id" not in rows or len(rows) >= 150000:
                raise ValueError(f"Invalid or possibly truncated response in {path}")
            return rows
        except (OSError, ValueError):
            path.unlink(missing_ok=True)
            if attempt == 3:
                raise
            time.sleep(attempt + 1)


def screen(host, cache):
    plx = host.parallax
    radius = min(0.5, max(60, 100 * plx) / 3600)
    adql = (f"SELECT {COLUMNS} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
            f"CIRCLE('ICRS',{host.ra},{host.dec},{radius:.9f})) AND parallax IS NOT NULL")
    path = cache / "cones" / (host.source_id + ".csv")
    neighbors = query(adql, path)
    a = neighbors.loc[neighbors.source_id != host.source_id].copy()
    dra, ddec = np.deg2rad(a.ra - host.ra), np.deg2rad(a.dec - host.dec)
    hav = np.sin(ddec/2)**2 + np.cos(np.deg2rad(host.dec)) * np.cos(np.deg2rad(a.dec)) * np.sin(dra/2)**2
    a["au"] = np.rad2deg(2 * np.arcsin(np.sqrt(hav))) * 3_600_000 / plx
    a["vt_kms"] = 4.74047 * np.hypot(a.pmra - host.pmra, a.pmdec - host.pmdec) / plx
    a["escape_kms"] = 42.122 * np.sqrt(2 / a.au)
    a["parallax_sigma"] = abs(a.parallax - plx) / np.hypot(np.hypot(a.parallax_error, host.parallax_error), .025)
    within = a.au <= NEIGHBOR_LIMIT_AU
    quality = ((plx >= 5 * host.parallax_error) & (a.parallax >= 5 * a.parallax_error)
               & (host.ruwe < 1.4) & (a.ruwe < 1.4) & (~a.duplicated_source) & (not host.duplicated_source))
    agreement = within & quality & (a.parallax_sigma <= 3)
    possible = agreement & (a.vt_kms <= a.escape_kms)
    receipt = dict(name=host.object_name, source_id=host.source_id, radius_deg=radius,
                   queried=len(a), within=int(within.sum()), quality=int((within & quality).sum()),
                   parallaxes=int(agreement.sum()), passing=int(possible.sum()))
    # Hash a stable ordering; TAP endpoints can reorder the same DR3 rows.
    receipt["id_sha256"] = hashlib.sha256("\n".join(sorted(neighbors.source_id)).encode()).hexdigest()
    a["speed_ratio"] = a.vt_kms / a.escape_kms
    a["host_name"] = host.object_name
    return receipt, a.loc[agreement, ["host_name", "source_id", "au", "vt_kms", "escape_kms", "speed_ratio"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.cache / "cones").mkdir(parents=True, exist_ok=True)
    f = args.cache / "catalog.fits"
    if not f.exists():
        download(CATALOG, f)
    if hashlib.md5(f.read_bytes()).hexdigest() != "f8979ea966d5230c02dd3573ba2659c2":
        raise ValueError("Zenodo catalog checksum mismatch")
    catalog = fits.getdata(f)
    names = {r["Gaia_ID"].split()[-1]: r["Name"] for r in catalog if r["Gaia_ID"]}
    ids = list(names)
    batches = [query(f"SELECT {COLUMNS} FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(ids[i:i+25])})",
                     args.cache / f"hosts{i:03}.csv") for i in range(0, len(ids), 25)]
    all_hosts = pd.concat(batches, ignore_index=True)
    hosts = all_hosts.loc[all_hosts.pmra.notna() & (all_hosts.parallax >= 2 * all_hosts.parallax_error)
                          & (all_hosts.parallax > 0)].copy()
    hosts["object_name"] = hosts.source_id.map(names)
    with ThreadPoolExecutor(6) as pool:
        results = list(pool.map(lambda row: screen(row, args.cache), hosts.itertuples()))
    counts = pd.DataFrame([r[0] for r in results]).sort_values("source_id")
    counts.to_csv(args.output / "fields.csv", index=False)
    matches = pd.concat([r[1] for r in results]).sort_values(["speed_ratio", "source_id"])
    matches.to_csv(args.output / "parallax_matches.csv", index=False, float_format="%.8g")
    print(len(catalog), len(names), len(hosts), counts[["queried","within","quality","parallaxes","passing"]].sum().to_dict())


if __name__ == "__main__":
    main()
