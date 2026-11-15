"""Build the compact measurements from public mission products."""
import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
BAT = "https://swift.gsfc.nasa.gov/results/transients/"
MAXI = "https://maxi.riken.jp/star_data/"
PRODUCTS = {
    "smcbat": BAT + "SMCX-1.lc.txt",
    "smcmaxi": MAXI + "J0117-734/J0117-734_g_lc_1day_all.dat",
    "lmcbat": BAT + "LMCX-4.lc.txt",
    "lmcmaxi": MAXI + "J0532-663/J0532-663_g_lc_1day_all.dat",
    "smcfermi": "https://gammaray.nsstc.nasa.gov/gbm/science/earth_occ/lightcurves/SMCX-1_GBM_occ_fluxhist_auto.txt",
    "smc_bat_orb": BAT + "SMCX-1.orbit.lc.txt",
    "smc_maxi_orb": MAXI + "J0117-734/J0117-734_g_lc_1orb_all.dat",
}
START, STOP = 59500, 61349


def digest(b):
    return hashlib.sha256(b).hexdigest()


def extract(raw_dir, dest):
    dest.mkdir(parents=True, exist_ok=True)
    retrieved_utc = datetime.now(timezone.utc).date().isoformat() if raw_dir is None else None
    receipts = {}
    with (dest / "measurements.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("series", "mjd", "rate", "error", "n_measurements"))
        for name, url in PRODUCTS.items():
            if raw_dir:
                raw = (raw_dir / name).read_bytes()
            else:
                with urlopen(url, timeout=120) as response:
                    raw = response.read()
            if name == "smcfermi":
                lines = [line for line in raw.splitlines() if line[:5].isdigit()]
                a = np.loadtxt(lines)
                t = (a[:, 0] + a[:, 1]) / 2
                rate, err = a[:, 2:4].T
            else:
                a = np.loadtxt(io.BytesIO(raw))
                t, rate, err = a[:, :3].T
            good = np.isfinite(rate) & np.isfinite(err) & (err > 0)
            if "bat" in name:
                # BAT daily dates identify the UTC day. For phase masking only,
                # pointing MET is sufficient to far better than the eclipse cut.
                if name.endswith("_orb"):
                    t = 51910 + t / 86400
                    good &= (a[:, 10] == 0) & (a[:, 6] > 0)
                else:
                    t = t + .5
                    good &= a[:, 7] == 0
            ceiling = 4 * np.median(err[good])
            good &= err < ceiling
            keep = good & (t < STOP)
            if name != "smcbat":
                keep &= t >= START
            receipts[name] = dict(url=url, sha256=digest(raw), bytes=len(raw),
                                  source_rows=len(a), error_ceiling=float(ceiling),
                                  eligible_raw_rows=int(keep.sum()))
            t, rate, err = t[keep], rate[keep], err[keep]
            if not name.endswith("_orb"):
                for date, y, e in zip(t, rate, err):
                    writer.writerow((name, date, y, e, 1))
                continue
            dt = t - 60505.7361918
            phase = np.abs((.25696017 * dt + .5 * 2.379e-9 * dt**2 + .5) % 1 - .5)
            for label, phase_mask in (("outside", phase >= .17), ("core", phase <= .055)):
                day = t[phase_mask].astype(int)
                weight = 1 / err[phase_mask]**2
                total = np.bincount(day, weight)
                weighted_rate = np.bincount(day, weight * rate[phase_mask])
                counts = np.bincount(day)
                for i in np.flatnonzero(total):
                    writer.writerow((name + "_" + label, i + .5,
                                     f"{weighted_rate[i]/total[i]:.10g}",
                                     f"{1/np.sqrt(total[i]):.10g}", counts[i]))
    manifest = dict(retrieved_utc=retrieved_utc, products=receipts,
                    measurements_sha256=digest((dest / "measurements.csv").read_bytes()))
    (dest / "sources.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", type=Path, help="Read exact mission files named by the keys in sources.json")
    p.add_argument("--output", type=Path, required=True, help="Destination for the regenerated CSV and manifest")
    args = p.parse_args()
    extract(args.raw_dir, args.output)
