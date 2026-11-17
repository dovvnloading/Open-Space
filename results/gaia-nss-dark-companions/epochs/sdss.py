"""Match a candidate CSV to every SDSS DR19 Astra visit, retaining native quality flags."""

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
from catalog_io import RangeStream, read_header, table_dtype, target_ids

SDSS = "https://dr19.sdss.org/sas/dr19/spectro/astra/0.6.0/summary/mwmAllVisit-0.6.0.fits.gz"
FIELDS = [
    "source_id",
    "survey",
    "instrument",
    "time_jd",
    "time_scale",
    "rv_kms",
    "err_kms",
    "quality",
    "spectrum_pk",
    "release",
    "telescope",
    "snr",
    "flags",
    "pipeline_flags",
    "n_components",
    "gaia_median_rv_kms",
    "sdss_id",
    "apogee_id",
]


def decode(x):
    return x.decode().strip() if isinstance(x, bytes) else str(x)


def normalize(t, inst):
    telescope = decode(t["telescope"])
    flags = int(t["zwarning_flags"] if inst == "BOSS" else t["spectrum_flags"])
    pipe = int(t["doppler_flags"]) if inst == "APOGEE" else 0
    components = int(t["n_components"]) if inst == "APOGEE" else -1
    if inst == "BOSS":
        tai0, tai1 = int(t["tai_beg"]), int(t["tai_end"])
        if min(tai0, tai1) <= 0 or tai1 < tai0:
            time_jd = float("nan")
        else:
            time_jd = 2400000.5 + (tai0 + tai1) / 172800.0
        scale = "TAI"
        rv, err = float(t["xcsao_v_rad"]), float(t["xcsao_e_v_rad"])
        ok = float(t["snr"]) > 10 and 0 < err < 30 and flags == 0
    else:
        scale = "UTC"
        start = decode(t["date_obs"])
        try:
            dt = datetime.fromisoformat(start).replace(tzinfo=timezone.utc)
            time_jd = 2440587.5 + (dt.timestamp() + float(t["exptime"]) / 2) / 86400
        except (OverflowError, ValueError):
            time_jd = float("nan")
        rv, err = float(t["v_rad"]), float(t["e_v_rel"])
        ok = (
            float(t["snr"]) > 10
            and 0 < err < 10
            and flags == 0
            and pipe == 0
            and components <= 1
        )
    ok = bool(
        ok and np.isfinite(rv) and np.isfinite(err) and 2400000 < time_jd < 2500000
    )
    return dict(
        source_id=str(t["gaia_dr3_source_id"]),
        survey="SDSS DR19",
        instrument=inst + "-" + telescope,
        time_jd=round(time_jd, 8),
        time_scale=scale,
        rv_kms=rv,
        err_kms=err,
        quality=int(ok),
        spectrum_pk=str(t["spectrum_pk"]),
        release=decode(t["release"]),
        telescope=telescope,
        snr=float(t["snr"]),
        flags=flags,
        pipeline_flags=pipe,
        n_components=components,
        gaia_median_rv_kms=float(t["gaia_v_rad"]),
        sdss_id=str(t["sdss_id"]),
        apogee_id=decode(t["sdss4_apogee_id"]),
    )


def main(targets, out):
    ids = target_ids(targets)
    out.mkdir(parents=True, exist_ok=True)
    ranges = RangeStream(SDSS, 934784647)
    rows = []
    counts = {}
    with ranges, gzip.GzipFile(fileobj=io.BufferedReader(ranges)) as f:
        while (h := read_header(f)) is not None:
            if h["NAXIS"] == "0":
                continue
            dt = table_dtype(h)
            n = int(h["NAXIS2"])
            inst = "BOSS" if "xcsao_v_rad" in dt.names else "APOGEE"
            assert inst not in counts
            for start in range(0, n, 20000):
                nt = min(20000, n - start)
                b = f.read(nt * dt.itemsize)
                assert len(b) == nt * dt.itemsize
                t = np.frombuffer(b, dtype=dt)
                sel = t[np.isin(t["gaia_dr3_source_id"], ids)]
                rows.extend(normalize(x, inst) for x in sel)
            f.read(-n * dt.itemsize % 2880)
            counts[inst] = n
    assert ranges.seen == 934784647 and set(counts) == {"BOSS", "APOGEE"}
    rows.sort(key=lambda r: (int(r["source_id"]), r["time_jd"], r["spectrum_pk"]))
    seen = set()
    for r in sorted(rows, key=lambda r: -r["snr"]):
        k = (r["source_id"], r["instrument"], r["time_jd"])
        if k in seen:
            r["quality"] = 0
        elif r["quality"]:
            seen.add(k)
    with (out / "sdss_visits.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    (out / "sdss_depth.json").write_text(
        json.dumps(
            dict(
                url=SDSS,
                sha256=ranges.digest.hexdigest(),
                input_bytes=ranges.seen,
                target_count=len(ids),
                target_file_sha256=hashlib.sha256(
                    Path(targets).read_bytes()
                ).hexdigest(),
                n_scanned=counts,
                n_matched=len(rows),
                n_quality=sum(r["quality"] for r in rows),
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", type=Path)
    parser.add_argument("out", type=Path)
    a = parser.parse_args()
    main(a.targets, a.out)
