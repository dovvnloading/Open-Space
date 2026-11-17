"""Inspect both RV estimators and the CCF decomposition in selected DR19 mwmVisit files."""

import argparse, csv, hashlib, io, json, sys, urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "epochs"))
from catalog_io import read_header, table_dtype

SOURCES = {
    "4661290764764683776",
    "1830677037678583680",
    "3664684869697065984",
    "3112097229257687680",
    "902085660651149440",
    "2130638481326861824",
    "2077800002648034944",
    "4658457941738020736",
    "2224573745414713984",
    "5713565792393224960",
}


def fetch(star, cache):
    sid, sdss_id = star
    name = f"mwmVisit-0.6.0-{sdss_id}.fits"
    url = f"https://dr19.sdss.org/sas/dr19/spectro/astra/0.6.0/spectra/visit/{sdss_id[-4:-2]}/{sdss_id[-2:]}/{name}"
    path = cache / name
    if not path.exists():
        with urllib.request.urlopen(url, timeout=90) as resp:
            path.write_bytes(resp.read())
    raw = path.read_bytes()
    f = io.BytesIO(raw)
    visits = []
    while h := read_header(f):
        length = int(h.get("NAXIS1", 0)) * int(h.get("NAXIS2", 0))
        if h.get("EXTNAME") == "APOGEE" and length:
            data = np.frombuffer(
                raw, table_dtype(h), count=int(h["NAXIS2"]), offset=f.tell()
            )
            for d in data:
                dr = float(d["v_rad"])
                xc = float(d["xcorr_v_rad"])
                # Some arrays carry the numeric value of the FITS null sentinel; test explicitly.
                good = np.isfinite(d["flux"]) & np.isfinite(d["ivar"]) & (d["ivar"] > 0)
                ratio = (
                    float(d["ccfwhm"]) / float(d["autofwhm"])
                    if float(d["autofwhm"]) > 0
                    else np.nan
                )
                visits.append(
                    dict(
                        source_id=sid,
                        sdss_id=sdss_id,
                        spectrum_pk=int(d["spectrum_pk_id"]),
                        mjd=int(d["mjd"]),
                        snr=float(d["snr"]),
                        v_rad_kms=dr,
                        xcorr_v_rad_kms=xc,
                        xcorr_minus_doppler_kms=xc - dr,
                        ccfwhm=float(d["ccfwhm"]),
                        autofwhm=float(d["autofwhm"]),
                        width_ratio=ratio,
                        n_components=int(d["n_components"]),
                        spectrum_flags=int(d["spectrum_flags"]),
                        doppler_flags=int(d["doppler_flags"]),
                        positive_ivar_flux_pixels=int(good.sum()),
                    )
                )
        f.seek((length + 2879) // 2880 * 2880, 1)
    return dict(
        source_id=sid, url=url, sha256=hashlib.sha256(raw).hexdigest(), nbytes=len(raw)
    ), visits


def main(cache, out):
    cache.mkdir(exist_ok=True, parents=True)
    root = Path(__file__).resolve().parents[1]
    selected = {}
    for p in (
        root / "epochs/sdss_visits.csv",
        root / "extension/epochs/sdss_visits.csv",
    ):
        with p.open() as f:
            for r in csv.DictReader(f):
                if r["source_id"] in SOURCES and r["instrument"].startswith("APOGEE"):
                    selected[r["source_id"]] = r["sdss_id"]
    with ThreadPoolExecutor(4) as e:
        results = list(e.map(lambda v: fetch(v, cache), sorted(selected.items())))
    records = []
    rows = []
    for rec, vis in results:
        records.append(rec)
        rows.extend(vis)
    rows.sort(key=lambda x: (x["source_id"], x["mjd"], x["spectrum_pk"]))
    (out / "apogee_spectrum_inputs.json").write_text(
        json.dumps(records, indent=2) + "\n"
    )
    with (out / "apogee_ccf.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    main(a.cache, a.out)
