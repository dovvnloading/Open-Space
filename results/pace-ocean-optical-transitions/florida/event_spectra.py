"""Retrieve native OCI bands for the eastern Florida May 2025 collocations."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
import numpy as np
import pandas as pd
from scipy.io import netcdf_file
from fetch import BASE, CACHE, download

ROOT = Path(__file__).resolve().parent


def main():
    info = pd.read_csv(ROOT / "usf_info.csv")
    bands = info.loc[info["Row Type"].eq("variable"), "Variable Name"]
    bands = [b for b in bands if b.startswith("Rrs_")]
    events = pd.read_csv(ROOT / "footprints.csv")
    events = events[
        (events.region == "florida_atlantic") & events.start.str.startswith("2025-05")
    ]
    receipts = []
    output = []
    for date, t in [("2025-05-13", 396), ("2025-05-17", 400), ("2025-05-21", 404)]:
        jobs = []
        for i in range(0, len(bands), 24):
            slab = f"[{t}][116:24:452][996:24:1092]"
            query = ",".join(v + slab for v in bands[i : i + 24])
            jobs.append((f"event_{t}_{i}.nc", BASE + quote(query, safe=",:")))
        with ThreadPoolExecutor(max_workers=3) as pool:
            chunks = list(pool.map(download, jobs))
        receipts.extend(chunks)
        for file, _, _, _ in chunks:
            with netcdf_file(CACHE / file, mmap=False) as ds:
                timestamp = pd.to_datetime(
                    ds.variables["time"].data[0], unit="s", utc=True
                )
                assert timestamp.strftime("%Y-%m-%d") == date
                lat, lon = [ds.variables[k].data for k in ("latitude", "longitude")]
                for (start, end), pairs in events.groupby(["start", "end"]):
                    if date not in (start, end):
                        continue
                    ys = np.abs(lat[:, None] - pairs.lat.to_numpy()).argmin(0)
                    xs = np.abs(lon[:, None] - pairs.lon.to_numpy()).argmin(0)
                    assert (
                        np.max(np.abs(lat[ys] - pairs.lat)) < 1e-5
                        and np.max(np.abs(lon[xs] - pairs.lon)) < 1e-5
                    )
                    for band in ds.variables:
                        if band.startswith("Rrs_"):
                            v = ds.variables[band].data[0, ys, xs]
                            valid = np.isfinite(v)
                            output.append(
                                (
                                    start,
                                    end,
                                    date,
                                    int(band[4:]),
                                    len(v),
                                    valid.sum(),
                                    np.nanmedian(v) if valid.any() else np.nan,
                                )
                            )
    cols = [
        "start",
        "end",
        "date",
        "wavelength_nm",
        "footprint_n",
        "valid_n",
        "median_Rrs_sr_inverse",
    ]
    pd.DataFrame(output, columns=cols).sort_values(cols[:4]).to_csv(
        ROOT / "event_spectra.csv", index=False, float_format="%.9g"
    )
    pd.DataFrame(receipts, columns=["file", "bytes", "sha256", "url"]).to_csv(
        ROOT / "event_receipts.csv", index=False
    )


if __name__ == "__main__":
    main()
