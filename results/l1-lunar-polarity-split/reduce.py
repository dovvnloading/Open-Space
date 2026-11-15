"""Rebuild one-minute event evidence from NASA SPDF HAPI CSVs (pandas, numpy)."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
FULL = ("2019-01-31T10:00:00Z", "2019-02-01T02:00:00Z")
LUNAR = ("2019-01-31T15:00:00Z", "2019-01-31T23:00:00Z")
GSE, GSM, POS = ["ex", "ey", "ez"], ["mx", "my", "mz"], ["x", "y", "z"]
DATA = {
    "ace16": ("AC_H0_MFI", "BGSEc,BGSM,SC_pos_GSE", GSE + GSM + POS, FULL),
    "wind1": ("WI_H0_MFI@0", "BGSM,BGSE,PGSE", GSM + GSE + POS, FULL),
    "dsc1s": ("DSCOVR_H0_MAG", "B1GSE,FLAG1", GSE + ["flag"], FULL),
    "dscpos": ("DSCOVR_ORBIT_PRE", "GSE_POS", POS, FULL),
    "aceplasma": ("AC_H0_SWE", "Np,Vp,Tpr,V_GSE", ["n", "v", "temp", "vx", "vy", "vz"], FULL),
    "omni": ("OMNI_HRO2_1MIN", "BZ_GSM,flow_speed,AE_INDEX,SYM_H", ["mz", "speed", "ae", "sym"], FULL),
    "omni_source": ("OMNI_HRO2_1MIN", "IMF", ["id"], LUNAR),
}
for sc in ["thb", "thc"]:
    DATA[sc] = (sc.upper() + "_L2_FGM@0", f"{sc}_fgm_fgs_quality,{sc}_fgs_gseQ,{sc}_fgs_gsmQ",
                ["flag"] + GSE + GSM, LUNAR)
    DATA[sc + "pos"] = (sc.upper() + "_OR_SSC", "XYZ_GSE", POS, LUNAR)
    DATA[sc + "esaQ"] = (sc.upper() + "_L2_ESA@0", f"{sc}_peif_data_quality,{sc}_peif_density,{sc}_peif_velocity_gse",
                         ["flag", "n", "vx", "vy", "vz"], LUNAR)


def reduce(cache, out):
    def read(name):
        dataset, parameters, columns, (start, end) = DATA[name]
        query = urllib.parse.urlencode(dict(id=dataset, parameters=parameters, format="csv",
                                            **{"time.min": start, "time.max": end}))
        url = "https://cdaweb.gsfc.nasa.gov/hapi/data?" + query
        path = cache / f"jan31_{name}.csv"
        if not path.exists():
            with urllib.request.urlopen(url, timeout=120) as response:
                content = response.read()
            if not content.startswith(b"2019-01-31T"):
                raise ValueError(f"Unexpected HAPI response from {dataset}")
            path.write_bytes(content)
        raw = path.read_bytes()
        table = pd.read_csv(path, header=None, names=["time"] + columns, index_col=0)
        table.index = pd.to_datetime(table.index, utc=True)
        if "flag" in table:
            table = table.loc[table.flag.eq(0)].drop(columns="flag")
        table = table.where(table.abs() < 1e8)
        if name == "omni":
            table["mz"] = table.mz.mask(table.mz.eq(9999.99))
        manifest = {"url": url, "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        return name, table.resample("min").mean(), manifest

    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        products = list(executor.map(read, DATA))
    tables = {name: data for name, data, _ in products}
    manifest = {name: meta for name, _, meta in products}

    # Earth-based GSE and GSM share +X. Wind publishes both vector representations.
    # Estimating their rotation uses Wind's basis, not its field as a proxy for DSCOVR.
    w = tables["wind1"]
    norm = w.ey**2 + w.ez**2
    sine = ((w.ey*w.mz - w.ez*w.my)/norm).rolling(11, center=True, min_periods=1).median()
    cosine = ((w.ey*w.my + w.ez*w.mz)/norm).rolling(11, center=True, min_periods=1).median()
    a = tables["ace16"]
    output = pd.DataFrame(index=w.index)
    output["ace_bz"] = a.mz
    output["wind_bz"] = w.mz
    ds = tables["dsc1s"]
    output["dscovr_bz"] = ds.ey*sine + ds.ez*cosine
    output["ace_rotation_residual"] = a.mz - a.ey*sine - a.ez*cosine
    for field in ["vx", "vy", "vz", "v"]:
        output["ace_" + field] = tables["aceplasma"][field]
    for sc, source, multiplier in [("ace", a, 1), ("wind", w, 6378),
                                   ("dscovr", tables["dscpos"], 1),
                                   ("thb", tables["thbpos"], 6378), ("thc", tables["thcpos"], 6378)]:
        for axis in POS:
            output[f"{sc}_{axis}_km"] = source[axis] * multiplier
    for sc in ["thb", "thc"]:
        output[f"{sc}_bz"] = tables[sc].mz
        for field in ["vx", "vy", "vz", "n"]:
            output[f"{sc}_{field}"] = tables[sc + "esaQ"][field]
    output["omni_bz"] = tables["omni"].mz
    output["omni_source"] = tables["omni_source"].id
    output.index.name = "utc"
    # Separate science series and orbit tables so each evidence artifact remains small.
    position = [column for column in output if column.endswith("_km")]
    for filename, frame in [("orbits.csv.gz", output[position]), ("minute.csv.gz", output.drop(columns=position))]:
        frame.to_csv(out / filename, float_format="%.1f" if filename == "orbits.csv.gz" else "%.5f",
                     compression={"method": "gzip", "mtime": 0})
    (out / "sources.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT)
    args = parser.parse_args()
    reduce(args.cache, args.output)
