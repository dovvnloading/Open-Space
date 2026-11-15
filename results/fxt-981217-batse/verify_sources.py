"""Verify frozen extracts against source bytes. Supply a local cache to avoid re-downloading."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np

BASE = Path(__file__).resolve().parent / "inputs"


def read_csv(filename):
    with (BASE / filename).open(newline="") as f:
        return list(csv.DictReader(f))


def binary_tables(data):
    """Return headers and data offsets from the specific no-heap BATSE FITS product."""
    tables = {}
    offset = 0
    while offset < len(data):
        pos = offset
        header = {}
        while data[pos:pos+8] != b"END     ":
            if pos + 80 > len(data):
                raise ValueError("Incomplete FITS header")
            card = data[pos:pos+80].decode("ascii")
            pos += 80
            if card[8:10] == "= ":
                header[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
        start = offset + ((pos-offset+80+2879)//2880)*2880
        size = int(header.get("NAXIS1",0))*int(header.get("NAXIS2",0))
        assert int(header.get("PCOUNT",0)) == 0, "unexpected FITS heap"
        tables[header.get("EXTNAME", "primary")] = (header, start)
        offset = start + ((size+2879)//2880)*2880
    return tables


def verify_batse(packed):
    data = gzip.decompress(packed)
    tables = binary_tables(data)
    h, off = tables["BATSE_CNTS"]
    assert h["TDIM7"] == "(8,6)" and h["TTYPE7"] == "COUNTS" and h["TTYPE9"] == "FLAGS"
    fmt = np.dtype([("mjd",">f8"),("spacecraft",">f4",5),("counts",">i2",(6,8)),
                    ("deadtime",">f4",8),("flags","u1",4),("attitude",">f4",4)])
    assert int(h["NAXIS1"]) == fmt.itemsize
    samples = np.frombuffer(data,fmt,int(h["NAXIS2"]),off)
    samples = samples[np.abs((samples["mjd"]-11164.079960)*86400)<850]
    local = read_csv("batse_lc.csv")
    assert len(local) == len(samples)
    for r, item in zip(local,samples,strict=True):
        assert abs(float(r["mjd"])-(item["mjd"]+40000)) < 1.e-10
        np.testing.assert_array_equal([int(r[f"d{d}_ch{c}"]) for c in range(1,7) for d in range(8)],
                                      item["counts"].ravel())
        np.testing.assert_array_equal([int(r[f"flag{i}"]) for i in range(1,5)], item["flags"])
    h, off = tables["BATSE_E_CALIB"]
    fmt = np.dtype([("det",">i2"),("time",">f8",2),("edges",">f4",5),("name","S16")])
    assert fmt.itemsize == int(h["NAXIS1"])
    energies = np.frombuffer(data, fmt, int(h["NAXIS2"]), off)
    for r, item in zip(read_csv("energy_edges.csv"),energies,strict=True):
        assert int(r["lad"]) == item["det"]
        np.testing.assert_allclose([float(r[f"edge{i}_keV"]) for i in range(1,6)],item["edges"],rtol=1e-7)


def verify_catalogs(stern, published):
    text = gzip.decompress(stern).decode()
    keys = re.search(r"^line\[1\] = (.+)",text,re.M)[1].split()
    original = [dict(zip(keys,r.split("|")[:-1],strict=True)) for r in text.split("<DATA>")[1].splitlines()
                if r.startswith("STK")]
    for saved, r in zip(read_csv("stern_catalog.csv"), original, strict=True):
        assert list(saved.values()) == [r[k] for k in
            ("name","time","trigger_num","peak_flux","ra","dec","error_radius","t90")]
    table = published.decode().split("#Table\tJ_A_A_711_A60_tablea2:")[1].split("#Table")[0]
    lines = table.splitlines()
    hdr = next(r for r in lines if r.startswith("recno\t")).split("\t")
    original = [dict(zip(hdr, [c.strip() for c in r.split("\t")],strict=True)) for r in lines
                if r.split("\t")[0].strip().isdigit()]
    assert original[49] == read_csv("published_wfc_row.csv")[0]
    for saved, r in zip(read_csv("wfc_catalog.csv"), original,strict=True):
        assert [saved[k].lower() for k in ("id","batse_id","class")] == [r[k].lower() for k in ("ID","BATSE","Type")]
        for s,k in [("mjd","MJDstart"),("ra","RAJ2000"),("dec","DEJ2000"),("r68_arcmin","ePos"),("t90_s","T90")]:
            if not saved[s] and not r[k]:
                continue
            # Preserve the numerical precision in each product independently.
            if abs(float(saved[s])-float(r[k])) > max(1.e-5, 0.0001*abs(float(r[k])) if s == "t90_s" else 0):
                raise ValueError(f"{saved['id']} {s}: {saved[s]} vs {r[k]}")


def verify_wfc_lightcurves(src):
    with tarfile.open(fileobj=io.BytesIO(src),mode="r:gz") as archive:
        directory = "WFC_FXT_OnlineData/"
        xray = np.loadtxt(archive.extractfile(directory+"WFC/981217_GRB_b10_LC"))
        saved = np.loadtxt(BASE/"wfc_lc.csv",delimiter=",",skiprows=1)
        np.testing.assert_array_equal(saved,xray[:,[2,3,7,8,4]])
        grbm = np.loadtxt(archive.extractfile(directory+"GRBM/981217.grbm"),skiprows=1)
        prior = grbm[grbm[:,0]<6908.544][-2:]
        after = grbm[grbm[:,0]>6908.544][:2]
        np.testing.assert_allclose(np.vstack((prior,after)), np.loadtxt(BASE/"grbm_nearest_samples.csv",delimiter=",",skiprows=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,help="optional directory for source files, named by provenance key")
    opt = parser.parse_args()
    manifest = json.loads((BASE/"provenance.json").read_text())
    source = {}
    for key in ["batse_lc","wfc_lc","stern_catalog","published_wfc_row"]:
        metadata = manifest[key]
        target = opt.cache/key if opt.cache else None
        if target and target.is_file():
            blob = target.read_bytes()
        else:
            with urlopen(metadata["url"],timeout=90) as response:
                blob = response.read()
            if target:
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes(blob)
        digest = hashlib.sha256(blob).hexdigest()
        if digest != metadata["sha256"] and key != "published_wfc_row":
            raise ValueError(f"Source version differs for {key}: {digest}")
        # VizieR responses carry volatile query metadata; compare all selected table values below.
        if digest != metadata["sha256"]:
            print(f"VizieR response hash changed: {digest}. Checking row values.")
        source[key]=blob
    verify_batse(source["batse_lc"])
    verify_wfc_lightcurves(source["wfc_lc"])
    verify_catalogs(source["stern_catalog"],source["published_wfc_row"])
    print("Verified BATSE light curve, calibration, WFC and GRBM light curves, and both complete screening populations.")
