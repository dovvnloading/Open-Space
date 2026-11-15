"""Reproduce the FASHI–DESI radio/optical search and its frequency controls."""
import argparse
import hashlib
import os
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ARCSEC = 180 * 3600 / np.pi
C = 299792.458
FAST = "https://nadc.china-vo.org/res/file_upload/download?id=1538184"
FAST_SHA = "9c612acac23776b99630f72330f818ca01b6e76517d0107a7872b056b703db4a"
TAP = "https://datalab.noirlab.edu/tap/sync"
SOURCE = 20260037137


def download(path, url, sha=None):
    if not path.exists():
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as temp:
            partial = Path(temp.name)
        try:
            urllib.request.urlretrieve(url, partial)
            if sha and hashlib.sha256(partial.read_bytes()).hexdigest() != sha:
                raise ValueError(f"Unexpected checksum for {url}")
            if not sha and partial.read_bytes()[:18] != b"targetid,ra,dec,z,":
                raise ValueError(f"Expected DESI CSV from {url}")
            os.replace(partial, path)
        finally:
            partial.unlink(missing_ok=True)
    if sha and hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        raise ValueError(f"Unexpected checksum for {path}")


def desi_query(sql, max_rows=500000):
    return TAP + "?" + urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
        "MAXREC": max_rows, "QUERY": sql})


def vectors(ra, dec):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    args.output.mkdir(exist_ok=True, parents=True)
    radio_path = args.cache/"dr2.csv"
    download(radio_path, FAST, FAST_SHA)
    radio = pd.read_csv(radio_path)
    if len(radio) != 156411:
        raise ValueError("Unexpected DR2 row count")
    optical_bins = []
    for index in range(17, 28):
        sql = (f"SELECT targetid,mean_fiber_ra as ra,mean_fiber_dec as dec,z,zerr "
               f"FROM desi_dr1.zpix WHERE z >= {index/100} AND z < {(index+1)/100} "
               "AND zwarn = 0 AND deltachi2 > 40 AND zcat_primary='true' AND spectype='GALAXY'")
        path = args.cache/f"desi{index}.csv"
        download(path, desi_query(sql))
        chunk = pd.read_csv(path, dtype={"targetid": str})
        if len(chunk) >= 500000:
            raise ValueError("DESI query hit maximum rows")
        optical_bins.append(chunk)
    optical = pd.concat(optical_bins, ignore_index=True)
    if len(optical) != 2454349 or optical.targetid.duplicated().any():
        raise ValueError("Unexpected DESI DR1 galaxy sample")
    radio_tree = cKDTree(vectors(radio.ra, radio.dec))
    optical_tree = cKDTree(vectors(optical.ra, optical.dec))
    pairs = radio_tree.sparse_distance_matrix(optical_tree, 2*np.sin(1.5/ARCSEC),
                                              output_type="ndarray")
    i, j = pairs["i"], pairs["j"]
    dv = C * np.log((1+optical.z.values[j])*radio.freq.values[i]/1667.35903)
    sep = 2 * ARCSEC * np.arcsin(pairs["v"]/2)
    controls = []
    for offset in np.r_[0, np.arange(-6000, -299, 300), np.arange(300, 6001, 300)]:
        hits, counts = np.unique(i[abs(dv-offset)<150], return_counts=True)
        controls.append((int(offset), len(hits), int(sum(counts>=2)), int(sum(counts))))
    pd.DataFrame(controls, columns=["shift_kms","hosts","double_optical_targets","target_rows"]
                 ).to_csv(args.output/"controls.csv", index=False)
    keep = abs(dv)<150
    out = optical.iloc[j[keep]].copy()
    out.columns = ["desi_targetid","spec_ra","spec_dec","spec_z","spec_zerr"]
    out.insert(0, "ID_FASHI", radio.fashi_id.values[i[keep]])
    out["sep_arcsec"] = sep[keep]
    out["v_optical_minus_oh_kms"] = dv[keep]
    out.sort_values("sep_arcsec").to_csv(args.output/"optical_matches.csv", index=False)
    radio[radio.fashi_id==SOURCE].to_csv(args.output/"fashi.csv", index=False)


if __name__ == "__main__":
    main()
