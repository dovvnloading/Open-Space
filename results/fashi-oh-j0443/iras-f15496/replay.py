"""Replay the infrared-bright subset using the companion campaign's frozen inputs."""
import argparse
import gzip
import hashlib
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def xyz(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    path = args.cache/"heasarc.raw"
    if not path.exists():
        urllib.request.urlretrieve("https://heasarc.gsfc.nasa.gov/FTP/heasarc/dbase/tdat_files/heasarc_irasfsc.tdat.gz", path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != "c31afc6cec2e4130e8a696fa54a31bbf3dc09a72d153bd6a542e79f44762d289":
        raise ValueError("IRAS input changed")
    with gzip.open(path, "rt") as handle:
        columns = []
        for line in handle:
            if line.strip() == "<DATA>":
                break
            match = re.match(r"field\[([^\]]+)\]", line)
            if match:
                columns.append(match[1])
        iras = pd.read_csv(handle, sep="|", names=columns+["end"], usecols=["name","ra","dec"],
                           header=None).dropna(subset=["ra","dec"])
    fast = pd.read_csv(args.cache/"dr2.csv")
    galaxies = pd.concat([pd.read_csv(args.cache/f"desi{n}.csv", dtype={"targetid": str})
                          for n in range(17,28)], ignore_index=True)
    if len(fast)!=156411 or len(galaxies)!=2454349:
        raise ValueError("Unexpected FASHI or DESI input row count; run the companion replay")
    a = xyz(fast.ra, fast.dec)
    b = xyz(galaxies.ra, galaxies.dec)
    p = xyz(iras.ra, iras.dec)
    arcsec = 180 * 3600 / np.pi
    ir_dist, ir_idx = cKDTree(p).query(a)
    r = cKDTree(a).sparse_distance_matrix(cKDTree(b), 2*np.sin(10/arcsec),output_type="ndarray")
    i, j = r["i"], r["j"]
    sep = 2*arcsec*np.arcsin(r["v"]/2)
    score = sep/np.maximum(6,fast.radec_err.values[i]*3600)
    offset = 299792.458*np.log((1+galaxies.z.values[j])*fast.freq.values[i]/1667.35903)
    host_ir = 2*arcsec*np.arcsin(np.linalg.norm(b[j]-p[ir_idx[i]],axis=1)/2)
    keep = (score<2)&(ir_dist[i]*arcsec<60)&(host_ir<30)
    i,j,offset,sep,host_ir = [v[keep] for v in (i,j,offset,sep,host_ir)]
    args.output.mkdir(exist_ok=True,parents=True)
    controls = []
    for v in np.r_[0, np.arange(-6000,-599,600), np.arange(600,6001,600)]:
        accepted = abs(offset-v)<230
        controls.append((int(v),len(np.unique(i[accepted]))))
        if v==0:
            ri,sj=i[accepted],j[accepted]
            out = pd.DataFrame({"fashi_id":fast.fashi_id.values[ri],
                                "desi_targetid":galaxies.targetid.values[sj],
                                "iras_name":iras.name.values[ir_idx[ri]],
                                "optical_z":galaxies.z.values[sj],"sep_arcsec":sep[accepted],
                                "dv_kms":offset[accepted],"optical_iras_sep_arcsec":host_ir[accepted]})
            out.sort_values("fashi_id").to_csv(args.output/"infrared_screen.csv",index=False)
    pd.DataFrame(controls,columns=["shift_kms","radio_sources"]).to_csv(args.output/"controls.csv",index=False)


if __name__ == "__main__":
    main()
