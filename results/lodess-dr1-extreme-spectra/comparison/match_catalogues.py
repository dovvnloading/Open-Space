"""Crossmatch the released LoDeSS positions against the native HEASARC and ASTRON tables.

Distances use ICRS great circles. Counters describe geometrical ambiguity.
An empty match or a coverage gap is not a flux limit. No resolution or scale
correction is applied. Statistical errors are described in SURVEYS.md.
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree

BANDS = [
    76,
    84,
    92,
    99,
    107,
    115,
    122,
    130,
    143,
    151,
    158,
    166,
    174,
    181,
    189,
    197,
    204,
    212,
    220,
    227,
]
# match radius (arcsec), integrated or peak intensity column, error or rms, other QA.
SURVEYS = {
    "gmrtas150m": (
        25,
        "int_flux_150_mhz",
        "int_flux_150_mhz_error",
        ["major_axis", "minor_axis", "structure_code", "flux_150_mhz", "rms_150_mhz"],
    ),
    "vlssr": (
        45,
        "flux_74_mhz",
        "flux_74_mhz_error",
        ["major_axis", "minor_axis", "rms_74_mhz"],
    ),
    "nvss": (
        30,
        "flux_20_cm",
        "flux_20_cm_error",
        ["major_axis", "minor_axis", "residual_code"],
    ),
    "rrs8c38mhz": (
        105,
        "int_flux_38_mhz",
        None,
        ["snr", "source_flag", "source_extent", "flux_38_mhz"],
    ),
    "wenss": (
        45,
        "int_flux_92_cm",
        "flux_92_cm_error",
        ["source_type", "fit_problem", "flux_92_cm"],
    ),
    "gleamegcat": (
        60,
        "int_flux_wide",
        "int_flux_wide_error",
        ["abs_flux_pct_error", "fit_flux_pct_error", "bii", "a_psf_wide", "b_psf_wide"],
    ),
}


def unit_xyz(ra, dec):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def separation(xyz, other):
    return (
        np.rad2deg(
            2 * np.arcsin(np.minimum(1, np.linalg.norm(xyz - other, axis=1) / 2))
        )
        * 3600
    )


def read_tdat(path, columns, min_dec):
    with gzip.open(path, "rt") as f:
        schema, rows = [], None
        for line in f:
            if line.startswith("#    TOTAL ROWS:"):
                rows = int(line.split(":")[-1].strip())
            if line.startswith("line[1] ="):
                schema = line.split(" = ")[1].split()
            if line.rstrip() == "<DATA>":
                break
        if not schema or rows is None:
            raise ValueError(f"Invalid HEASARC table: {path}")
        # Column names in the header are definitive; the TDAT ends in a trailer.
        result = []
        count = 0
        for chunk in pd.read_csv(
            f,
            sep="|",
            names=schema + ["_end"],
            usecols=columns,
            chunksize=35000,
            low_memory=False,
        ):
            entries = chunk.loc[chunk["dec"].notna()].copy()
            count += len(entries)
            result.append(entries[entries.dec > min_dec])
        if count != rows:
            raise ValueError(f"{path}: read {count} rows, expected {rows}")
        return pd.concat(result, ignore_index=True)


def survey_match(name, rawdir, src, out):
    radius, flux, err, extra = SURVEYS[name]
    cols = ["name", "ra", "dec", flux] + ([err] if err else []) + extra
    if name == "gleamegcat":
        for band in BANDS:
            cols.extend(
                [
                    f"int_flux_{band}_mhz",
                    f"int_flux_{band}_mhz_error",
                    f"local_rms_{band}_mhz",
                ]
            )
    table = read_tdat(rawdir / f"{name}.tdat.gz", cols, src.dec.min() - 1)
    positions = unit_xyz(table.ra, table.dec)
    tree = cKDTree(positions)
    targets = unit_xyz(src.ra, src.dec)
    distances, which = tree.query(targets, k=2)
    arcsec = np.rad2deg(2 * np.arcsin(distances / 2)) * 3600
    matches = arcsec[:, 0] <= radius
    candidates = table.iloc[which[matches, 0]].copy()
    candidates.index = src.index[matches]
    nearest = pd.DataFrame(index=src.index)
    nearest["n_match"] = tree.query_ball_point(
        targets, 2 * np.sin(np.deg2rad(radius / 3600) / 2), return_length=True
    )
    nearest.loc[matches, "sep_as"] = arcsec[matches, 0]
    nearest.loc[matches, "next_sep_as"] = arcsec[matches, 1]
    for col in ("name", "ra", "dec", *extra):
        nearest[col] = candidates[col]
    nearest["flux_jy"] = candidates[flux] / 1000
    if err:
        nearest["err_jy"] = candidates[err] / 1000
    else:
        # 8C S/N is peak/rms; for unresolved sources peak/SNR estimates local rms.
        nearest["err_jy"] = candidates["flux_38_mhz"] / candidates["snr"] / 1000
    if name == "gleamegcat":
        phot = []
        for band in BANDS:
            phot.append(
                pd.DataFrame(
                    {
                        "row": candidates.index,
                        "freq_mhz": band,
                        "flux_jy": candidates[f"int_flux_{band}_mhz"] / 1000,
                        "fit_err_jy": candidates[f"int_flux_{band}_mhz_error"] / 1000,
                        "rms_jy_beam": candidates[f"local_rms_{band}_mhz"] / 1000,
                        "external_scale_pct": candidates.abs_flux_pct_error,
                        "internal_scale_pct": candidates.fit_flux_pct_error,
                    }
                )
            )
        pd.concat(phot, ignore_index=True).to_csv(
            out / "gleam_subbands.csv.gz", index=False
        )
    nearest.add_prefix(name + "_").to_csv(
        out / f"match_{name}.csv.gz", index_label="row"
    )
    print(name, int(matches.sum()), flush=True)


def lotss_match(rawdir, src, out):
    parts = sorted(rawdir.glob("lotss_*.csv.gz"))
    expected = [f"lotss_{a:05d}.csv.gz" for a in range(0, len(src), 500)]
    if [p.name for p in parts] != expected:
        raise ValueError("Missing LoTSS chunks")
    neighbours = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    xyz = unit_xyz(src.ra, src.dec)
    neighbours["sep_as"] = separation(
        xyz[neighbours.rowid], unit_xyz(neighbours.ra, neighbours.dec)
    )
    neighbours.sort_values(["rowid", "sep_as"], inplace=True)
    valid = neighbours[neighbours.sep_as <= 25]
    nearest = valid.drop_duplicates("rowid").set_index("rowid").reindex(src.index)
    nearest.rename(
        columns={
            "source_name": "name",
            "total_flux": "flux_jy",
            "e_total_flux": "err_jy",
            "peak_flux": "peak_jy_beam",
            "isl_rms": "rms_jy_beam",
        },
        inplace=True,
    )
    for col in ("flux_jy", "err_jy", "peak_jy_beam", "rms_jy_beam"):
        nearest[col] /= 1000
    for radius in [25, 45, 75, 135, 180]:
        subset = neighbours[neighbours.sep_as <= radius].groupby("rowid")
        nearest[f"n{radius}"] = subset.size().reindex(src.index, fill_value=0)
        nearest[f"sum{radius}_jy"] = subset.total_flux.sum().reindex(src.index) / 1000
    nearest.add_prefix("lotss_").to_csv(out / "match_lotss.csv.gz", index_label="row")
    print("lotss", int(nearest.name.notna().sum()), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("lodess_fits")
    p.add_argument("raw_dir", type=Path)
    p.add_argument("out_dir", type=Path)
    p.add_argument(
        "--surveys",
        nargs="+",
        default=list(SURVEYS) + ["lotss"],
        choices=list(SURVEYS) + ["lotss"],
    )
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    source = Table.read(a.lodess_fits)
    positions = pd.DataFrame(
        {"ra": np.array(source["RA"]), "dec": np.array(source["DEC"])}
    )
    for survey in a.surveys:
        if survey == "lotss":
            lotss_match(a.raw_dir, positions, a.out_dir)
        else:
            survey_match(survey, a.raw_dir, positions, a.out_dir)
