"""External VLSSr parent and forced measurements in the central 4 degrees of good fields.

The parent uses >1 Jy peak and convolved axes <90 arcsec at 73.8 MHz. A
field-level nonmatch is no released LoDeSS source centroid within 60 arcsec.
The local rms is measured empirically for each signed forced PSF amplitude.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree

from match_catalogues import read_tdat, unit_xyz
from photometry import RadioImage


def run(data, comp):
    vl = read_tdat(
        comp / "raw/vlssr.tdat.gz",
        [
            "name",
            "ra",
            "dec",
            "flux_74_mhz",
            "flux_74_mhz_error",
            "major_axis",
            "minor_axis",
        ],
        0,
    )
    vl = vl[
        (vl.flux_74_mhz > 1000) & (vl.major_axis < 90) & (vl.minor_axis < 90)
    ].copy()
    fields = pd.read_csv(data / "fields.csv").query("quality=='good'")
    catalog = Table.read(data / "sources.fits")
    p = unit_xyz(catalog["RA"], catalog["DEC"])
    memberships = pd.read_csv(data / "source_fields.csv").set_index("row")
    field_rows = []
    measurements = []
    for f in fields.itertuples(index=False):
        center = unit_xyz([f.ra_deg], [f.dec_deg])[0]
        pos = unit_xyz(vl.ra, vl.dec)
        distances = np.rad2deg(
            2 * np.arcsin(np.minimum(1, np.linalg.norm(pos - center, axis=1) / 2))
        )
        area = distances < 4
        entries = vl[area].copy()
        index = memberships.index[memberships.field == f.field]
        tree = cKDTree(p[index])
        d, j = tree.query(pos[area])
        entries["closest_row"] = index[j]
        entries["closest_as"] = np.rad2deg(2 * np.arcsin(d / 2)) * 3600
        entries["field"] = f.field
        entries["r_deg"] = distances[area]
        field_rows.append(entries)
        with RadioImage(data / "images" / f"{f.field}-I.fits") as image:
            for e in entries.query("closest_as>60").itertuples():
                m = image.measure(e.ra, e.dec)
                if m:
                    measurements.append(
                        dict(
                            name=e.name,
                            ra=e.ra,
                            dec=e.dec,
                            field=f.field,
                            r_deg=e.r_deg,
                            s74_jy=e.flux_74_mhz / 1000,
                            rms74_jy=e.flux_74_mhz_error / 1000,
                            closest_as=e.closest_as,
                            **m,
                        )
                    )
    out = comp / "results"
    parent = pd.concat(field_rows, ignore_index=True)
    parent.to_csv(out / "vlssr_parent.csv.gz", index=False)
    phot = pd.DataFrame(measurements)
    phot["limit95_jy"] = np.maximum(phot.psf_jy, 0) + 1.645 * phot.psf_noise_jy
    phot["s74_to_limit95"] = phot.s74_jy / phot.limit95_jy
    phot.sort_values("s74_to_limit95", ascending=False, inplace=True)
    phot.to_csv(out / "vlssr_nondetections.csv", index=False)
    unique = phot.drop_duplicates("name").reset_index(drop=True)
    unique[["name", "ra", "dec"]].to_csv(
        out / "reverse_positions.csv", index_label="row"
    )
    Table([unique.ra, unique.dec], names=["RA", "DEC"]).write(
        out / "reverse_positions.fits", overwrite=True
    )
    print(
        "parent",
        len(parent),
        parent.name.nunique(),
        "nonmatches",
        len(phot),
        phot.name.nunique(),
    )
    print(
        phot[
            [
                "name",
                "field",
                "r_deg",
                "s74_jy",
                "psf_jy",
                "psf_noise_jy",
                "s74_to_limit95",
            ]
        ]
        .head(24)
        .to_string(index=False)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("survey_data", type=Path)
    p.add_argument("comparison", type=Path)
    args = p.parse_args()
    run(args.survey_data, args.comparison)
