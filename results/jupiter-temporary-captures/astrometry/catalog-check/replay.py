"""Evaluate catalog centroids at orbital predictions for 47 DECam exposures."""

import argparse
import json
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from skyfit import Gravity, SkyFit


def planned_sql(rows):
    conditions = []
    for row in rows.itertuples():
        base, tag = row.filename.removesuffix(".fits.fz").rsplit("_", 1)
        names = sorted({tag, "ls9", "v1", "v2", "vx", "ls10"})
        aliases = ", ".join(f"'{base}_{value}'" for value in names)
        conditions.append(
            f"(exposure IN ({aliases}) AND ra BETWEEN {row.pred_ra - 0.003} AND "
            f"{row.pred_ra + 0.003} AND dec BETWEEN {row.pred_dec - 0.003} "
            f"AND {row.pred_dec + 0.003})"
        )
    return (
        "SELECT measid,objectid,exposure,ra,dec,mjd,mag_auto,magerr_auto,flags "
        "FROM nsc_dr2.meas WHERE " + " OR ".join(conditions)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kernel", type=Path, help="Recompute predictions from the fitted orbit"
    )
    parser.add_argument(
        "--live", action="store_true", help="Repeat the NOIRLab catalog queries"
    )
    options = parser.parse_args()
    rows = pd.read_csv(HERE / "exposures.csv")
    if options.kernel:
        model = SkyFit(rows.assign(site="W84"), Gravity(options.kernel))
        orbit = json.loads((HERE.parent / "measurements/fit-summary.json").read_text())[
            "all"
        ]
        error, _ = model.prediction(np.array(orbit["barycentric_state_au_day"]))
        east, north = error.reshape(-1, 2).T
        np.testing.assert_allclose(
            rows.pred_ra,
            rows.ra + east / (3600 * np.cos(np.deg2rad(rows.dec))),
            atol=2e-8,
            rtol=0,
        )
        np.testing.assert_allclose(
            rows.pred_dec, rows.dec + north / 3600, atol=2e-8, rtol=0
        )
    detections = pd.read_csv(HERE / "catalog.csv")
    if options.live:
        import io

        pieces = []
        for i in range(0, len(rows), 3):
            sql = planned_sql(rows.iloc[i : i + 3])
            url = "https://datalab.noirlab.edu/query/query?" + urlencode(
                {"sql": sql, "ofmt": "csv", "async": "False"}
            )
            req = Request(url, headers={"X-DL-AuthToken": "anonymous.0.0.anon_access"})
            with urlopen(req, timeout=150) as response:
                pieces.append(pd.read_csv(io.BytesIO(response.read())))
        detections = pd.concat(pieces, ignore_index=True)
    results = []
    for row in rows.itertuples():
        hits = detections.loc[detections.exposure.str[:17] == row.filename[:17]]
        separations = (
            np.hypot(
                (hits.ra - row.pred_ra) * np.cos(np.deg2rad(row.pred_dec)),
                hits.dec - row.pred_dec,
            )
            * 3600
        )
        results.append(
            {
                "utc": row.utc,
                "catalog_rows": len(hits),
                "nearest_arcsec": min(separations, default=None),
            }
        )
    print(pd.DataFrame(results).to_csv(index=False, float_format="%.4f"), end="")


if __name__ == "__main__":
    main()
