"""Re-measure externally supported reverse-search nondetections in overlapping fields."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from match_catalogues import unit_xyz
from photometry import RadioImage


def run(survey_data, comparison):
    out = comparison / "results"
    cases = pd.read_csv(out / "reverse_counterparts.csv.gz")
    selected = cases[
        (cases.psf_jy < 3 * cases.psf_noise_jy)
        & (
            (cases.nvss_flux_jy > 0.02)
            | (cases.gmrtas150m_flux_jy > 0.08)
            | (cases.lotss_flux_jy > 0.08)
        )
    ]
    pos = unit_xyz(selected.ra, selected.dec)
    fields = pd.read_csv(survey_data / "fields.csv")
    measures = []
    for f in fields.itertuples(index=False):
        path = survey_data / "images" / f"{f.field}-I.fits"
        distances = np.rad2deg(
            2
            * np.arcsin(
                np.minimum(
                    1,
                    np.linalg.norm(pos - unit_xyz([f.ra_deg], [f.dec_deg]), axis=1) / 2,
                )
            )
        )
        use = distances < 10
        if not path.exists() or not np.any(use):
            continue
        with RadioImage(path) as image:
            for row, sep in zip(selected[use].itertuples(index=False), distances[use]):
                m = image.measure(row.ra, row.dec)
                if m:
                    measures.append(
                        dict(
                            name=row.name,
                            field=f.field,
                            quality=f.quality,
                            r_deg=sep,
                            date_obs=image.hdus[0].header["DATE-OBS"],
                            s74_jy=row.s74_jy,
                            **m,
                        )
                    )
    t = pd.DataFrame(measures)
    t["snr"] = t.psf_jy / t.psf_noise_jy
    t.to_csv(out / "vlssr_repeat.csv", index=False)
    print(t.pivot(index="name", columns="field", values="psf_jy").round(3))
    print(t.groupby("field").snr.agg(["min", "median", "max"]))
    print(
        "Independent >=5 sigma detections:",
        t.query("field!='P136+36' and snr>=5").name.nunique(),
        "out of",
        len(selected),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("survey_data", type=Path)
    p.add_argument("comparison", type=Path)
    args = p.parse_args()
    run(args.survey_data, args.comparison)
