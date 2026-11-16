"""Export collocated OCI observations for independent checks of ranked date pairs."""

import numpy as np
import pandas as pd
from analyze import BANDS, REGIONS, ROOT, quality, read, window


def select_pairs():
    a, lat, lon, depth = read()
    dates = (
        pd.to_datetime(a["time"], unit="s", utc=True).strftime("%Y-%m-%d").to_numpy()
    )
    valid = quality(a)
    frames = []
    ranking = pd.read_csv(ROOT / "regional_changes.csv").query("coverage_pass")
    targeted = ranking.sort_values("delta_C", key=abs, ascending=False)
    for row in targeted.itertuples():
        i, j = [dates.tolist().index(d) for d in (row.start, row.end)]
        mask = window(lat, lon, REGIONS[row.region]) & (depth >= 200)
        mask &= valid[i] & valid[j]
        iy, ix = np.nonzero(mask)
        f = pd.DataFrame(
            dict(
                region=row.region,
                start=row.start,
                end=row.end,
                lat=lat[iy],
                lon=lon[ix],
                depth_m=depth[mask],
            )
        )
        for name, t in (("start", i), ("end", j)):
            for b in BANDS:
                f[f"{name}_Rrs_{b}"] = a[f"Rrs_{b}"][t, iy, ix]
        frames.append(f)
    pd.concat(frames).to_csv(ROOT / "footprints.csv", index=False, float_format="%.9g")


if __name__ == "__main__":
    select_pairs()
