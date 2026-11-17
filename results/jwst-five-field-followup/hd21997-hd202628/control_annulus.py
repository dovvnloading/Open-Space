"""Count positive/negative matched-filter maxima outside the science annulus."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import maximum_filter


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    sys.path.insert(0, str(a.merged / "nircam"))
    import nircam_reduce as rd

    # Original annular threshold estimator, moved outside the discovery zone.
    rd.INNER, rd.OUTER = 6.1, 8.5
    res = []
    for field in ["HD21997", "HD202628"]:
        meta = json.loads((a.out / (field + "_measurement.json")).read_text())
        z = np.load(a.out / (field + "_rdi_maps.npz"))
        amp = z["amplitude"].astype(float)
        scale = meta["scale_arcsec"]
        cx, cy = meta["center_xy0"]
        y, x = np.indices(amp.shape)
        radius = np.hypot(x - cx, y - cy) * scale
        cut, snr, rings = rd.thresholds(amp, radius, [cx, cy], scale)
        valid = np.isfinite(cut) & np.isfinite(amp)
        c_area = valid.sum() * scale**2
        s_area = ((radius >= 2.5) & (radius < 5.5) & np.isfinite(amp)).sum() * scale**2
        pos = (amp > cut) & (
            amp == maximum_filter(np.nan_to_num(amp, nan=-1e9), size=5)
        )
        neg = np.zeros_like(pos)
        for ring in rings:
            neg |= (
                (radius >= ring["lo_arcsec"])
                & (radius < ring["hi_arcsec"])
                & (amp < 2 * ring["mean"] - ring["cutoff"])
            )
        neg &= (-amp) == maximum_filter(np.nan_to_num(-amp, nan=-1e9), size=5)
        res.append(
            dict(
                field=field,
                control_inner_arcsec=rd.INNER,
                control_outer_arcsec=rd.OUTER,
                control_area_arcsec2=c_area,
                search_area_arcsec2=s_area,
                positive=int(pos.sum()),
                negative=int(neg.sum()),
                scaled_positive=pos.sum() * s_area / c_area,
            )
        )
    rd.write_table(a.out / "F444W_background.csv", list(res[0]), res)


if __name__ == "__main__":
    main()
