"""Forced two-band fluxes at the Gaia background star and detector controls."""

import argparse
import csv
from pathlib import Path

import numpy as np
from miri_reduce import calibrators, load, matched, median, write_table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("."))
    a = ap.parse_args()
    target_offset = np.array([4.019374664255579, 3.661244990280572])
    out = []
    for band, sci, ref, bg, ref_bg, cal in [
        ("F1065C", "38", "35", "40", "33", "25"),
        ("F1550C", "37", "36", "39", "34", "27"),
    ]:
        science, h, s0 = load(next(a.inputs.glob("jw090560" + sci + "*fits")))
        cd = np.array([[h["CD1_1"], h["CD1_2"]], [h["CD2_1"], h["CD2_2"]]]) * 3600
        center = np.array([h["CRPIX1"] - 1, h["CRPIX2"] - 1])
        x, y = center + np.linalg.solve(cd, target_offset)
        ix, iy = round(x), round(y)
        data = np.load(a.out / (band.lower() + "_miri_maps.npz"))
        amp = data["amplitude"]
        trans = data["relative_transmission"]
        snr = data["snr"]
        cut = data["threshold"]
        # Report the nearest blind-grid point; no positional refit or peak selection.
        f = amp[iy, ix] * h["PIXAR_SR"] * 1e12 / trans[iy, ix]
        noise_rows = list(
            csv.DictReader((a.out / (band.lower() + "_miri_noise.csv")).open())
        )
        r = np.hypot(x - center[0], y - center[1]) * np.sqrt(abs(np.linalg.det(cd)))
        rr = next(
            q for q in noise_rows if float(q["lo_arcsec"]) <= r < float(q["hi_arcsec"])
        )
        err = float(rr["sigma"]) * h["PIXAR_SR"] * 1e12 / trans[iy, ix]
        psf, _, _, _ = calibrators(sorted(a.inputs.glob("jw015380" + cal + "*fits")))
        kernel = psf[11:20, 11:20]
        for kind, obs in [
            ("science", sci),
            ("reference", ref),
            ("target background", bg),
            ("reference background", ref_bg),
        ]:
            for path in sorted(a.inputs.glob("jw090560" + obs + "*fits")):
                z, hw, _ = load(path)
                raw = matched(median(z, axis=0), kernel)[iy, ix] * hw["PIXAR_SR"] * 1e12
                out.append(
                    dict(
                        band=band,
                        data_kind=kind,
                        filename=path.name,
                        detector_x0=ix,
                        detector_y0=iy,
                        raw_template_uJy=raw,
                        final_forced_uJy=float(f),
                        ring_noise_uJy=err,
                        final_snr=float(snr[iy, ix]),
                        in_search=int(np.isfinite(cut[iy, ix])),
                        final_relative_transmission=float(trans[iy, ix]),
                    )
                )
    write_table(a.out / "miri_source_checks.csv", list(out[0]), out)
    for row in out:
        if row["data_kind"] == "science":
            print(row)


if __name__ == "__main__":
    main()
