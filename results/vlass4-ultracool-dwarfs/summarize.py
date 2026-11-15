"""Calculate the archive search denominator and a conditional bright-source limit."""

import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import beta

# Flux density in mJy at a distance in parsecs, expressed in erg/s/Hz.
LUM_FACTOR = 4 * np.pi * (3.085677581491367e18**2) * 1e-26


def summarize(screen, objects):
    valid = screen.dropna(subset=["snr"]).copy()
    valid["snr"] = valid.snr.astype(float).round(3)
    controls = valid[[f"control{i}" for i in range(8)]].to_numpy()
    plx = objects.reindex(valid.id).set_axis(valid.index)
    sample = valid[
        (plx.plx_formula > 0)
        & (plx.plxerr_formula > 0)
        & (plx.plx_formula / plx.plxerr_formula >= 10)
        & valid.sp.ge(20)
    ].copy()
    sample["d_pc"] = 1000 / plx.loc[sample.index, "plx_formula"]
    # Treat co-located catalog components as one sky position. Use only one
    # row per group within one maximum nominal beam size (4 arcsec).
    kept = []
    for r in sample.sort_values("id").itertuples():
        if any(
            np.hypot(
                ((r.ra_obs - a.ra_obs + 180) % 360 - 180)
                * np.cos(np.deg2rad(r.dec_obs)),
                r.dec_obs - a.dec_obs,
            )
            * 3600
            < 4
            for a in kept
        ):
            continue
        kept.append(r)
    sample = sample.loc[[r.Index for r in kept]]
    limits = []
    for threshold in (1e14, 1e15):
        expected = threshold / (LUM_FACTOR * sample.d_pc**2)
        eligible = sample[expected >= 10 * sample.rms_mjy]
        n = len(eligible)
        k = int((eligible.snr >= 6).sum())
        limits.append(
            dict(
                luminosity_erg_s_hz=threshold,
                eligible_parallax_TY_positions=n,
                threshold_crossings=k,
                upper_95_binomial_if_independent=round(
                    float(beta.ppf(0.95, k + 1, n - k)), 5
                )
                if n > k
                else None,
            )
        )
    return dict(
        footprint_positions=len(screen),
        measured=len(valid),
        median_mad_rms_mjy=float(valid.rms_mjy.median()),
        positive_at_6=int((valid.snr >= 6).sum()),
        controls=int(np.isfinite(controls).sum()),
        controls_at_6=int((controls >= 6).sum()),
        fourth_epoch_weak_lead=dict(
            threshold_in_pixel_noise_units=3.681,
            on_target_at_or_above=int((valid.snr >= 3.681).sum()),
            offset_positions_at_or_above=int((controls >= 3.681).sum()),
            T8_or_later_targets=int(valid.sp.ge(28).sum()),
            T8_or_later_targets_at_or_above=int(
                (valid.sp.ge(28) & valid.snr.ge(3.681)).sum()
            ),
        ),
        conditional_limits=limits,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    args = p.parse_args()
    c = args.cache
    print(
        json.dumps(
            summarize(
                pd.read_csv(c / "screen.csv"),
                pd.read_csv(c / "ucd.csv", low_memory=False),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
