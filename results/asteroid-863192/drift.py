"""Check sensitivity to an imposed relative semimajor-axis drift."""

import argparse
import json
from pathlib import Path
import numpy as np
from orbits import Gravity
from encounter import encounter_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kernel", type=Path)
    parser.add_argument("replay_directory", type=Path)
    args = parser.parse_args()
    summary = json.loads((args.replay_directory / "summary.json").read_text())
    rates = np.array([-3, -1, 0, 1, 3]) * 1e-4
    states = [
        np.tile(summary[k]["barycentric_icrf_au_au_per_day"], (len(rates), 1))
        for k in ("2014 SD315", "2025 UC18")
    ]
    result, _ = encounter_grid(
        Gravity(args.kernel), *states, step_days=2, drift_au_myr=rates
    )
    perturbed, _ = encounter_grid(
        Gravity(args.kernel),
        states[0][:1],
        states[1][:1],
        step_days=1,
        massive_asteroids=True,
    )
    center = perturbed[0, 1]
    local = np.arange((center - 4) * 365.25, (center + 4) * 365.25, 1.0) / 365.25
    fine, _ = encounter_grid(
        Gravity(args.kernel),
        states[0][:1],
        states[1][:1],
        step_days=0.25,
        massive_asteroids=True,
        lookbacks=local,
    )
    np.savetxt(
        args.replay_directory / "massive-asteroids.csv",
        np.vstack((perturbed, fine)),
        delimiter=",",
        header="sampled_min_km,lookback_years,relative_m_s",
        comments="",
        fmt="%.7f",
    )
    window = np.arange(1832 * 365.25, 1841 * 365.25, 1.0) / 365.25
    refinements = []
    for step in (0.5, 0.25, 0.125):
        best, _ = encounter_grid(
            Gravity(args.kernel),
            states[0][:1],
            states[1][:1],
            step_days=step,
            lookbacks=window,
        )
        refinements.append([step, *best[0]])
    np.savetxt(
        args.replay_directory / "refined.csv",
        refinements,
        delimiter=",",
        header="step_days,sampled_min_km,lookback_years,relative_m_s",
        comments="",
        fmt="%.7f",
    )
    np.savetxt(
        args.replay_directory / "drift.csv",
        np.column_stack((rates, result)),
        delimiter=",",
        header="secondary_drift_au_myr,sampled_min_km,lookback_years,relative_m_s",
        comments="",
        fmt="%.7f",
    )


if __name__ == "__main__":
    main()
