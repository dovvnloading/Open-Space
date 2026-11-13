"""Plan and measure independent later-epoch TESS observations of selected variables."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from .tess import assess, blind, enough, flatten, load_product
from .tess_archive import DEFAULT_PLAN, atomic_bytes, cached_product, inventory


def sinusoid(curve):
    """Return the strongest independently fit 0.2–8 day sinusoid in this light curve."""
    import numpy as np
    from astropy.timeseries import LombScargle

    time, flux, error = curve
    span = time[-1] - time[0]
    low, high = 1 / 8, 1 / 0.2
    trials = int(np.ceil((high - low) * span * 15)) + 1
    grid = np.linspace(low, high, trials)
    model = LombScargle(time, flux, error, fit_mean=True, center_data=True)
    coarse = model.power(grid, method="fast", assume_regular_frequency=True)
    peak = int(np.argmax(coarse))
    fine = np.linspace(grid[max(0, peak - 2)], grid[min(trials - 1, peak + 2)], 101)
    powers = model.power(fine, method="chi2")
    index = int(np.argmax(powers))
    return {"period_days": float(1 / fine[index]), "standard_power": float(powers[index]),
            "coarse_trials": trials}


def make_later_plan(tics, sectors, inventory_dir, cache):
    """Choose one later observation per TIC using only product-list membership."""
    if not tics or any(tic <= 0 for tic in tics) or len(tics) != len(set(tics)):
        raise ValueError("Supply unique, positive TIC IDs")
    if not sectors or sectors != sorted(set(sectors)):
        raise ValueError("Supply sectors in strictly increasing order")
    listings = [inventory(sector, inventory_dir) for sector in sectors]
    targets = []
    for tic in tics:
        for names, manifest in listings:
            if tic in names:
                name = names[tic]
                product = {"filename": name, "sector": manifest["sector"],
                           "uri": "mast:TESS/product/" + name}
                _, pinned = cached_product(product, cache)
                targets.append({"tic_id": tic, "product": pinned})
                break
        else:
            targets.append({"tic_id": tic, "product": None})
    return {"schema_version": 1,
            "selection": "Earliest available listed sector for each input TIC, using only inventory membership.",
            "inventories": [metadata for _, metadata in listings], "targets": targets}


def run(later_plan, discovery_plans, cache):
    from astropy.timeseries import LombScargle

    early = {target["tic_id"]: target["products"][0]
             for plan in discovery_plans for target in plan["targets"]}
    cases = []
    for record in later_plan["targets"]:
        tic = record["tic_id"]
        if record["product"] is None:
            cases.append({"tic_id": tic, "outcome": "no_later_product", "sources": []})
            continue
        if record["product"]["sector"] <= early[tic]["sector"]:
            raise ValueError(f"TIC {tic} needs an independent later sector")
        first_data, first_source = load_product(early[tic], cache, tic)
        last_data, last_source = load_product(record["product"], cache, tic)
        sources = [first_source, last_source]
        first, last = flatten(first_data), flatten(last_data)
        if not enough(first) or not enough(last):
            cases.append({"tic_id": tic, "outcome": "insufficient_coverage", "sources": sources})
            continue
        fit, revisit = blind(first), blind(last)
        expected, measured = fit["period_days"], revisit["period_days"]
        first_sine, later_sine = sinusoid(first), sinusoid(last)
        first_frequency = 1 / first_sine["period_days"]
        last_model = LombScargle(*last)
        period_difference = later_sine["period_days"] / first_sine["period_days"] - 1
        check, independent = assess(last, fit), assess(last, revisit)
        check.pop("events")
        independent.pop("events")
        cases.append({
            "tic_id": tic, "outcome": "measured", "sources": sources,
            "baseline_start_to_start_days": float(last[0][0] - first[0][0]),
            "discovery_bls": fit, "later_blind_bls": revisit,
            "later_fixed_ephemeris_depth": check, "later_blind_depth": independent,
            "blind_period_fractional_difference": measured / expected - 1,
            "discovery_blind_sinusoid": first_sine, "later_blind_sinusoid": later_sine,
            "later_sinusoid_power_at_first_frequency": float(last_model.power(first_frequency)),
            "later_sinusoid_power_at_twice_first_frequency": float(last_model.power(2 * first_frequency)),
            "sinusoid_period_fractional_difference": period_difference,
            "sinusoid_period_agrees_within_2pct": abs(period_difference) <= 0.02,
        })
    return {"schema_version": 1, "experiment": "later_sector_variability_recheck_v1",
            "sample_role": "selected_followup", "cases": cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan", help="Select and checksum a later epoch for each requested TIC")
    plan.add_argument("--tics", type=int, nargs="+", required=True)
    plan.add_argument("--sectors", type=int, nargs="+", required=True)
    plan.add_argument("--inventories", type=Path, required=True)
    measure = commands.add_parser("run", help="Search a later observation independently and at the original phase")
    data = DEFAULT_PLAN.parent
    measure.add_argument("--later-plan", type=Path, default=data / "tess-third-plan.json")
    measure.add_argument("--discovery-plans", type=Path, nargs="+",
                         default=[DEFAULT_PLAN, data / "tess-expansion-plan.json"])
    for command in (plan, measure):
        command.add_argument("--cache", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "plan":
        output = make_later_plan(args.tics, args.sectors, args.inventories, args.cache)
    else:
        snapshots = [path.read_bytes() for path in (args.later_plan, *args.discovery_plans)]
        output = run(json.loads(snapshots[0]), [json.loads(v) for v in snapshots[1:]], args.cache)
        output["input_plan_sha256"] = [sha256(v).hexdigest() for v in snapshots]
    atomic_bytes(args.output, (json.dumps(output, indent=2, allow_nan=False) + "\n").encode())


if __name__ == "__main__":
    main()
