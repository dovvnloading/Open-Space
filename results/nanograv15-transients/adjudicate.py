"""Join the measured post-discovery controls to every original surviving fit."""

import csv
import json
from collections import Counter
from pathlib import Path


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def number(value, digits=4):
    return "—" if value is None or value == "" else f"{float(value):.{digits}g}"


def probability(row, key, threshold):
    return row.get(key) not in (None, "") and float(row[key]) <= threshold


def audit(directory, singles, array):
    evidence = Path(directory)
    coherent = read_rows(evidence / "coherence/events.csv")
    dm = read_rows(evidence / "dm-control/fixed_controls.csv")
    survivors = array["survivors"]
    if len(coherent) != len(survivors) or len(dm) != 6 * len(survivors):
        raise ValueError("post-discovery ledger is incomplete")
    controls = {(int(r["survivor_ordinal"]), r["control"]): r for r in dm}
    review = []
    lines = ["# Review of surviving timing fits", "",
             "Each array row is a fit at a selected epoch. Several shapes and sky directions can describe one excursion. "
             "The tables give every fit that passed the original internal checks. The complete participant lists, "
             "amplitude estimates and uncertainties are linked below.", "",
             "## Individual pulsars", ""]
    individual = [(p["pulsar"], r) for p in singles for r in p["survivors"]]
    if individual:
        lines.extend(["| Pulsar | Epoch (MJD TDB) | Template | Width (days) | Coefficient ± SE (µs) | Global p bound | EPTA directional p | Spectrum |",
                      "| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |"])
        for psr, event in individual:
            ext = event["checks"]["epta"]
            lines.append(f"| {psr} | {number(event['epoch_mjd'], 9)} | {event['template']} | "
                         f"{number(event['width_d'])} | {number(event['amplitude_us'])} ± {number(event['se_us'])} | "
                         f"{number(event['p_global'])} | {number(ext.get('directional_p'))} | {event['spectrum']} |")
        lines.extend(["", "Every listed individual fit passed the fixed wideband and independent EPTA checks and the "
                      "noise and visit-deletion checks in its [individual record](single)."])
    else:
        lines.append("No individual fit passed all of the preregistered checks. The [per-pulsar dispositions](individual-dispositions.csv) retain primary crossings and reasons for withholding them.")
    lines.extend(["", "## Common timing patterns", "",
                  "The empirical p uses the maximum over the entire bank after randomizing pulsar signs. "
                  "EPTA p values use Holm adjustment across all internally surviving common fits. "
                  "The DM column gives the coefficient of the achromatic signal on the discovery axis after fitting "
                  "an independent dispersion delay on each usable telescope-day. For a negative clock signal, "
                  "positive alignment in this control means agreement with that negative discovery sign. "
                  "A dash means the fit was unavailable.", "",
                  "| ID | Epoch (MJD TDB) | Pattern and shape | Discovery coefficient | Units | Pulsars | Sign-test p | WB | EPTA adjusted p | DM control ± SE |",
                  "| ---: | ---: | --- | ---: | --- | ---: | ---: | --- | ---: | ---: |"])
    tallies = Counter()
    for ordinal, (event, r) in enumerate(zip(survivors, coherent), 1):
        if (int(r["survivor_ordinal"]) != ordinal or r["family"] != event["family"]
                or r["shape"] != event["shape"] or float(r["epoch_mjd"]) != event["epoch_mjd"]):
            raise ValueError("post-discovery rows differ from the frozen array fits")
        free = controls[ordinal, "dm_free_all"]
        raw = controls[ordinal, "matched_raw_all"]
        is_sign = probability(r, "mc_high", .005)
        subset = r["pass_both_sign_and_subset_checks"].lower() == "true"
        external = probability(r, "epta_p_holm", .05)
        wideband = r["wideband_check"] == "same-direction support"
        chrom = probability(free, "p_fixed_one_sided", .05)
        all_support = subset and external and wideband and chrom
        tallies.update({"full_sign": int(is_sign), "noise_subset_and_sign": int(subset),
                        "external_testable": int(bool(r.get("epta_p_holm"))), "external_holm": int(external),
                        "dm_testable": int(bool(free.get("p_fixed_one_sided"))), "dm_support": int(chrom),
                        "dm_expected_snr_ge3": int(float(free.get("expected_snr_at_discovery_amplitude") or 0) >= 3),
                        "matched_raw_support": int(probability(raw, "p_fixed_one_sided", .05)),
                        "wb_support": int(wideband), "all_support": int(all_support)})
        wb = "supports" if wideband else r["wideband_check"]
        lines.append(f"| {ordinal} | {number(event['epoch_mjd'], 9)} | {r['family']} {r['shape']} | "
                     f"{number(event['amplitude'])} | {r['units']} | {r['full_n']} | {number(r['p_sign_max_bank'])} | "
                     f"{wb} | {number(r.get('epta_p_holm'))} | {number(free.get('signed_amplitude'))} ± {number(free.get('se'))} |")
        review.append({"survivor_ordinal": ordinal, "family": r["family"], "shape": r["shape"],
                       "epoch_mjd": event["epoch_mjd"], "sky_id": event["sky_id"], "pulsars": event["pulsars"],
                       "ra_deg": event["sky_ra_deg"], "dec_deg": event["sky_dec_deg"], "units": r["units"],
                       "amplitude": event["amplitude"], "components": event["components"],
                       "component_se": event["component_se"], "p_joint_primary": event["p_joint"],
                       "parametric_leave_one": event["internal_check"], "coherence": r, "dm_free": free,
                       "matched_visits_before_dm_removal": raw,
                       "good_noise_subset_and_sign": subset, "external_support_after_holm": external,
                       "dm_support": chrom, "all_followups_support": all_support})
    if not survivors:
        lines.append("| — | — | No internally surviving fit | — | — | — | — | — | — | — |")
    lines.extend(["", "Detailed results include the signed NB, WB and EPTA amplitudes, their uncertainties, per-pulsar "
                  "contributions, sky coordinates, full participant lists, diagnostic-subset checks, and separate "
                  "Arecibo and Green Bank fits. A local p value for a follow-up fit is descriptive. Sparse frequency or sky "
                  "coverage can make a split inconclusive. [Full row-by-row record](event-review.json) · "
                  "[Per-pulsar contributions](coherence/pulsars.csv) · [Telescope checks](coherence/telescope_checks.csv) · "
                  "[Matched dispersion controls](dm-control/fixed_controls.csv) · [Interpretation protocol](../ADJUDICATION.md)", ""])
    (evidence / "EVENTS.md").write_text("\n".join(lines))
    (evidence / "event-review.json").write_text(json.dumps(review, indent=2, allow_nan=False) + "\n")
    summary = {key: tallies[key] for key in ("full_sign", "noise_subset_and_sign", "external_testable",
               "external_holm", "dm_testable", "dm_support", "dm_expected_snr_ge3", "matched_raw_support",
               "wb_support", "all_support")}
    summary.update(array_internal=len(survivors), individual=len(individual))
    (evidence / "review-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
