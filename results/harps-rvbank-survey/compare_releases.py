"""Compare frozen SERVAL periods reported under the identical RVBank series name."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from runtime import require_survey_workflow


def read(path):
    return pd.read_csv(path, keep_default_na=False, dtype={"star": str}, low_memory=False)


def run(full, corrected, supplement):
    require_survey_workflow()
    inherited = read(full/"vetting.csv")
    updated = read(supplement/"all-activity.csv")
    inherited = inherited.merge(updated[["star","signal_index","new_rv_status","delta_bic"]],
                                on=["star","signal_index"],validate="one_to_one")
    series = inherited.groupby("star",sort=False)
    partial = read(corrected/"vetting.csv")
    reference_span = read(full/"stars.csv").set_index("star").baseline_d
    corrected_span = read(corrected/"stars.csv").set_index("star").baseline_d
    records = []
    for peak in partial.itertuples(index=False):
        row = dict(star=peak.star, signal_index=peak.signal_index, period_d=peak.period_d,
                   k_ms=peak.k_ms, fap=peak.fap, rv_status=peak.rv_status,
                   external_status=peak.external_status, prior_signal=peak.prior_signal)
        records.append(row)
        if peak.star not in series.groups:
            row["comparison"] = "no same-name fitted series in full release"
            continue
        family = series.get_group(peak.star)
        freq = 1. / float(peak.period_d)
        diff = (1. / family.period_d - freq).abs()
        span = min(float(corrected_span[peak.star]), float(reference_span[peak.star]))
        same = diff.le(min(.05*freq, 1./span))
        selected = family[same] if any(same) else family
        counterpart = selected.loc[diff.loc[selected.index].idxmin()]
        row.update(full_signal_index=int(counterpart.signal_index), full_period_d=counterpart.period_d,
                   full_k_ms=counterpart.k_ms, full_fap=counterpart.fap,
                   full_rv_status=counterpart.new_rv_status, full_external_status=counterpart.external_status,
                   full_prior_signal=counterpart.prior_signal,
                   separation_rayleigh=float(diff.loc[counterpart.name])*span)
        if not any(same):
            row["comparison"] = "distinct from fitted full-release periods"
        elif counterpart.fap<=.01 and counterpart.delta_bic>=10:
            row["comparison"] = "same period above thresholds in full release"
        else:
            row["comparison"] = "same period below thresholds in full release"
    table = pd.DataFrame(records)
    table.to_csv(corrected/"reduction-comparison.csv", index=False)
    table[table.rv_status.eq("pass")].to_csv(corrected/"survivor-comparison.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["full","corrected","supplement"]:
        p.add_argument("--"+name, required=True, type=Path)
    a = p.parse_args()
    run(a.full, a.corrected, a.supplement)
