"""Compare native 2024 HFR spectra with seven-second antenna-bias telemetry."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import cdflib
import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
AEB_URL = ("https://research.ssl.berkeley.edu/data/psp/data/sci/fields/l2/aeb/"
           "2024/11/psp_fld_l2_aeb_20241106_v01.cdf")
AEB_SHA = "60ad4d18afc61471b6585ea8400f817bbb8d515b770ca2fa5cdbd40ace1f741e"
CENTER = np.datetime64("2024-11-06T18:44:00")
BANDS = [(1.3, 2), (2, 4), (4, 8), (8, 16)]


def get_cdf(path, url, digest):
    data = path.read_bytes() if path.exists() else urlopen(url, timeout=120).read()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Checksum mismatch for {path.name}")
    if not path.exists():
        path.write_bytes(data)
    return cdflib.CDF(path)


def seconds(cdf, epoch):
    return (cdflib.cdfepoch.to_datetime(cdf.varget(epoch)) - CENTER) / np.timedelta64(1, "s")


def native_power(cdf, antenna):
    key = "hfr_auto_averages_" + antenna
    var = "psp_fld_l2_rfs_" + key
    t = seconds(cdf, "epoch_" + key)
    mask = ((abs(t) < 600) & (cdf.varget(var+"_gain") == 1) &
            (cdf.varget(var+"_hl") == 0) & (cdf.varget(var+"_nsum") == 40))
    f = cdf.varget("frequency_"+key)[mask] / 1e6
    if not np.all(f == f[0]):
        raise ValueError("Frequency grid changed during the comparison")
    flux = cdf.varget(var)[mask].astype(float)
    flux[(flux <= 0) | (flux > cdf.varattsget(var)["VALIDMAX"])] = np.nan
    bands = []
    for lo, hi in BANDS:
        selected = flux[:, (f[0] >= lo) & (f[0] < hi)]
        db = np.nanmedian(10*np.log10(selected), axis=1)
        db[np.isfinite(selected).sum(axis=1) < selected.shape[1]/2] = np.nan
        bands.append(db)
    return t[mask], np.column_stack(bands)


def rank_effect(x, powers):
    return np.array([spearmanr(x, band).statistic for band in powers.T])


def r4(x):
    return np.round(x, 4).tolist()


def step_differences(bt, current, times, power, shift):
    # AEB records report settings when packets arrive. Use pairs of adjoining
    # seven-second bins, each delimited by two actual housekeeping timestamps.
    good = (abs(np.diff(bt) - 7) < .15) & (current[:-1] > 0)
    good &= (bt[:-1] > -600) & (bt[1:] < 600)
    means = np.full((len(good), 4), np.nan)
    for i in np.where(good)[0]:
        mask = (times >= bt[i]+shift) & (times < bt[i+1]+shift)
        if mask.any():
            means[i] = np.mean(power[mask], axis=0)
    ix = np.where(np.all(np.isfinite(means[:-1]) & np.isfinite(means[1:]), axis=1))[0]
    x = np.diff(current)[ix]
    y = np.diff(means, axis=0)[ix]
    return x, y


def compare_steps(bt, current, times, power, shift):
    x, y = step_differences(bt, current, times, power, shift)
    rho = rank_effect(x, y)
    xc = x-x.mean()
    slope = xc @ y / (xc @ xc)
    # Descriptive phase-shift diagnostic. The retained pairs can have gaps and
    # the commanded staircase recurs, so these ranks are not calibrated p-values.
    alternatives = np.array([rank_effect(np.roll(x, k), y) for k in range(len(x))])
    return {
        "adjacent_pairs": len(x), "spearman_rho": r4(rho), "db_per_uA_OLS": r4(slope),
        "phase_shift_ranks_abs_ge": (abs(alternatives) >= abs(rho)-1e-10).sum(axis=0).tolist(),
        "phase_shifts_total": len(x),
        "halves_spearman_rho": [r4(rank_effect(x[g], y[g])) for g in np.array_split(np.arange(len(x)), 2)],
    }


def command_overlap(times, current):
    # Match exact DAC settings. Core and short outside windows replicate the radio
    # contrast above. A second split tests inner versus outer shadow step plateaus.
    offset = abs(times)
    center = current[offset <= 180]
    pre = current[(times >= -600) & (times < -240)]
    after = current[(times > 240) & (times <= 600)]
    dt = np.diff(times)
    valid = (abs(dt-7) < .15) & (current[:-1] > 0)
    middle = abs(times[:-1]+dt/2)
    core = current[:-1][valid & (middle <= 60)]
    shoulders = current[:-1][valid & (middle >= 120) & (middle <= 240)]
    return {
        "core_3min_packets": len(center),
        "left_and_right_4_to_10min_packets": [len(pre), len(after)],
        "same_command_in_core_and_both_flanks": len(np.intersect1d(np.intersect1d(center, pre), after)),
        "inner_1min_positive_plateaus": len(core),
        "outer_2_to_4min_positive_plateaus": len(shoulders),
        "shared_DAC_values_inner1_outer2to4": len(np.intersect1d(core, shoulders)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aeb", type=Path, required=True, help="local AEB CDF cache path")
    parser.add_argument("--hfr", type=Path, required=True, help="local HFR CDF cache path")
    args = parser.parse_args()
    event = next(e for e in json.loads((HERE/"manifest.json").read_text()) if e["date"] == "2024-11-06")
    aeb = get_cdf(args.aeb, AEB_URL, AEB_SHA)
    hfr = get_cdf(args.hfr, event["cdf"], event["sha256"])
    out = {"aeb_url": AEB_URL, "aeb_sha256": AEB_SHA, "hfr_sha256": event["sha256"],
           "bands_MHz": BANDS, "nominal_packet_offset_s": 0, "antennas": {}}
    for i, antenna in enumerate(("ch0_V1V2", "ch1_V3V4"), start=1):
        prefix = f"psp_fld_l2_aeb{i}_"
        t = seconds(aeb, f"epoch_aeb{i}")
        current = aeb.varget(prefix+f"V{2*i-1}_WHIP_CURR_DAC").astype(float)
        other = aeb.varget(prefix+f"V{2*i}_WHIP_CURR_DAC")
        measured = aeb.varget(prefix+f"V{2*i-1}_WHIP_CURR")
        native_t, native_p = native_power(hfr, antenna)
        mask = abs(t) < 600
        if not np.all(current[mask] == other[mask]):
            raise ValueError("Dipole endpoints received different bias commands")
        out["antennas"][antenna] = {
            "records_within_10_min": int(mask.sum()),
            "telemetry_flags": np.unique(aeb.varget(prefix+"TELEM_FLAG")[mask]).tolist(),
            "distinct_command_settings": int(np.unique(current[mask]).size),
            "command_range_uA": r4([current[mask].min(),current[mask].max()]),
            "median_abs_command_minus_measured_uA": round(float(np.median(abs(current[mask]-measured[mask]))),4),
            "exact_command_overlap": command_overlap(t, current),
            "packet_offset_sensitivity_s": {str(s): compare_steps(t, current, native_t, native_p, s)
                                           for s in (-7, -3.5, 0, 3.5, 7)},
        }
    print(json.dumps(out, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
