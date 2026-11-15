"""Reproduce the March 2018 GLEAM-X J1627 pulse-spacing measurement."""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from scipy.stats import rankdata, spearmanr


HERE = Path(__file__).resolve().parent
SOURCE = "https://raw.githubusercontent.com/nhurleywalker/GLEAM-X_Periodic_Transient/653038358ec146f4eba7128e3c97d319bd55ca72/"
PERIOD = 1091.1708  # authors' zero-period-derivative stack, seconds
GPS = datetime(1980, 1, 6, tzinfo=timezone.utc)


def inputs(cache):
    cache.mkdir(parents=True, exist_ok=True)
    paths = []
    for remote, digest in json.loads((HERE / "inputs.json").read_text()).items():
        path = cache / Path(remote).name
        if not path.exists():
            with urlopen(SOURCE + remote, timeout=60) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != digest:
                raise ValueError(f"Unexpected upstream content at {remote}")
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected cached content at {path}")
        paths.append(path)
    return paths


def fit_pair(x, flux, flat_base=False):
    height = flux.max()
    def residual(p):
        a, left, w, b, right, s, base, grad = p
        leading = a * np.exp(-0.5 * ((x - left) / w) ** 2)
        trailing = b * np.exp(-0.5 * ((x - right) / s) ** 2)
        return leading + trailing + base + grad * (x - 15) - flux
    smoothed = gaussian_filter1d(flux, 1.5)
    start = [0.7 * height, 12, 2, 0.4 * height, 19, 2, smoothed.min(), 0]
    grad_bound = height * 0.2 if not flat_base else 1e-12
    lower = [0, 9, 0.5, 0, 16, 0.5, -height, -grad_bound]
    upper = [2*height, 15, 5, 2*height, 22, 7, height, grad_bound]
    start = np.clip(start, lower, upper)
    p = least_squares(residual, start, bounds=(lower, upper)).x
    return p[1], p[4], p[4] - p[1]


def peaks(x, flux):
    grid = np.arange(7, 23.01, 0.05)
    smooth = gaussian_filter1d(np.interp(grid, x, flux), 15)
    positions = []
    for lo, hi in [(8, 16), (16, 22)]:
        inside = (grid > lo) & (grid < hi)
        positions.append(grid[inside][smooth[inside].argmax()])
    between = (grid>=positions[0]) & (grid<=positions[1])
    dip = min(np.interp(positions,grid,smooth))-smooth[between].min()
    return positions[1] - positions[0], dip


def measure(path):
    obsid = int(path.name[:10])
    time = GPS + timedelta(seconds=obsid-18)  # GPS - UTC in 2018
    frequency = float(re.search(r"frequency ([\d.]+) MHz", path.read_text()).group(1))
    row = dict(obsid=obsid, utc_start=time.strftime("%Y-%m-%dT%H:%M:%SZ"), MHz=frequency)
    t, flux = np.loadtxt(path).T
    phase = (obsid + t + PERIOD/2) % PERIOD - 0.48*PERIOD
    finite = np.isfinite(flux)
    off = finite & ((phase < -35) | (phase > 33))
    on = finite & (phase > 4) & (phase < 27)
    row.update(on_samples=int(on.sum()), off_samples=int(off.sum()), status="coverage")
    if off.sum() < 6 or on.sum() < 37:
        return row
    baseline = np.median(flux[off])
    signal = flux[on] - baseline
    noise = 1.4826 * np.median(np.abs(flux[off] - baseline))
    row["peak_snr"] = signal.max()/noise if noise > 0 else 0
    if row["peak_snr"] < 5:
        row["status"] = "low_snr"
        return row
    x = phase[on]
    order = np.argsort(x)
    x, signal = x[order], signal[order]
    if time.date().isoformat() < "2018-03-11" or time.date().isoformat() > "2018-03-16":
        row["status"] = "outside_comparison_interval"
        return row
    left, right, gap = fit_pair(x, signal)
    direct_gap, dip = peaks(x,signal)
    row.update(status="measured", leading_s=left, trailing_s=right, spacing_s=gap,
               flat_baseline_spacing_s=fit_pair(x, signal, flat_base=True)[2],
               smoothed_peak_spacing_s=direct_gap, dip_over_native_noise=dip/noise)
    return row


def quarters(cache):
    out = []
    for band, obsid in [(85,1205008192),(119,1205007112),(154,1205011432),
                         (185,1205010352),(215,1205009272)]:
        a = np.loadtxt(cache / f"{band}MHz_dedisper.txt")
        nu = np.loadtxt(cache / f"{band}MHz_freq.txt") * 1e-6
        good = np.isfinite(a).all(axis=1) & (np.ptp(a, axis=1)>0)
        a, nu = a[good], nu[good]
        # Both gates have complete samples across all five files. Other
        # parts of the dynamic spectra contain zeros as padding.
        if np.any(a[:,64:110] == 0) or np.any(a[:,112:127] == 0):
            raise ValueError("Missing samples in the common spectral gates")
        a -= np.median(a[:,112:127],axis=1)[:,None]
        # The aligned release uses a different arbitrary origin from the profile
        # files. The 28 s translation centers the two components in the fit.
        phase = 0.5*np.arange(a.shape[1])-28
        gate = (phase>3.8) & (phase<26.8)
        for j, indices in enumerate(np.array_split(np.arange(len(a)),4), 1):
            y = a[indices].mean(axis=0)
            left,right,gap = fit_pair(phase[gate],y[gate])
            out.append(dict(obsid=obsid,band_file=band,quarter=j,channels=len(indices),
                            smoothed_peak_spacing_s=peaks(phase[gate],y[gate])[0],
                            MHz=nu[indices].mean(),leading_s=left,trailing_s=right,spacing_s=gap))
    return out


def pretty(value):
    if isinstance(value, (np.floating, float)):
        return round(float(value), 6)
    return value


def save_csv(p, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with p.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=keys)
        writer.writeheader()
        writer.writerows({k:pretty(v) for k,v in row.items()} for row in rows)


def rank_permutation(nu, gap):
    # Ordering by spacing reduces the exact Spearman test to integer dot
    # products of the assigned tied frequency ranks with consecutive ranks.
    idx = np.argsort(gap)
    frequency_ranks = (2*rankdata(nu)).astype(int)
    spacing_ranks = 2*np.arange(1,len(gap)+1)-(len(gap)+1)
    ranks = frequency_ranks-(len(gap)+1)
    threshold = abs(int(spacing_ranks @ ranks[idx]))
    counts = Counter(ranks.tolist())
    values = list(counts)
    def count(i, total):
        if i == len(gap):
            return int(abs(total) >= threshold)
        extreme = 0
        for rank in values:
            if counts[rank]:
                counts[rank] -= 1
                extreme += count(i+1,total+int(spacing_ranks[i])*rank)
                counts[rank] += 1
        return extreme
    extreme = count(0,0)
    n = math.factorial(len(gap)) // math.prod(math.factorial(v) for v in counts.values())
    return dict(extreme=extreme,assignments=n,p_two_sided=extreme/n)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--out",type=Path,default=HERE)
    args=parser.parse_args()
    paths=inputs(args.cache)
    rows=[measure(p) for p in paths if p.suffix==".dat"]
    q=quarters(args.cache)
    core=[r for r in rows if r["utc_start"].startswith("2018-03-13") and r["status"]=="measured"]
    nu=np.array([r["MHz"] for r in core])
    spacing=np.array([r["spacing_s"] for r in core])
    ids=np.array([r["obsid"] for r in core])
    design=np.column_stack([np.ones(len(ids)),np.log(nu/150),(ids-ids.mean())/3600])
    result=dict(screened_profiles=len(rows),march13_count=len(core),band_medians={})
    for freq in np.unique(nu):
        members=[r for r in core if r["MHz"]==freq]
        result["band_medians"][f"{freq:.2f}"]=dict(n=len(members),**{
            col:pretty(np.median([r[col] for r in members])) for col in (
                "leading_s","trailing_s","spacing_s","flat_baseline_spacing_s","smoothed_peak_spacing_s")})
    for col in ("leading_s","trailing_s","spacing_s","flat_baseline_spacing_s","smoothed_peak_spacing_s"):
        vals=[r[col] for r in core]
        coefs=np.linalg.lstsq(design,vals,rcond=None)[0]
        result[col]=dict(spearman=pretty(spearmanr(nu,vals).statistic),
                         log_frequency_coefficient=pretty(coefs[1]),per_hour_coefficient=pretty(coefs[2]))
    fig_rotations={r["obsid"] for r in q}
    for key, selected in (
        ("five_dynamic_spectrum_rotations",np.array([i in fig_rotations for i in ids])),
        ("other_six_rotations",np.array([i not in fig_rotations for i in ids]))):
        subset=design[selected]
        fit=np.linalg.lstsq(subset,spacing[selected],rcond=None)[0]
        quadratic=np.column_stack([subset,subset[:,2]**2])
        result[key]=dict(n=int(selected.sum()),rho=pretty(spearmanr(nu[selected],spacing[selected]).statistic),
                         slope_log_MHz_linear_time=pretty(fit[1]),
                         slope_log_MHz_quadratic_time=pretty(np.linalg.lstsq(quadratic,spacing[selected],rcond=None)[0][1]))
    loo=[]
    for omitted in range(len(ids)):
        keep=np.arange(len(ids))!=omitted
        loo.append(np.linalg.lstsq(design[keep],spacing[keep],rcond=None)[0][1])
    result["leave_one_rotation_out_log_MHz_coefficient"]=[pretty(min(loo)),pretty(max(loo))]
    result["rank_exchangeability_check"]=rank_permutation(nu,spacing)
    higher = nu > 100
    result["omitting_88_MHz"] = dict(
        n=int(higher.sum()), rho=pretty(spearmanr(nu[higher],spacing[higher]).statistic),
        **rank_permutation(nu[higher],spacing[higher]))
    resolved=np.array([r["dip_over_native_noise"]>=2 for r in core])
    direct=[r["smoothed_peak_spacing_s"] for r in core]
    result["valley_at_least_twice_native_noise"]=dict(n=int(resolved.sum()),
        fitted_rho=pretty(spearmanr(nu[resolved],spacing[resolved]).statistic),
        smoothed_maxima_rho=pretty(spearmanr(nu[resolved],np.array(direct)[resolved]).statistic))
    observed=abs(spearmanr(nu,spacing).statistic)
    shift_rhos=[spearmanr(nu,np.roll(spacing,i)).statistic for i in range(len(nu))]
    result["cyclic_alignment"] = dict(
        as_extreme=sum(bool(abs(x) >= observed - 1e-12) for x in shift_rhos),
        possible=len(shift_rhos), correlations=[pretty(v) for v in shift_rhos])
    result["adjacent_band_changes"] = dict(
        opposite=int(np.sum(np.diff(nu)*np.diff(spacing)<0)),
        total=int(np.count_nonzero(np.diff(nu))))
    args.out.mkdir(parents=True,exist_ok=True)
    save_csv(args.out/"profiles.csv",rows)
    save_csv(args.out/"simultaneous_quarters.csv",q)
    (args.out/"summary.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
