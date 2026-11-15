"""Test differential Na emission at Jupiter around Io's independently dated eruption."""
import argparse
import hashlib
import json
import warnings
from pathlib import Path
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RECORD = "https://zenodo.org/records/21420670/files/"
FILES = {"IoIO_Na_nebula_2017-2026.ecsv": "397ad9c6f6e07c3b1ceb4b81f1338120",
         "IoIO_Torus_2017-2026.ecsv": "293dfce1cd1542e8a113d4ff70d7c362"}
ROOT = Path(__file__).resolve().parent
OCT = pd.Timestamp("2024-10-22")
DEC = pd.Timestamp("2024-12-27")
END = DEC + pd.Timedelta(days=112)  # exclusive
DAY = pd.Timedelta(days=1)
COLUMNS = ["Na3", "Na6", "Na12", "Na24", "dawn", "dusk"]
WINDOWS = [("reference", OCT-42*DAY, OCT), ("onset bracket", OCT, DEC)] + [
    (f"PJ68 +{start} to +{start+27} d", DEC+start*DAY, DEC+(start+28)*DAY)
    for start in range(0,112,28)]


def read_release(cache):
    tables = []
    for filename, expected in FILES.items():
        path = cache / filename
        if not path.exists():
            cache.mkdir(parents=True, exist_ok=True)
            with urlopen(RECORD+filename, timeout=90) as stream:
                content = stream.read()
            if hashlib.md5(content).hexdigest() != expected:
                raise ValueError(f"Unexpected downloaded bytes for {filename}")
            path.write_bytes(content)
        if hashlib.md5(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Checksum mismatch for {path}")
        t = pd.read_csv(path, sep=" ", comment="#")
        t = t.loc[t["mask"].ne(True)].copy()  # pre-2023 missing masks retained
        t["date"] = pd.to_datetime(t.tavg, utc=True).dt.tz_localize(None).dt.normalize()
        if "sun_angle" in t:
            for rad in (3, 6, 12, 24):
                t[f"Na{rad}"] = t[f"largest_sub_annular_sb_{rad:.1f}_jupiterRad"]
        else:
            t = t.rename(columns={"ansa_left_surf_bright":"dawn","ansa_right_surf_bright":"dusk"})
        tables.append(t)
    return tables


def summarize(n, t):
    cols = COLUMNS[:4]+["sun_angle","airmass","rayleigh_conversion"]
    return n.groupby("date")[cols].median().join(
        t.groupby("date")[["dawn","dusk"]].median(), how="outer").sort_index()


def cut(daily, start, stop):
    return daily.loc[(daily.index >= start) & (daily.index < stop)]


def window_table(daily):
    records = []
    for label, start, stop in WINDOWS:
        d = cut(daily, start, stop)
        for channel in COLUMNS:
            x = d[channel].dropna()
            records.append(dict(window=label, start=start.date().isoformat(),
                                end_inclusive=(stop-DAY).date().isoformat(), channel=channel,
                                nights=x.size, median_R=x.median(), q25_R=x.quantile(.25),
                                q75_R=x.quantile(.75)))
    return pd.DataFrame(records)


def scan(series, first, end, nday):
    # One trial for EVERY calendar start. Slice in days, not in observation count.
    continuous = series.reindex(pd.date_range(first,end-DAY))
    slices = np.lib.stride_tricks.sliding_window_view(continuous.values, nday)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(slices, axis=1)
    return pd.DataFrame({"nights":np.isfinite(slices).sum(axis=1),"median":med},
                        index=continuous.index[:len(med)])


def historical_controls(daily):
    na = daily.dropna(subset=["Na6"])
    season = (na.index.to_series().diff() > 75*DAY).cumsum()
    rows = []
    for _, d in na.groupby(season):
        opp = d.sun_angle.idxmax()
        # Offset fixed by observed 2024 season's maximum elongation (Dec 8).
        anchor = opp + 19*DAY
        pre = cut(d, anchor-42*DAY, anchor)
        tail = cut(d, anchor+28*DAY, anchor+70*DAY)
        if d.sun_angle.max()<170 or min(len(pre),len(tail))<5:
            continue
        rows.append({"anchor":str(anchor.date()), "pre_nights":len(pre), "late_nights":len(tail),
                     "pre_R":pre.Na6.median(),"late_R":tail.Na6.median(),
                     "relative":tail.Na6.median()/pre.Na6.median()})
    return pd.DataFrame(rows)


def checks(n, t, daily):
    out = {}
    series = daily.Na6.dropna()
    base = cut(series, OCT-42*DAY, OCT).median()
    for width in (14,21,28,42,56):
        win = scan(series, OCT, END, width)
        out[f"slide_{width}d"]={
            "n_windows":len(win), "min_nights":int(win.nights.min()),
            "maximum_R":win["median"].max(), "minimum_R":win["median"].min(),
            "peak_start":str(win["median"].idxmax().date())}
    # Calendar-aligned injection into already observed event-season measurements.
    # A constant multiplicative increment that lasts exactly one 28-day window
    # scales its nightly median by that factor. Other overlapping windows can only
    # increase the recovery fraction. Missing dates are never filled.
    win = scan(series,OCT,END,28)
    out["injection"]={str(factor):{"trials":len(win),"recovered_above_1.5x_baseline":
                       int(((win["median"]*factor)>1.5*base).sum())}
                      for factor in (1.,1.5,2.,2.5,3.)}
    out["full_recovery_multiplier_at_1.5x_threshold"]=1.5*base/win["median"].min()
    s = scan(series,pd.Timestamp("2022-08-01"),pd.Timestamp("2022-12-01"),28)
    s = s.loc[s.nights>=7]
    out["2022_positive_control"]={"peak_R":s["median"].max(),
             "peak_start":str(s["median"].idxmax().date()),
             "nights":int(s.loc[s["median"].idxmax(),"nights"])}
    in_range = cut(series,OCT,END)
    dates = pd.Series(in_range.index,index=in_range.index)
    gap = dates.diff().idxmax()
    out["largest_sampling_gap"]={"last_before":str(dates.shift().loc[gap].date()),
              "first_after":str(gap.date()),"empty_dates":int((dates.diff().loc[gap]/DAY)-1)}
    out["airmass_sensitivity"]={}
    for limit in (1.3,1.5,2.,3.):
        # Restrict individual exposures before median aggregation.
        good = n.loc[n.airmass <= limit].groupby("date").Na6.median()
        ref = cut(good,OCT-42*DAY,OCT); late = cut(good,DEC+28*DAY,DEC+84*DAY)
        out["airmass_sensitivity"][str(limit)]={"reference_R":ref.median(),"late_R":late.median(),
                                                "n_ref":len(ref),"n_late":len(late)}
    out["elongation_matched"]={}
    for lower,upper in [(80,110),(110,140),(140,165)]:
        mask=(daily.sun_angle>=lower)&(daily.sun_angle<upper)
        a=cut(daily[mask],"2024-08-01","2024-12-07").dropna(subset=["Na6"])
        b=cut(daily[mask],"2024-12-27",END).dropna(subset=["Na6"])
        out["elongation_matched"][f"{lower}-{upper}"]={
            "pre_R":a.Na6.median(),"later_R":b.Na6.median(),
            "n_pre":len(a),"n_later":len(b)}
    # Test the sodium trend within fixed System III longitude bins at low airmass.
    geometry = n.loc[n.airmass.le(1.5)].copy()
    geometry["lonbin"] = (geometry.Jupiter_PDObsLon % 360 // 45).astype(int)
    matched = []
    for bin_id, group in geometry.groupby("lonbin"):
        by_day = group.groupby("date").Na6.median()
        a, b = cut(by_day, OCT-42*DAY, OCT), cut(by_day, DEC+28*DAY, DEC+84*DAY)
        if min(len(a),len(b)) >= 3:
            matched.append({"longitude_start_deg":int(bin_id*45), "n_reference":len(a),
                            "n_late":len(b), "reference_R":a.median(),"late_R":b.median()})
    out["Na_by_viewing_longitude_at_airmass_le_1.5"] = matched
    # Independent sulfur-line chronology; November is exploratory. Match viewing
    # longitude for each ansa and re-aggregate exposures within each night/bin.
    nov = cut(daily,"2024-11-01","2024-11-29")
    out["november_1_28"] = {col:{"n":int(nov[col].count()),"R":nov[col].median()}
                            for col in COLUMNS}
    out["same_night_three_channels"] = {}
    for label, subset in [("reference",cut(daily,OCT-42*DAY,OCT)),("november",nov)]:
        paired = subset[["Na6","dawn","dusk"]].dropna()
        out["same_night_three_channels"][label] = {"n":len(paired), **paired.median().to_dict(),
               "median_dawn_over_Na6":(paired.dawn/paired.Na6).median(),
               "median_dusk_over_Na6":(paired.dusk/paired.Na6).median()}
    sulfur=t.loc[t.airmass.le(1.5)].copy()
    sulfur["lonbin"] = (sulfur.Jupiter_PDObsLon%360//45).astype(int)
    out["november_torus_by_longitude"] = {}
    for channel in ("dawn","dusk"):
        matches = []
        for bin_id, group in sulfur.groupby("lonbin"):
            by_day=group.groupby("date")[channel].median().dropna()
            ref=cut(by_day,OCT-42*DAY,OCT)
            peak=cut(by_day,"2024-11-01","2024-11-29")
            if min(len(ref),len(peak))>=3:
                matches.append({"longitude_start_deg":int(bin_id*45),
                                "reference_n":len(ref),"november_n":len(peak),
                                "difference_R":peak.median()-ref.median()})
        out["november_torus_by_longitude"][channel] = matches
    # Equal-weight resampling of seven-day UTC blocks preserves correlations
    # within a week, including uneven numbers of nightly data points.
    def resample_blocks(window, start, rng):
        block_id = np.asarray((window.index-start)/DAY, dtype=int) // 7
        blocks = [v.to_numpy() for _,v in window.groupby(block_id)]
        indices = rng.integers(len(blocks), size=(4000,len(blocks)))
        return np.array([np.median(np.concatenate([blocks[j] for j in row])) for row in indices])
    ref = cut(series,OCT-42*DAY,OCT)
    later = cut(series,DEC+28*DAY,DEC+84*DAY)
    rng = np.random.default_rng(27122024)
    boot = resample_blocks(later, DEC+28*DAY, rng)/resample_blocks(ref,OCT-42*DAY,rng)
    out["late_over_reference"] = {"ratio":later.median()/ref.median(),
        "week_block_bootstrap_95":np.quantile(boot, [.025,.975]).tolist(),
        "reference_nights":len(ref), "late_nights":len(later)}
    return out


def figure(daily, path):
    start = OCT-42*DAY
    d = cut(daily,start,END)
    fig,axs=plt.subplots(2,1,figsize=(9,5.5),sharex=True,layout="constrained")
    for ax, names, colors, labels in [
        (axs[0],["Na6","Na12"],["#146c80","#ab5926"],["6 Jupiter radii","12 Jupiter radii"]),
        (axs[1],["dawn","dusk"],["#614b87","#157458"],["Dawn","Dusk"])]:
        for col,color,label in zip(names,colors,labels):
            ax.scatter(d.index,d[col],color=color,s=8,alpha=.35)
            # Display consecutive non-overlapping 14-day bins with >=3 nights.
            gr=d[col].resample("14D",origin=DEC); meds=gr.median().where(gr.count()>=3)
            ax.plot(meds.index+7*DAY,meds,color=color,lw=2,marker="o",ms=3,label=label)
        ax.axvspan(OCT,DEC,color="#a4b1bd",alpha=.17)
        ax.axvline(DEC,color="#333333",lw=1,ls="--")
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(loc="upper right",ncols=2,frameon=False)
        ax.grid(axis="y",alpha=.2)
    axs[0].set_title("Jovian sodium nebula during Io’s 2024 southern eruption",loc="left")
    axs[0].set_ylabel("Na brightness (R)")
    axs[1].set_ylabel("[S II] brightness (R)")
    axs[1].set_xlabel("UTC date   •   Shading brackets eruption onset between Juno visits")
    fig.savefig(path,dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=ROOT/".cache")
    parser.add_argument("--out",type=Path,default=ROOT)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    n,t=read_release(args.cache)
    daily=summarize(n,t)
    daily.to_csv(args.out/"nightly.csv",float_format="%.4f",date_format="%Y-%m-%d")
    window_table(daily).to_csv(args.out/"windows.csv",index=False,float_format="%.4f")
    historical_controls(daily).to_csv(args.out/"historical.csv",index=False,float_format="%.4f")
    result=checks(n,t,daily)
    (args.out/"checks.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    figure(daily,args.out/"brightness.png")


if __name__=="__main__":
    main()
