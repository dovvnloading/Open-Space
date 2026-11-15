"""Measure Perseverance sky radiometry at Phobos conjunctions and at shifted epochs."""
import argparse
import hashlib
from pathlib import Path
import urllib.request
import warnings

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CHANNELS = ["RDS_TOP_7"] + [f"RDS_LAT_{i}" for i in range(2, 9)]
AXIS = np.arange(-90, 91)
CENTRAL = np.abs(AXIS[2:-2]) <= 40


def fetch(source, cache):
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / source["product"]
    if not path.exists():
        with urllib.request.urlopen(source.url, timeout=90) as response:
            payload = response.read()
        if hashlib.sha256(payload).hexdigest() != source.sha256:
            raise ValueError(f"Unexpected archive version: {source['product']}")
        path.write_bytes(payload)
    elif hashlib.sha256(path.read_bytes()).hexdigest() != source.sha256:
        raise ValueError(f"Corrupt or updated archive product: {source['product']}")
    return path


def read_series(source, cache):
    raw = pd.read_csv(fetch(source, cache), usecols=["SCLK", *CHANNELS])
    raw = raw.drop_duplicates("SCLK").set_index("SCLK").sort_index()[CHANNELS]
    return raw.reindex(np.arange(raw.index.min(), raw.index.max() + 1))


def measure(values):
    if values.shape != (181, len(CHANNELS)):
        return None
    valid = np.isfinite(values).all(axis=1)
    flanks = (np.abs(AXIS) >= 45) & valid
    if (not valid[np.abs(AXIS) <= 42].all()
            or min(sum(flanks & (AXIS < 0)), sum(flanks & (AXIS > 0))) < 12
            or min(values[valid, 0]) < 10):
        return None
    slope, intercept = np.polyfit(AXIS[flanks], values[flanks], 1)
    baseline = AXIS[:, None] * slope + intercept
    deficits = 100 * (1 - values / baseline)
    smooth = np.lib.stride_tricks.sliding_window_view(deficits, 5, axis=0).mean(axis=-1)[CENTRAL]
    # Contemporaneous median of seven independent directions suppresses common darkening.
    contrast = smooth[:, 1:] - np.median(smooth[:, 1:], axis=1)[:, None]
    peaks = contrast.max(axis=0)
    # Only genuinely adjacent sensors count. The missing LAT1 leaves LAT8 and LAT2 separated by 90 deg.
    pair = np.minimum(peaks[:-1], peaks[1:])
    strongest = int(np.argmax(pair))
    return dict(score=float(pair[strongest]), pair=f"{strongest+2}/{strongest+3}",
                absolute=smooth.max(axis=0), times=smooth.argmax(axis=0)-40,
                contrast=peaks, contrast_times=contrast.argmax(axis=0)-40)


def window(series, center):
    return series.reindex(center + AXIS).to_numpy()


def screen(cache, output):
    events = pd.read_csv(HERE / "geometry.csv")
    products = pd.read_csv(HERE / "sources.csv").set_index("sol")
    outcomes, controls, detail, fragments = [], [], [], []
    for sol, subset in events.groupby("sol"):
        source = products.loc[sol]
        available = pd.notna(source["url"])
        data = read_series(source, cache) if available else None
        for event in subset.itertuples():
            center = event.sclk
            measurement = measure(window(data, center)) if available else None
            row = dict(sol=sol, utc=event.utc, geometry=event.geometry)
            if measurement is None:
                row["coverage"] = "no_listed_RDS_product" if not available else "insufficient_daylight_or_sampling"
                outcomes.append(row)
                continue
            comparison = []
            # The event and trial windows share a sol and grid phase; trial windows never overlap.
            for t in range(center - 181*((center-data.index.min())//181), data.index.max()+1, 181):
                if min(abs(t - other) for other in subset.sclk) <= 181:
                    continue
                score = measure(window(data, t))
                if score is not None:
                    comparison.append(score["score"])
                    controls.append(dict(sol=sol,event_sclk=center,control_sclk=t,score_pct=score["score"]))
            row.update(coverage="usable", adjacent_pair=measurement["pair"],
                       score_pct=measurement["score"], controls=len(comparison),
                       controls_ge_event=sum(x >= measurement["score"] for x in comparison))
            outcomes.append(row)
            for idx, channel in enumerate(CHANNELS):
                detail.append(dict(sol=sol,utc=event.utc,channel=channel,
                    depth_pct=measurement["absolute"][idx],peak_offset_s=measurement["times"][idx],
                    contrast_pct=measurement["contrast"][idx-1] if idx else np.nan,
                    contrast_peak_offset_s=measurement["contrast_times"][idx-1] if idx else np.nan))
            if sol in [401,403,622,1061,1062,1722]:
                block = data.reindex(center + AXIS).reset_index(names="sclk")
                block.insert(0,"offset_s",AXIS); block.insert(0,"sol",sol);fragments.append(block)
    for name, data in [("outcomes",pd.DataFrame(outcomes)), ("channel_depths",pd.DataFrame(detail)),
                       ("control_scores",pd.DataFrame(controls)),("selected_windows",pd.concat(fragments))]:
        data.to_csv(output / (name+".csv"),index=False,float_format="%.9g")


def verify_windows():
    samples = pd.read_csv(HERE/"selected_windows.csv")
    expected = pd.read_csv(HERE/"outcomes.csv").set_index("sol")
    for sol, block in samples.groupby("sol"):
        measured = measure(block.sort_values("offset_s")[CHANNELS].to_numpy())
        row = expected.loc[sol]
        if isinstance(row,pd.DataFrame):
            row = row.iloc[0]
        assert np.isclose(measured["score"],row.score_pct,rtol=5e-6)
        assert measured["pair"] == row.adjacent_pair
    print(f"Verified {samples.sol.nunique()} event windows")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full",action="store_true",help="Rebuild all measurements and controls from NASA PDS")
    parser.add_argument("--cache",type=Path,default=Path("meda-cache"))
    parser.add_argument("--output",type=Path,default=Path("."))
    args=parser.parse_args()
    if args.full:
        args.output.mkdir(parents=True,exist_ok=True)
        screen(args.cache,args.output)
    else:
        verify_windows()
