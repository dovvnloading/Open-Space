"""Reproduce the Pool Malebo candidate screen and fixed-node checks from archived responses."""

import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, theilslopes

ROOT = Path(__file__).resolve().parent
LOWER = "13217000061"
UPPER = "13217000091"
STAGE = "13219000041"  # distinct Chenal reach well above the pool
KEY = ["cycle_id", "pass_id"]
EVENT = "2024-01-08"
LOW = "2024-07-13"


def read_feature(item):
    data = (ROOT / item["filename"]).read_bytes()
    assert hashlib.sha256(data).hexdigest() == item["sha256"]
    table = pd.read_csv(io.StringIO(json.loads(data)["results"]["csv"])).replace(
        -999999999999.0, np.nan
    )
    assert len(table) == item["rows"]
    assert table.sword_version.eq(16).all()
    table["date"] = pd.to_datetime(
        table.time_str, format="%Y-%m-%dT%H:%M:%SZ", errors="coerce", utc=True
    ).dt.strftime("%Y-%m-%d")
    return table


def screen(table):
    return table[
        table.date.notna()
        & table.reach_q.isin([0, 1])
        & (table.partial_f == 0)
        & (table.obs_frac_n >= 0.8)
        & table.wse_u.between(0, 0.30)
        & table.wse.notna()
    ].copy()


def block_interval(frame, x, y, repetitions=1500, block=3):
    # Consecutive observed passes in blocks. Gap lengths differ; this is a sensitivity interval.
    xval, yval = frame[x].to_numpy(), frame[y].to_numpy()
    n = len(frame)
    rng = np.random.default_rng(13061)
    boot = []
    for _ in range(repetitions):
        starts = rng.integers(0, n, (n + block - 1) // block)
        idx = ((starts[:, None] + np.arange(block)) % n).ravel()[:n]
        if len(np.unique(xval[idx])) >= 2:
            boot.append(theilslopes(yval[idx], xval[idx]).slope)
    return np.quantile(boot, [0.025, 0.975]).tolist()


def spherical_m(lat1, lon1, lat2, lon2):
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    term = (
        np.sin(dlat / 2) ** 2
        + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2) ** 2
    )
    return 2 * 6371008.8 * np.arcsin(np.sqrt(term))


def main():
    reaches = {
        Path(i["filename"]).stem: read_feature(i)
        for i in json.loads((ROOT / "manifest.json").read_text())
    }
    screened = {rid: screen(table) for rid, table in reaches.items()}
    counts = []
    for rid, table in reaches.items():
        good = screened[rid]
        counts.append(
            dict(
                reach_id=rid,
                raw=len(table),
                dated=table.date.notna().sum(),
                valid_consensus=table.sos_consensus_q.notna().sum(),
                screened=len(good),
                screened_consensus=good.sos_consensus_q.notna().sum(),
                good_flag=(good.reach_q == 0).sum(),
            )
        )
    pd.DataFrame(counts).to_csv(ROOT / "coverage.csv", index=False)
    lower = screened[LOWER]
    upper = screened[UPPER]
    pair = lower[
        KEY
        + [
            "date",
            "wse",
            "wse_u",
            "slope",
            "width",
            "sos_consensus_q",
            "loc_offset",
            "crid",
        ]
    ].merge(
        upper[KEY + ["wse", "wse_u", "loc_offset", "crid"]],
        on=KEY,
        suffixes=("_lower", "_upper"),
        validate="one_to_one",
    )
    assert pair.crid_lower.eq(pair.crid_upper).all()
    pair["drop_m"] = pair.wse_upper - pair.wse_lower
    # Reach 13219000041 is a separate stage indicator on track 111 only. Its observations share the sensor.
    pair = pair.merge(
        screened[STAGE][KEY + ["wse"]].rename(columns={"wse": "chenal_wse"}),
        on=KEY,
        how="left",
        validate="one_to_one",
    )
    pair.sort_values("date").to_csv(ROOT / "paired-reach.csv", index=False)
    stats = []
    for (ps, crid), rows in list(pair.groupby(["pass_id", "crid_lower"])) + [
        ((ps, "all"), rows) for ps, rows in pair.groupby("pass_id")
    ]:
        for predictor in ["chenal_wse", "wse_upper"]:
            g = rows[["date", predictor, "drop_m"]].dropna().sort_values("date")
            if len(g) < 3:
                continue
            x = g[predictor]
            y = g.drop_m
            ci = (
                block_interval(g, predictor, "drop_m")
                if len(g) >= 8
                else [np.nan, np.nan]
            )
            stats.append(
                dict(
                    pass_id=ps,
                    crid=crid,
                    predictor=predictor,
                    n=len(g),
                    spearman_r=spearmanr(x, y).statistic,
                    theil_drop_per_m=theilslopes(y, x).slope,
                    block95_low=ci[0],
                    block95_high=ci[1],
                )
            )
    pd.DataFrame(stats).to_csv(ROOT / "rank-and-slope.csv", index=False)
    # Assess departures from the relationship fitted without the independently dated flood window.
    data = (
        pair[(pair.pass_id == 111) & pair.chenal_wse.notna()].sort_values("date").copy()
    )
    flood = data.date.between("2023-12-01", "2024-01-31")
    calibration = data[~flood]
    fit = theilslopes(calibration.drop_m, calibration.chenal_wse, method="joint")
    data["event"] = flood
    data["predicted_drop_m"] = fit.intercept + fit.slope * data.chenal_wse
    data["residual_m"] = data.drop_m - data.predicted_drop_m
    data["outside_calibration_range"] = ~data.chenal_wse.between(
        calibration.chenal_wse.min(), calibration.chenal_wse.max()
    )
    data[
        [
            "date",
            "event",
            "chenal_wse",
            "drop_m",
            "predicted_drop_m",
            "residual_m",
            "outside_calibration_range",
        ]
    ].to_csv(ROOT / "event-residuals.csv", index=False)
    node_rows = []
    for i in json.loads((ROOT / "node-manifest.json").read_text()):
        d = read_feature(i)
        g = d[d.date.isin([EVENT, LOW]) & d.node_q.isin([0, 1])].set_index("date")
        assert len(g) == 2
        a, b = g.loc[EVENT], g.loc[LOW]
        node_rows.append(
            dict(
                node_id=int(a.node_id),
                lat_high=a.lat,
                lon_high=a.lon,
                lat_low=b.lat,
                lon_low=b.lon,
                wse_high=a.wse,
                wse_low=b.wse,
                delta_m=a.wse - b.wse,
                u_high=a.wse_u,
                u_low=b.wse_u,
                displacement_m=spherical_m(a.lat, a.lon, b.lat, b.lon),
            )
        )
    pd.DataFrame(node_rows).to_csv(ROOT / "fixed-nodes.csv", index=False)
    print(pd.DataFrame(counts).to_string(index=False))
    print(pd.DataFrame(stats).round(3).to_string(index=False))
    print(pd.DataFrame(node_rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
