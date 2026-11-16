"""Check both releases at the 15 profile nodes plus two laser comparison nodes."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd

from profile import BASE, FIELDS, URL, WINDOWS


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def fetch_node(record):
    node = record["v17_node_id"]
    params = dict(feature="Node", feature_id=node, output="csv",
                  start_time="2023-12-12T00:00:00Z", end_time="2024-08-21T00:00:00Z",
                  collection_name="SWOT_L2_HR_RiverSP_D", fields=FIELDS)
    url = URL + "?" + urlencode(params)
    with urlopen(url, timeout=90) as response:
        raw = response.read()
    decoded = json.loads(raw)
    if decoded["status"] != "200 OK":
        raise ValueError(f"{node}: {decoded['status']}")
    data = decoded["results"]["csv"].encode()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if len(rows) != decoded["hits"] or {r["node_id"] for r in rows} != {node}:
        raise ValueError(f"Unexpected rows for {node}")
    return dict(v17_node_id=node, url=url, response_sha256=sha(raw),
                csv_sha256=sha(data), rows=len(rows)), data.decode()


def retrieve():
    translations = json.loads((BASE / "d-node-translation.json").read_text())["rows"]
    with ThreadPoolExecutor(max_workers=5) as pool:
        results = list(pool.map(fetch_node, translations))
    receipts, payloads = zip(*results)
    (BASE / "d-node-receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    # Keep the unmodified CSV payload of each query, including the original newlines.
    with (BASE / "d-node-payloads.json.gz").open("wb") as dest:
        with gzip.GzipFile(fileobj=dest, mode="wb", filename="", mtime=0) as zipped:
            zipped.write(json.dumps(payloads, separators=(",", ":")).encode())


def archived_nodes():
    receipts = json.loads((BASE / "d-node-receipts.json").read_text())
    with gzip.open(BASE / "d-node-payloads.json.gz", "rt") as src:
        payloads = json.load(src)
    tables = []
    for receipt, payload in zip(receipts, payloads, strict=True):
        if sha(payload.encode()) != receipt["csv_sha256"]:
            raise ValueError(f"CSV checksum mismatch for {receipt['v17_node_id']}")
        table = pd.read_csv(io.StringIO(payload))
        if len(table) != receipt["rows"] or {str(i) for i in table.node_id} != {receipt["v17_node_id"]}:
            raise ValueError(f"Unexpected archived rows for {receipt['v17_node_id']}")
        tables.append(table)
    return pd.concat(tables).replace(-999999999999.0, np.nan)


def summarize():
    mapping = pd.DataFrame(json.loads((BASE / "d-node-translation.json").read_text())["rows"])
    if (len(mapping) != 17 or not mapping.v16_node_id.is_unique or not mapping.v17_node_id.is_unique
            or not (mapping[["shift_flag", "boundary_flag"]] == "0").all().all()):
        raise ValueError("Expected 17 direct, unshifted, one-to-one SWORD mappings")
    mapping = mapping[["v17_node_id", "v16_node_id", "v16_reach_id"]].astype(int)
    source = archived_nodes()
    current = source.merge(mapping, left_on="node_id", right_on="v17_node_id", validate="many_to_one")
    previous = pd.read_csv(BASE / "measurements.csv").replace(-999999999999.0, np.nan)
    previous = previous.merge(mapping, left_on="node_id", right_on="v16_node_id", validate="many_to_one")
    # Independent water-surface observations at two additional nodes on 13 December.
    laser = pd.read_csv(BASE.parent.parent / "laser" / "node_matches.csv")
    laser = laser[(laser.laser_day == "2023-12-13") & (laser.swot_day == laser.laser_day)]
    day = current[(current.time_str.str[:10] == "2023-12-13") & current.pass_id.eq(570)]
    controls = laser.merge(day, left_on="node_id", right_on="v16_node_id", validate="many_to_one")
    controls["d_accepted"] = (controls.node_q.isin([0, 1]) &
                              controls.wse_u.between(0, .30, inclusive="right") & controls.wse.notna())
    controls["laser_minus_d_m"] = controls.laser_median_m - controls.wse
    v2_december = previous[previous.time_str.str[:10].eq("2023-12-13")].set_index("v16_node_id")
    lat_shift = controls.lat - controls.v16_node_id.map(v2_december.lat)
    lon_shift = controls.lon - controls.v16_node_id.map(v2_december.lon)
    controls["centroid_shift_km"] = np.hypot(lat_shift * 110.574, lon_shift * 111.32 * np.cos(np.radians(controls.lat)))
    cols = ["laser_day", "v16_node_id", "v17_node_id", "radius_km", "n_segments", "laser_median_m",
            "swot_m", "wse", "node_q", "wse_u", "d_accepted", "centroid_shift_km", "laser_minus_d_m"]
    controls[cols].rename(columns={"swot_m":"v2_wse_m","wse":"d_wse_m","node_q":"d_node_q","wse_u":"d_wse_u_m"}).to_csv(
        BASE / "d-laser-check.csv", index=False, float_format="%.5f")
    results = []
    retained = []
    for release, data in (("2.0", previous), ("D", current)):
        data = data.copy()
        data["revision"] = data.granuleUR.str.extract(r"_(\d+)\.zip$").astype(int)
        data = data.sort_values("revision").drop_duplicates(["node_id", "cycle_id", "pass_id"], keep="last")
        for label, high, low, track in (("target", WINDOWS["target_high"], WINDOWS["target_low"], 111),
                                        ("check", WINDOWS["check_high"], WINDOWS["check_low"], 570)):
            visits = data[data.pass_id.eq(track)]
            good = visits[visits.node_q.isin([0, 1]) & visits.wse_u.between(0, .30, inclusive="right") & visits.wse.notna()]
            high_rows = good[good.time_str.str[:10] == high]
            low_rows = good[good.time_str.str[:10] == low]
            both = high_rows.merge(low_rows, on="v16_node_id", suffixes=("_hi", "_lo"), validate="one_to_one")
            both["rise"] = both.wse_hi - both.wse_lo
            for reach, lo, hi in ((13217000061, 51, 55), (13217000061, 56, 60), (13217000091, 1, 5)):
                ids = mapping[mapping.v16_reach_id.eq(reach)].v16_node_id
                ids = ids[((ids // 10) % 1000).between(lo, hi)]
                group = both[both.v16_node_id.isin(ids)]
                for row in group.itertuples():
                    retained.append(dict(release=release, pair=label, v16_node_id=row.v16_node_id,
                                         v16_reach_id=reach, first_node=lo, high_wse=row.wse_hi,
                                         low_wse=row.wse_lo, rise_m=row.rise))
                results.append(dict(release=release, pair=label, pass_id=track, high=high, low=low,
                                    v16_reach_id=reach, first_node=lo, last_node=hi,
                                    n_high=high_rows.v16_node_id.isin(ids).sum(),
                                    n_low=low_rows.v16_node_id.isin(ids).sum(), n_paired=len(group),
                                    median_rise_m=group.rise.median() if len(group) else np.nan))
    pd.DataFrame(results).to_csv(BASE / "d-node-check.csv", index=False, float_format="%.5f")
    paired = pd.DataFrame(retained)
    common = paired[paired.release.eq("D")].merge(
        paired[paired.release.eq("2.0")], on=["pair", "v16_reach_id", "first_node", "v16_node_id"],
        suffixes=("_d", "_v2"), validate="one_to_one"
    ).drop(columns=["release_d", "release_v2"])
    common["d_minus_v2_high_m"] = common.high_wse_d - common.high_wse_v2
    common["d_minus_v2_low_m"] = common.low_wse_d - common.low_wse_v2
    common.to_csv(BASE / "d-node-common.csv", index=False, float_format="%.5f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retrieve", action="store_true")
    if parser.parse_args().retrieve:
        retrieve()
    summarize()
