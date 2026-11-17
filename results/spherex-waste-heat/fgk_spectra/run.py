"""Extend the Gaia-matched DR1 measurements to the other 48 Contardo & Hogg stars.

The original five FGK positions and comparisons remain in ../spectra. Run the
four stages in order. A repeat of the cube stage reuses complete cached stamps.
"""
import argparse
import csv
import json
import math
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SPX = ROOT.parent / "spectra"
sys.path.insert(0, str(SPX))
import select_fields
from check_controls import retrieve
from photometry import measure, save
from read_cubes import fetch_tile


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    tmp = path.with_suffix(path.suffix + ".new")
    with tmp.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def write_json(path, value):
    tmp = path.with_suffix(path.suffix + ".new")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def retry(operation, *args):
    for attempt in range(3):
        try:
            return operation(*args)
        except Exception:
            if attempt == 2:
                raise
            time.sleep(3 * 2 ** attempt)


def select():
    catalog = read_rows(ROOT.parent / "literature/contardo_hogg.csv")
    prior = {r["source_id"] for r in read_rows(SPX / "positions_fgk.csv") if "_cmp" not in r["label"]}
    extra = [r["source_id"] for r in catalog if r["source_id"] not in prior]
    if len(set(extra)) != 48 or len(prior) != 5:
        raise ValueError("Published membership or the original five changed")
    labels = {source_id: f"FGK{i}" for i, source_id in enumerate(extra, 6)}
    query = f"SELECT {select_fields.COLUMNS} FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(extra)})"
    stars = retry(select_fields.gaia, query)
    if {r["source_id"] for r in stars} != set(extra):
        raise ValueError("Incomplete Gaia query")
    for star in stars:
        star["label"] = labels[star["source_id"]]
    with ThreadPoolExecutor(4) as pool:
        groups = list(pool.map(lambda r: retry(select_fields.field, r), stars))
    published = set(extra) | prior
    rows = [row for group in groups for row in group
            if "_cmp" not in row["label"] or row["source_id"] not in published]
    write_csv(ROOT / "positions.csv", sorted(rows, key=lambda r: r["label"]))
    write_json(ROOT / "gaia_queries.json", sorted(select_fields.RECEIPTS, key=lambda r:r["url"]))
    print(f"{len(extra)} targets and {len(rows)-len(extra)} comparison positions selected", flush=True)


def extract(cache):
    rows = read_rows(ROOT / "positions.csv")
    tiles = defaultdict(list)
    for row in rows:
        tiles[int(row["tile"])].append(row)
    cache.mkdir(parents=True, exist_ok=True)
    receipts = ROOT / "cube_receipts"
    receipts.mkdir(exist_ok=True)
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(fetch_tile, [(tile, entries, cache, receipts) for tile, entries in tiles.items()]))
    paths = [cache / (r["label"]+".npz") for r in rows]
    spectra, features, coverage, colors = measure(paths)
    if len(spectra) != 102*len(rows) or len(coverage) != len(rows):
        raise ValueError("Unexpected channel or position count")
    for name, data in (("channel_spectra",spectra),("feature_tests",features),
                       ("coverage",coverage),("continuum",colors)):
        save(ROOT / (name+".csv"),data)
    print(f"{len(tiles)} cubes, {len(rows)} positions, {len(spectra)} channels", flush=True)


def control_check():
    peers = [r for r in read_rows(ROOT / "positions.csv") if "_cmp" in r["label"]]
    # One receipt per queried star records its actual PM-adjusted search center.
    with ThreadPoolExecutor(6) as pool:
        pairs = list(pool.map(lambda r:retry(retrieve,r),peers))
    write_csv(ROOT / "comparison_wise.csv", [result for result,_ in pairs])
    write_json(ROOT / "comparison_queries.json", [receipt for _,receipt in pairs])


def compare(star, color, valid_peers, mask):
    target = color.get((star,mask),{})
    if not target.get("W2core_over_H") or not target["ratio_sky_sigma"]:
        return None,None,0
    ratio=float(target["W2core_over_H"])
    deltas,zs=[],[]
    for label in valid_peers:
        other=color.get((label,mask),{})
        if not other.get("W2core_over_H") or not other["ratio_sky_sigma"]:
            continue
        delta=ratio-float(other["W2core_over_H"])
        uncertainty=math.hypot(float(target["ratio_sky_sigma"]),float(other["ratio_sky_sigma"]))
        if uncertainty>0:
            deltas.append(delta)
            zs.append(delta/uncertainty)
    return (min(deltas),min(zs),len(zs)) if zs else (None,None,0)


def summarize(prior_controls):
    source = {r["source_id"]:r for r in read_rows(ROOT.parent/"literature/contardo_hogg.csv")}
    all_tables=[(SPX, "positions_fgk.csv"),(ROOT, "positions.csv")]
    combined=[]
    for directory,name in all_tables:
        positions=read_rows(directory/name)
        colors={(r["label"],r["extraction"]):r for r in read_rows(directory/"continuum.csv")}
        quality_file=prior_controls if directory==SPX and prior_controls else directory/"comparison_wise.csv"
        quality={r["label"]:r["status"] for r in read_rows(quality_file)}
        coverage={r["label"]:r for r in read_rows(directory/"coverage.csv")}
        for star in positions:
            label=star["label"]
            if "_cmp" in label: continue
            peers=[r["label"] for r in positions if r["label"].startswith(label+"_cmp")]
            good=[p for p in peers if quality.get(p,"").startswith(("W12_","W123_"))]
            target=source[star["source_id"]]
            row=dict(label=label,source_id=star["source_id"],distance_inv_parallax_pc=round(1000/float(star["parallax"]),2),
                     published_bb_k=target["bb_temperature_k"],published_lir_ratio=target["lir_over_lstar"],
                     strict_channels=coverage[label]["strict_ch"],outlier_included_channels=coverage[label]["available_ch"],
                     n_matched=len(peers),n_W12_quality=len(good),n_W123_quality=sum(quality.get(p)=="W123_photospheric" for p in peers))
            passes=[]
            can_test=True
            for short,mask in [("strict","strict"),("inclusive","include_OUTLIER")]:
                measured=colors[(label,mask)]
                delta,z,n=compare(label,colors,good,mask)
                row[short+"_ratio"]=measured["W2core_over_H"]
                row[short+"_sky_sigma"]=measured["ratio_sky_sigma"]
                row[short+"_n_usable_peers"]=n
                row[short+"_delta_vs_reddest_peer"]="" if delta is None else round(delta,6)
                row[short+"_min_sky_z"]="" if z is None else round(z,4)
                can_test &= n >= 2
                passes.append(z is not None and z>=3 and n>=2)
            row["continuum_check"]=("redder_than_both_in_both_masks" if all(passes) else
                                    "below_threshold" if can_test else "insufficient_controls_or_coverage")
            combined.append(row)
    if len({r["source_id"] for r in combined})!=53 or len(combined)!=53:
        raise ValueError("Combined summary does not contain all 53 published IDs")
    combined.sort(key=lambda r:int(r["label"][3:]))
    write_csv(ROOT/"all_53.csv",combined)
    counts=defaultdict(int)
    for row in combined: counts[row["continuum_check"]]+=1
    write_json(ROOT/"counts.json",dict(published_target_count=53,classification_counts=dict(counts),
                strict_channels_min=min(int(r["strict_channels"]) for r in combined),
                strict_channels_max=max(int(r["strict_channels"]) for r in combined),
                total_matched=sum(int(r["n_matched"]) for r in combined),
                clean_W12_matched=sum(int(r["n_W12_quality"]) for r in combined),
                clean_W123_matched=sum(int(r["n_W123_quality"]) for r in combined)))
    print(dict(counts),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=("select","extract","controls","summarize"))
    parser.add_argument("--cache",type=Path)
    parser.add_argument("--prior-controls",type=Path)
    args=parser.parse_args()
    if args.stage=="select": select()
    elif args.stage=="extract":
        if args.cache is None: parser.error("--cache is required for extract")
        extract(args.cache)
    elif args.stage=="controls": control_check()
    else: summarize(args.prior_controls)
