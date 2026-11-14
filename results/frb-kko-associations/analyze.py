"""Reproduce the Cat2/KKO crossmatch and its right-ascension controls."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from matplotlib.path import Path as Polygon

ROOT = Path(__file__).resolve().parent
BASE = "https://assets.chime-frb.ca/catalogs/"
CAT = BASE + "catalog2/releases/20251231/data/chimefrbcat2.csv"
KKO = BASE + "outriggers1/releases/2502.11217v2/data/outriggers1.csv"
MAPS = ("https://www.canfar.net/storage/vault/file/AstroDataCitationDOI/"
        "CISTI.CANFAR/25.0066/data/localizations/data/")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def number(value):
    return float(value) if value else math.nan


def read_input(url, name, cache, hashes):
    path = cache / name
    if not path.exists():
        with urlopen(url, timeout=120) as response:
            data = response.read()
        if digest(data) != hashes[name]:
            raise ValueError(f"Archive checksum mismatch for {name}")
        path.write_bytes(data)
    data = path.read_bytes()
    if digest(data) != hashes[name]:
        raise ValueError(f"Checksum mismatch for {path}")
    return data


def build_inputs(cache, hashes):
    import h5py

    cat = list(csv.DictReader(io.StringIO(read_input(CAT, "chime.csv", cache, hashes).decode())))
    kko = list(csv.DictReader(io.StringIO(read_input(KKO, "outriggers.csv", cache, hashes).decode())))
    events = {r["tns_name"]: r for r in cat if r["sub_num"] == "0"}
    if len(events) != len({r["tns_name"] for r in cat}):
        raise ValueError("Each event must have one component zero")
    early = [r for r in events.values() if not r["repeater_name"]
             and r["sidelobe_flag"] == "0" and r["excluded_flag"] == "0"
             and all(math.isfinite(number(r[k])) for k in ("ra", "dec", "dm_fitb"))]
    anchors = [a for a in kko if not a["repeater_name"]]
    pairs = []
    for a in anchors:
        counterpart = events.get(a["tns_name"], {})
        dm = number(counterpart.get("dm_fitb"))
        fitted = math.isfinite(dm)
        if not fitted:
            dm = float(a["dm"])
        for b in early:
            if int(b["event_id"]) == int(a["event_id"]) or b["tns_name"] == a["tns_name"]:
                continue
            ddm = abs(dm - float(b["dm_fitb"]))
            if ddm > 5 or abs(float(a["dec"]) - float(b["dec"])) > 0.5:
                continue
            name = b["tns_name"] + "_localization.h5"
            with h5py.File(io.BytesIO(read_input(MAPS + name, name, cache, hashes))) as h:
                if h.attrs["tns_name"] != b["tns_name"]:
                    raise ValueError("Localization carries a different event identifier")
                contours = {cl: {island: coords[:].T.tolist() for island, coords in group.items()}
                            for cl, group in h["contours"].items()}
            pairs.append(dict(anchor=a["tns_name"], earlier=b["tns_name"],
                              anchor_id=a["event_id"], earlier_id=b["event_id"],
                              anchor_ra=float(a["ra"]), anchor_dec=float(a["dec"]),
                              anchor_dm=dm, anchor_catalog_dm=float(a["dm"]),
                              dm_source="Cat2 fit" if fitted else "KKO catalog",
                              earlier_dm=float(b["dm_fitb"]), ddm=ddm,
                              earlier_ra=float(b["ra"]), earlier_dec=float(b["dec"]),
                              contours=contours))
    return dict(n_rows=len(cat), n_events=len(events), n_early=len(early),
                n_kko=len(kko), n_anchors=len(anchors), pairs=pairs)


def matches(pair, shifts, level):
    # Unwrap each island at the localization centre to handle the RA seam.
    origin = pair["earlier_ra"]
    ra = (pair["anchor_ra"] + shifts - origin + 180) % 360 - 180
    dec = np.full(ra.shape, pair["anchor_dec"])
    points = np.column_stack([ra, dec])
    hit = np.zeros(len(points), dtype=bool)
    for vertices in pair["contours"][str(level)].values():
        xy = np.array(vertices)
        xy[:, 0] = (xy[:, 0] - origin + 180) % 360 - 180
        in_box = ((ra >= xy[:, 0].min()) & (ra <= xy[:, 0].max()) &
                  (dec >= xy[:, 1].min()) & (dec <= xy[:, 1].max()))
        hit[in_box] |= Polygon(xy).contains_points(points[in_box])
    return hit


def plot_candidates(inputs, filename):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 5.2), constrained_layout=True)
    by_name = {p["earlier"]: p for p in inputs["pairs"]}
    colors = {"95": "#bbcad4", "68": "#146585", "50": "#13465a"}
    for ax, name in zip(axes, ("FRB20210408G", "FRB20220716D")):
        p = by_name[name]
        cos_dec = math.cos(math.radians(p["anchor_dec"]))
        def offset(ra, dec):
            return ((ra - p["anchor_ra"] + 180) % 360 - 180) * cos_dec, dec - p["anchor_dec"]
        for level, color in colors.items():
            for index, coords in enumerate(p["contours"][level].values()):
                ra, dec = np.asarray(coords).T
                ax.plot(*offset(ra, dec), color=color, lw=1.25,
                        label=f"{level}% region" if index == 0 else None)
        ax.scatter(*offset(p["earlier_ra"], p["earlier_dec"]), marker="+", s=95,
                   color="#53565a", zorder=4, label="Catalog 2 position")
        ax.scatter(0, 0, marker="*", s=140, color="#bd401a", zorder=5, label="KKO position")
        ax.set_title(f'{name} with KKO burst {p["anchor"]}', loc="left", fontsize=11)
        ax.set_ylabel("Dec. offset (deg)")
        ax.grid(alpha=0.15)
        ax.invert_xaxis()
        ax.set_aspect("equal", adjustable="datalim")
    axes[-1].set_xlabel("Right ascension offset × cos(dec.) (deg)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center",
               ncols=3, frameon=False, fontsize=9)
    fig.savefig(filename, dpi=155)
    plt.close(fig)


def evaluate(inputs):
    # Rigid rotations keep the KKO sky pattern and all declinations, DMs and maps.
    # Delta RA = 0 is the real sky and is excluded from control counts.
    shifts = np.arange(36000) / 100
    rows = []
    pairs = inputs["pairs"]
    tests = []
    for pair in pairs:
        levels = [cl for cl in (50, 68, 90, 95) if matches(pair, shifts[:1], cl)[0]]
        rows.append({k: v for k, v in pair.items() if k != "contours"} |
                    {"smallest_contour_pct": min(levels) if levels else None})
    for limit in [0.05, 0.5, 2, 5]:
        for level in (50, 68, 95):
            hits = np.zeros(len(shifts), dtype=int)
            for p in pairs:
                if p["ddm"] <= limit:
                    hits += matches(p, shifts, level)
            counts = hits[1:]
            tests.append(dict(max_delta_dm=limit, contour_pct=level,
                              observed=int(hits[0]), n_control=len(counts),
                              mean_control=float(np.mean(counts)),
                              n_ge_one=int(np.count_nonzero(counts >= 1)),
                              n_ge_two=int(np.count_nonzero(counts >= 2))))
    # Coarse local control allows for RA-dependent sky rates and Galactic foreground.
    # The retained shifts displace ordinary fields beyond the typical beam islands.
    near = ((shifts >= 6) & (shifts <= 30)) | ((shifts >= 330) & (shifts <= 354))
    local = np.sum([matches(p, shifts, 68) for p in pairs if p["ddm"] <= 2], axis=0)[near]
    result = {k:v for k,v in inputs.items() if k != "pairs"}
    result.update(n_dm_dec_pairs=len(rows), rotations=tests,
                  local_6_to_30_deg=dict(n_control=len(local), mean_control=float(local.mean()),
                                         n_ge_one=int(sum(local >= 1)), n_ge_two=int(sum(local >= 2))),
                  candidate_rows=[r for r in rows if r["smallest_contour_pct"] is not None],
                  screened_rows=rows)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-archives", type=Path, metavar="CACHE",
                        help="Re-extract each archived contour and verify the compact input")
    parser.add_argument("--plot", type=Path, metavar="FIGURE", help="Draw both localization maps")
    args = parser.parse_args()
    frozen = (ROOT / "contours.json.gz").read_bytes()
    hashes = json.loads((ROOT / "sha256.json").read_text())
    if digest(frozen) != hashes["contours.json.gz"]:
        raise ValueError("Frozen contours checksum mismatch")
    inputs = json.loads(gzip.decompress(frozen))
    if args.from_archives:
        args.from_archives.mkdir(exist_ok=True, parents=True)
        downloaded = build_inputs(args.from_archives, hashes)
        if downloaded != inputs:
            raise ValueError("Archive replay differs from retained input")
    if args.plot:
        plot_candidates(inputs, args.plot)
    print(json.dumps(evaluate(inputs), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
