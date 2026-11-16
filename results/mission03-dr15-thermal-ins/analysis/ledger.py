"""Join the unique source queues to independent cone results and provisional dispositions."""

import argparse
import csv
import json
from pathlib import Path
from summarize import summarize


def csvrows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def run(results, overlap, cache):
    catalog = {}
    queues = {}
    prior = {}
    for queue in ["measured", "boundary", "rescue", "bright", "difference"]:
        for r in csvrows(results / (queue + ".csv")):
            sid = r["SRCID"]
            catalog[sid] = r
            queues.setdefault(sid, []).append(queue)
        for r in csvrows(overlap / (queue + "-prior.csv")):
            if r["within_review_radius"] == "1":
                prior.setdefault(r["SRCID"], set()).add(r["prior"])
    sids = {
        v["SRCID"]: v
        for fn in ["simbad.json", "difference_simbad.json"]
        for v in json.loads((results / fn).read_text())["objects"]
    }
    erass = {
        v["SRCID"]: v
        for fn in ["erass1.json", "difference_erass1.json"]
        for v in json.loads((results / fn).read_text())["objects"]
    }
    # Every deduplicated target gets its own row regardless of follow-up outcome.
    allq = cache / "all-targets.csv"
    write(
        allq,
        [
            {
                k: r.get(k, "")
                for k in [
                    "SRCID",
                    "IAUNAME",
                    "RA",
                    "DEC",
                    "F_SOFT",
                    "RADEC_ERR",
                    "MJD_FIRST",
                    "MJD_LAST",
                ]
            }
            for r in catalog.values()
        ],
    )
    phot = list(summarize(allq, cache))
    write(results / "counterparts.csv", phot)
    ledger = []
    for p in phot:
        sid = p["srcid"]
        sim = sids[sid]["matches"]
        e = erass[sid]
        names = sorted(prior.get(sid, []))
        close = []
        for c, sep, mag in [
            ("gaia", "gaia_sep", "gaia_mag"),
            ("legacy", "legacy_sep", "legacy_mag"),
            ("ps1", "ps1_r_snr5_sep", "ps1_r_snr5_mag"),
            ("allwise", "allwise_sep", "allwise_mag"),
        ]:
            if p[sep] != "" and p[sep] <= p["r99_arcsec"] and p[mag] != "":
                close.append(c)
        if names:
            reason = "prior-position review"
        elif close:
            reason = "possible optical/IR association"
        elif sim:
            reason = "SIMBAD neighbor requires review"
        else:
            reason = "no positional veto; imaging/coverage required"
        if sid in ("3086526020100024", "3078029060100001"):
            reason = "rejected: original events concentrate in detector pixels"
        if sid == "3030258180100013":
            reason = "original EPIC and deep CFHT reviewed; ambiguous ultrasoft source"
        if sid == "3010946500100005":
            reason = (
                "original pn reviewed; optical limit weak and long-term flux unresolved"
            )
        if sid == "3012772050100005":
            reason = "nova wings and camera edge compromise detection"
        if sid == "3093239080100009":
            reason = (
                "nearby red G=15.16 star at 3 arcsec; counterpart remains plausible"
            )
        if sid == "3084132010100270":
            reason = "known [HP99] 159 white-dwarf binary field"
        if reason.startswith("no positional veto"):
            obsid = sid[1:11]
            if obsid in {
                "0111090301",
                "0077340101",
                "0077340201",
                "0503220101",
                "0112870401",
                "0113050601",
            }:
                reason = (
                    "supernova-remnant field: unresolved diffuse-emission contamination"
                )
        ledger.append(
            dict(
                srcid=sid,
                iauname=p["name"],
                queues="|".join(queues[sid]),
                prior_names="|".join(names),
                possible_counterpart_surveys="|".join(close),
                simbad_neighbors_15arcsec=len(sim),
                simbad_nearest=sim[0]["main_id"] if sim else "",
                erass1_in_public_hemisphere=e["in_western_hemisphere"],
                erass1_neighbors_30arcsec=len(e["matches"]),
                disposition=reason,
            )
        )
    write(results / "review-ledger.csv", ledger)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("results", "overlap", "cache"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    run(a.results, a.overlap, a.cache)
