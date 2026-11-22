"""Measure the response to published timing events after the blind search."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import t as student

from search import fit_pulsar, measurable, score, waveform

EVENTS = (("J1713+0747", 54750.), ("J1713+0747", 57510.), ("J1643-1224", 57074.))
TEMPLATES = ("dm_tail", "profile_like_tail", "achromatic_pulse")
WIDTH = 60.


def run(data, evidence):
    rows = []
    for pulsar, epoch in EVENTS:
        primary = json.loads((evidence / "single" / (pulsar + ".json")).read_text())
        for dataset in ("nb", "wb", "epta"):
            path = data / dataset / (pulsar + ".npz")
            projection = fit_pulsar(path, fit_soft=False) if path.is_file() else None
            for name in TEMPLATES:
                item = dict(pulsar=pulsar, published_epoch_mjd=epoch, dataset=dataset, template=name, width_d=WIDTH)
                rows.append(item)
                if projection is None:
                    item["status"] = "no accepted timing reduction"
                    continue
                p = projection
                delta = p.t - epoch
                before = np.unique(np.floor(p.t[(delta < 0) & (delta >= -WIDTH)]))
                after = np.unique(np.floor(p.t[(delta >= 0) & (delta <= WIDTH)]))
                item.update(days_before=len(before), days_after=len(after))
                if not len(before) or not len(after):
                    item["status"] = "no bracketing observations within one width"
                    continue
                raw = waveform(p, name, epoch, WIDTH)
                norms, _, tval, amp, se = score(p, p.whiten(raw))
                if not measurable(p, raw, norms)[0]:
                    item["status"] = "unidentifiable"
                    continue
                item.update(status="tested", amplitude_us=float(amp[0]), se_us=float(se[0]),
                            t=float(tval[0]), fixed_two_sided_p=float(2*student.sf(abs(tval[0]), p.dof-1)))
                if dataset == "nb":
                    same = [r for r in primary.get("primary_crossings", []) if r["template"] == name and
                            r["width_d"] == WIDTH and abs(r["epoch_mjd"]-epoch) <= 30]
                    item["blind_crossings_within_30d"] = len(same)
                    if same:
                        nearest = min(same, key=lambda r: abs(r["epoch_mjd"]-epoch))
                        item["nearest_blind_mjd"] = nearest["epoch_mjd"]
                        item["nearest_blind_global_p"] = nearest["p_global"]
                        item["nearest_blind_failed_checks"] = ",".join(nearest["unpassed"])
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (evidence / "known-event-controls.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--evidence", required=True, type=Path)
    args = parser.parse_args()
    run(args.data, args.evidence)
