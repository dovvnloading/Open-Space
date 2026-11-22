"""Search a prespecified moderate-count cohort first, preserving the full parent survey."""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import tempfile

import campaign
import report
from archive import dump

GROUPS = 8
OFFSET = 192
PER_STRATUM = 32
UPPER = (8000, 24000, 72000, 200000)
COHORT = "moderate-count-v1"


def select():
    evidence = campaign.EVIDENCE
    inventory = evidence / "inventory/inventory.json"
    freeze = campaign.read(evidence / "inventory/campaign.json")
    if hashlib.sha256(inventory.read_bytes()).hexdigest() != freeze["inventory_sha256"]:
        raise ValueError("Inventory no longer matches the frozen sample")
    rows = campaign.read(inventory)
    groups = [[] for _ in range(GROUPS)]
    strata = []
    lower = -1
    for upper in UPPER:
        eligible = [r for r in rows if r.get("event_count") is not None
                    and lower < r["event_count"] <= upper]
        # The order uses no photons, sky positions, source classes or trigger outcomes.
        eligible.sort(key=lambda r: hashlib.sha256(
            f"chandra-acis-fxt:{COHORT}:{r['obsid']}".encode("ascii")).digest())
        selected = [r["obsid"] for r in eligible[:PER_STRATUM]]
        for i, oid in enumerate(selected):
            groups[i % GROUPS].append(oid)
        strata.append({"count_above": lower, "count_at_most": upper, "available": len(eligible),
                       "selected_in_hash_order": selected})
        lower = upper
    return {"cohort": COHORT, "inventory_sha256": freeze["inventory_sha256"],
            "rank": "ascending SHA256 of ASCII chandra-acis-fxt:moderate-count-v1:<decimal ObsID>",
            "strata": strata, "groups": groups}


def run(group, output):
    selection = select()
    if group not in range(GROUPS):
        raise ValueError("Early group outside the fixed selection")
    evidence = campaign.EVIDENCE
    ids = selection["groups"][group]
    pending = []
    for oid in ids:
        path = evidence / "observations" / f"{oid}.json.gz"
        if not path.exists() or campaign.read_obs(path)["status"] != "complete":
            pending.append(oid)
    # Only this ephemeral worker sees the alternate assignments. The offset
    # also gives its cached responses distinct paths when the full run resumes.
    assignments = [[] for _ in range(OFFSET + GROUPS)]
    assignments[OFFSET + group] = pending
    dump(evidence / "inventory/shards.json", assignments)
    campaign.shard(OFFSET + group, output)
    dump(output / f"moderate-group-{group}.json", {"selected": ids, "assigned": pending,
                                                  "git_commit": os.environ["GITHUB_SHA"]})


def deposit():
    evidence = campaign.EVIDENCE
    selection = select()
    ids = [oid for group in selection["groups"] for oid in group]
    rel = Path("provisional") / COHORT / (os.environ["GITHUB_RUN_ID"] + "-" + os.environ["GITHUB_RUN_ATTEMPT"])
    dest = evidence / rel
    dest.mkdir(parents=True, exist_ok=True)
    dump(dest / "selection.json", selection)
    with tempfile.TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as tmp:
        view = Path(tmp)
        (view / "inventory").symlink_to(evidence / "inventory", target_is_directory=True)
        (view / "observations").mkdir()
        for oid in ids:
            file = evidence / "observations" / f"{oid}.json.gz"
            if file.is_file():
                (view / "observations" / file.name).symlink_to(file)
        report.EVIDENCE = view
        try:
            report.run()
        finally:
            report.EVIDENCE = evidence
        for item in view.iterdir():
            if item.is_file():
                shutil.copyfile(item, dest / item.name)
    page = campaign.HERE / "README.md"
    text = page.read_text().replace("(evidence/", f"(evidence/{rel}/")
    headline = text.splitlines()[2].strip("*")
    scope = (f"The early cohort selected {len(ids)} distinct observations in four catalog event-count "
             "strata before opening their event lists. It excludes observations above 200,000 catalog "
             "counts and observations without a count estimate. The denominator includes only completed "
             "product scans, and all pending observations remain assigned in the full survey.")
    checks = (f"[Selected ObsIDs and selection rule](evidence/{rel}/selection.json) · "
              "[Per-observation ledgers](evidence/observations) · [Field notes](FIELD_NOTES.md)")
    page.write_text(text + "\n" + scope + "\n\n" + checks + "\n")
    dump(dest / "scope.json", {"selected": len(ids), "scope": scope})
    summary = campaign.read(dest / "summary.json")
    if summary["label"]:
        report.update_catalog("Early count-stratified cohort. " + headline, summary["label"])
    else:
        raise RuntimeError("The early sample yielded no supported catalog measurement")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(required=True, dest="command")
    worker = sub.add_parser("group")
    worker.add_argument("index", type=int)
    worker.add_argument("out", type=Path)
    sub.add_parser("deposit")
    args = p.parse_args()
    if args.command == "group":
        run(args.index, args.out)
    else:
        deposit()
