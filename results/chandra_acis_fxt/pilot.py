"""Checkpoint three observations from the frozen inventory using the production search."""

import hashlib
import os
import re
import tempfile
from pathlib import Path

from campaign import (EVIDENCE, N_CONTROL, N_INJECT, dump, identity, process,
                      read, read_obs, save_obs)

COUNT_TARGETS = (10_000, 30_000, 60_000)


def select_entries(rows):
    eligible = [row for row in rows if (row["exposure_s"] or 0) > 5000 and (row["event_count"] or 0) > 0]
    selected = []
    for target in COUNT_TARGETS:
        if not eligible:
            break
        entry = min(eligible, key=lambda row: (abs(row["event_count"] - target), int(row["obsid"])))
        selected.append({"obsid": entry["obsid"], "count_target": target,
                         "catalog_event_count": entry["event_count"], "catalog_exposure_s": entry["exposure_s"]})
        eligible.remove(entry)
    return selected


def clean_error(error):
    # Emit only a short single line; strip any signed query in a diagnostic URL.
    msg = re.sub(r"(https?://[^\s?#]+)[?#][^\s]+", r"\1", str(error))
    msg = re.sub(r"[^ -~]+", " ", msg)
    return " ".join(msg.split())[:220]


def scored_products(ledger):
    scored = []
    for product in ledger.get("products", []):
        scan = product.get("search", {})
        if (scan.get("exposure_s", 0) > 0 and len(scan.get("controls", [])) == N_CONTROL
                and len(scan.get("injections", [])) == N_INJECT):
            scored.append({"filename": product["input"]["filename"], "sha256": product["input"]["sha256"],
                           "exposure_60s_start_support_s": scan["exposure_s"],
                           "real_triggers": len(scan["triggers"]), "control_trials": N_CONTROL,
                           "injection_trials": N_INJECT, "followup_processed": bool(product.get("complete"))})
    return scored


def main():
    ref = "refs/heads/research/chandra-acis-fast-transients"
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != ref
            or os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("The early checkpoint requires the research branch's Actions campaign")
    run_id = os.environ["GITHUB_RUN_ID"]
    attempt = os.environ["GITHUB_RUN_ATTEMPT"]
    temp_root = Path(os.environ["RUNNER_TEMP"])
    os.environ.setdefault("RETAIN_EVENTS", str(temp_root / "candidate-inputs"))
    os.environ.setdefault("CIAO_PREFIX", str(temp_root / "ciao"))
    inv_path = EVIDENCE / "inventory" / "inventory.json"
    freeze = read(inv_path.with_name("campaign.json"))
    if (freeze["implementation"] != identity()
            or freeze["inventory_sha256"] != hashlib.sha256(inv_path.read_bytes()).hexdigest()
            or freeze["n_controls"] != N_CONTROL or freeze["n_injections_per_product"] != N_INJECT):
        raise ValueError("Checkpoint inputs or search settings differ from the frozen campaign")

    selected = select_entries(read(inv_path))
    checkpoint = {"inventory_sha256": freeze["inventory_sha256"], "run_id": run_id, "run_attempt": attempt,
                  "git_commit": os.environ["GITHUB_SHA"], "count_targets": list(COUNT_TARGETS),
                  "selection": "Distinct positive-count entries with catalog exposure >5000 s; "
                               "nearest count to each target in order, then lowest ObsID.",
                  "selected": selected, "observations": []}
    checkpoint_path = EVIDENCE / "pilot.json"
    dump(checkpoint_path, checkpoint)
    cache = EVIDENCE / "catalogs" / f"{run_id}-{attempt}-pilot"
    cache.mkdir(parents=True, exist_ok=True)
    for entry in selected:
        obsid = entry["obsid"]
        path = EVIDENCE / "observations" / f"{obsid}.json.gz"
        prior = read_obs(path) if path.is_file() else {}
        reused = prior.get("status") == "complete"
        if not reused:
            print(f"Checkpoint ObsID {obsid} starts at {entry['catalog_event_count']} catalog events", flush=True)
            try:
                with tempfile.TemporaryDirectory(dir=temp_root) as tmp:
                    process(obsid, EVIDENCE, Path(tmp), cache)
            except Exception as exc:
                ledger = read_obs(path) if path.is_file() else {"obsid": obsid, "products": [], "errors": []}
                ledger["status"] = "error"
                ledger["errors"].append(clean_error(f"{type(exc).__name__}: {exc}"))
                path.parent.mkdir(parents=True, exist_ok=True)
                save_obs(path, ledger)
        ledger = read_obs(path)
        scored = scored_products(ledger)
        errors = [clean_error(e) for e in ledger.get("errors", [])]
        checkpoint["observations"].append({"obsid": obsid, "status": ledger["status"], "reused": reused,
                                           "scored_products": scored, "errors": errors})
        dump(checkpoint_path, checkpoint)
        print(f"Checkpoint ObsID {obsid}: {ledger['status']}, {len(scored)} scored products, "
              f"{sum(p['real_triggers'] for p in scored)} real triggers", flush=True)
        for error in errors[:2]:
            print(f"Checkpoint ObsID {obsid} diagnostic: {error}", flush=True)
    if not any(row["scored_products"] for row in checkpoint["observations"]):
        raise RuntimeError("No selected observation yielded a scored product with positive search exposure")


if __name__ == "__main__":
    main()
