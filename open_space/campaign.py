"""Advance a durable, bounded SDSS galaxy research campaign on one plate visit."""

import argparse
from contextlib import closing, contextmanager
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys
import tempfile

from . import cohort, ledger
from .dataio import receipt_paths
from .archive import PILOT, write_immutable
from .review import config_from_env
from .recurrence import summarize
from .runner import run_work


@contextmanager
def coordinator(root):
    """An OS lock is released by process exit even after a crash."""
    try:
        import fcntl
    except ImportError as exc:
        raise RuntimeError("This coordinator requires a POSIX filesystem supporting flock") from exc
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".coordinator.lock").open("a") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another coordinator is active in this campaign directory") from exc
        yield


def checkpoint(path, data):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as staged:
        name = staged.name
        try:
            staged.write(cohort.encode(data))
            staged.flush()
            os.fsync(staged.fileno())
        except BaseException:
            Path(name).unlink(missing_ok=True)
            raise
    try:
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def prepare_page(root, cursor, spec):
    page_dir = root / "pages" / f"after-{cursor:04d}"
    page_dir.mkdir(parents=True, exist_ok=True)
    path = page_dir / "plan.json"
    if not path.exists():
        cohort.plan(spec["plate"], spec["mjd"], cursor, spec["limit"], None, path)
    plan_bytes = path.read_bytes()
    plan = json.loads(plan_bytes)
    for key, expected in (("plate", spec["plate"]), ("mjd", spec["mjd"]),
                          ("after_fiber", cursor), ("limit", spec["limit"])):
        if plan[key] != expected:
            raise ValueError(f"Saved page violates campaign setting {key}")
    digest = hashlib.sha256(plan_bytes).hexdigest()
    objects = page_dir / "observations"
    receipt = objects / f"cohort-{digest[:16]}.json"
    if not receipt.exists():
        receipt = cohort.fetch(path, objects)
    if json.loads(receipt.read_bytes())["plan_sha256"] != digest:
        raise ValueError("Fetched receipt differs from saved plan")
    paths = receipt_paths(receipt)
    expected_ids = [cohort.archive_source(t["plate"], t["mjd"], t["fiber"])[0]
                    for t in plan["targets"]]
    expected_paths = [f"sdss-{key}.json" for key in expected_ids]
    if [p.name for p in paths] != expected_paths:
        raise ValueError("Fetched receipt must include every planned target in order")
    return page_dir, paths, plan["targets"]


def run_campaign(root, plate, mjd, limit, pages, max_calls, config, code):
    cohort.selection_query(plate, mjd, 0, limit)
    if pages < 1 or max_calls < 0:
        raise ValueError("--pages must be positive and --max-calls nonnegative")
    # The visit has at most 1000 fibers and each page has at least one target
    # except for the terminal query. No choice of --pages can grow the survey.
    spec = {
        "schema_version": 1, "plate": plate, "mjd": mjd, "limit": limit,
        "software_sha256": code, "model_config": config,
    }
    identity = hashlib.sha256(cohort.encode(spec)).hexdigest()
    state_file = root / "state.json"
    calls = jobs_done = pages_done = 0
    with coordinator(root):
        config_file = root / "campaign.json"
        if config_file.exists() and json.loads(config_file.read_bytes()) != spec:
            raise ValueError("Saved software, selection or model settings differ. Use the original settings or a new campaign directory")
        write_immutable(config_file, cohort.encode(spec))
        if state_file.exists():
            state = json.loads(state_file.read_bytes())
            if state["campaign_sha256"] != identity:
                raise ValueError("Cursor belongs to a different campaign")
        else:
            state = {"campaign_sha256": identity, "after_fiber": 0, "exhausted": False}
            checkpoint(state_file, state)
        status = "exhausted" if state["exhausted"] else "page_limit"
        current_job = None
        with closing(ledger.connect(root / "ledger.sqlite")) as db:
            while pages_done < pages and not state["exhausted"]:
                page, paths, targets = prepare_page(root, state["after_fiber"], spec)
                reports = page / "reports"
                reports.mkdir(exist_ok=True)
                for manifest in paths:
                    job = ledger.enqueue(db, manifest, config, code)
                    existing = ledger.evidence(db, job)
                    if existing["state"] == "failed":
                        # One retry per job per invocation, preserving completed stages.
                        ledger.retry(db, job, code)
                    if existing["state"] != "completed":
                        stats = run_work(db, code, max_jobs=1, max_calls=max_calls-calls, only_job=job)
                        calls += stats["provider_calls_reserved"]
                        jobs_done += stats["completed"]
                        if stats["paused_for_call_budget"]:
                            status = "paused_for_call_budget"
                    packet = ledger.evidence(db, job)
                    checkpoint(reports / f"{manifest.stem}.json", packet)
                    if packet["state"] != "completed":
                        current_job = job
                        if status != "paused_for_call_budget":
                            status = "failed" if packet["state"] == "failed" else "waiting_for_worker"
                        break
                if current_job is not None:
                    break
                # Publish the cross-target screen before committing the page cursor.
                # On replay the same per-target reports regenerate the same result.
                packets = (json.loads(p.read_bytes()) for p in
                           sorted((root / "pages").glob("*/reports/*.json")))
                checkpoint(root / "recurrence.json", summarize(packets))
                if targets:
                    state["after_fiber"] = targets[-1]["fiber"]
                state["exhausted"] = len(targets) < limit or state["after_fiber"] == 1000
                checkpoint(state_file, state)
                pages_done += 1
                status = "exhausted" if state["exhausted"] else "page_limit"
    return {
        **state, "status": status, "pages_committed": pages_done,
        "jobs_completed_this_invocation": jobs_done,
        "provider_calls_reserved_this_invocation": calls,
        "current_job": current_job, "root": str(root),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("runs/campaign"))
    parser.add_argument("--plate", type=int, default=PILOT[0])
    parser.add_argument("--mjd", type=int, default=PILOT[1])
    parser.add_argument("--limit", type=int, default=3, help="Maximum spectra per catalog page")
    parser.add_argument("--pages", type=int, default=1, help="Maximum completed pages per invocation")
    parser.add_argument("--models", action="store_true", help="Enable model-led follow-up")
    parser.add_argument("--max-calls", type=int, default=6, help="Total model-call reservations per invocation")
    args = parser.parse_args(argv)
    try:
        result = run_campaign(args.root, args.plate, args.mjd, args.limit, args.pages, args.max_calls,
                              config_from_env() if args.models else None, ledger.fingerprint())
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, sqlite3.Error) as exc:
        print(f"open-space campaign: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return int(result["status"] in ("failed", "waiting_for_worker"))


if __name__ == "__main__":
    raise SystemExit(main())
