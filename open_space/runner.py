"""Run deterministic detection, model-selected measurements and adversarial review."""

import json

from . import ledger
from .followup import measure
from .review import run_role
from .spectra import analyze


def run_work(db, code, max_jobs, max_calls, role=run_role, only_job=None):
    if max_jobs < 1 or max_calls < 0:
        raise ValueError("max-jobs must be positive and max-calls must be nonnegative")
    stats = {"completed": 0, "failed": 0, "paused_for_call_budget": 0,
             "lost_leases": 0, "provider_calls_reserved": 0}
    for _ in range(max_jobs):
        job = ledger.claim(db, code, only_job)
        if job is None:
            break
        try:
            manifest, rows = ledger.source(db, job)
            config = json.loads(job["config"])
            saved = ledger.read_stages(db, job["id"])
            if "detection" not in saved:
                scan = analyze(rows, manifest["redshift"])
                scan["source_context"] = {"kind": "observed", "observation_id": manifest["observation_id"]}
                ledger.checkpoint(db, job, "detection", scan)
                saved["detection"] = scan
            scan = saved["detection"]
            if config and scan["features"]:
                for stage in ("analyst", "measurements", "reviewer"):
                    if stage in saved:
                        continue
                    if stage == "measurements":
                        result = measure(saved["analyst"]["requests"], scan["features"], rows)
                    else:
                        if stats["provider_calls_reserved"] >= max_calls:
                            ledger.finish(db, job, "queued")
                            stats["paused_for_call_budget"] += 1
                            return stats
                        ledger.reserve_call(db, job)
                        stats["provider_calls_reserved"] += 1
                        previous = None if stage == "analyst" else {
                            "analyst": saved["analyst"], "followup_measurements": saved["measurements"]}
                        result = role(stage, scan, previous, config)
                    ledger.checkpoint(db, job, stage, result)
                    saved[stage] = result
            ledger.finish(db, job, "completed")
            stats["completed"] += 1
        except ledger.LostLease:
            stats["lost_leases"] += 1
        except Exception as exc:
            # Keep successful, billable work. Explicit retry resumes the missing role.
            try:
                ledger.finish(db, job, "failed", str(exc)[:700])
            except ledger.LostLease:
                stats["lost_leases"] += 1
            else:
                stats["failed"] += 1
    return stats
