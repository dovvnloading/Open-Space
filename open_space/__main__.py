"""Command-line entry point for campaigns, workers and evidence exports."""

import argparse
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import sys

from . import ledger
from .dataio import receipt_paths
from .review import config_from_env
from .runner import run_work
from .spectra import benchmark


PILOT = Path("data/sdss-6497-56329-0139.json")


def emit(value):
    print(json.dumps(value, indent=2, allow_nan=False))


def parser():
    p = argparse.ArgumentParser(prog="open-space", description=__doc__)
    p.add_argument("--db", type=Path, default=Path("runs/ledger.sqlite"))
    cmds = p.add_subparsers(dest="command", required=True)
    cmds.add_parser("pilot", help="Verify and analyze the bundled real galaxy spectrum offline")
    add = cmds.add_parser("enqueue", help="Pin observed inputs and create idempotent jobs")
    add.add_argument("manifests", nargs="*", type=Path)
    add.add_argument("--receipt", type=Path, help="Cohort fetch receipt")
    add.add_argument("--models", action="store_true", help="Enable model-led follow-up for these jobs")
    work = cmds.add_parser("work", help="Drain a bounded slice of the queue, then exit")
    work.add_argument("--max-jobs", type=int, default=1)
    work.add_argument("--max-calls", type=int, default=2)
    cmds.add_parser("status", help="List jobs, states, errors and call counts")
    report = cmds.add_parser("report", help="Export JSON evidence, including partial results")
    report.add_argument("job_id")
    retry = cmds.add_parser("retry", help="Requeue a failed job and retain successful stages")
    retry.add_argument("job_id")
    bench = cmds.add_parser("benchmark", help="Evaluate deterministic synthetic injections and nulls")
    bench.add_argument("--seeds", type=int, default=12)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "benchmark":
            emit(benchmark(args.seeds))
            return 0
        code = ledger.fingerprint()
        with closing(ledger.connect(args.db)) as db:
            if args.command == "pilot":
                job = ledger.enqueue(db, PILOT, None, code)
                # Pilot execution is isolated from the general queue and never spends tokens.
                if ledger.evidence(db, job)["state"] != "completed":
                    run_work(db, code, 1, 0, only_job=job)
                packet = ledger.evidence(db, job)
                emit(packet)
                return 0 if packet["state"] == "completed" else 1
            if args.command == "enqueue":
                paths = args.manifests + (receipt_paths(args.receipt) if args.receipt else [])
                if not paths and not args.receipt:
                    raise ValueError("Supply observation manifests or a cohort receipt")
                config = config_from_env() if args.models else None
                emit([ledger.enqueue(db, path, config, code) for path in paths])
            elif args.command == "work":
                stats = run_work(db, code, args.max_jobs, args.max_calls)
                emit(stats)
                return int(stats["failed"] > 0 or stats["lost_leases"] > 0)
            elif args.command == "status":
                emit(ledger.summary(db))
            elif args.command == "retry":
                ledger.retry(db, args.job_id, code)
                emit({"job_id": args.job_id, "state": "queued"})
            elif args.command == "report":
                emit(ledger.evidence(db, args.job_id))
            return 0
    except (OSError, ValueError, TypeError, KeyError, sqlite3.Error) as exc:
        print(f"open-space: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
