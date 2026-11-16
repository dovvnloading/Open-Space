"""Verify the repository's own claims from its committed bytes.

Examples:

  python -m open_space.verify                       check the last commit
  python -m open_space.verify --worktree            check uncommitted work as it would be committed
  python -m open_space.verify --result NAME         check one result folder
  python -m open_space.verify --replay              also run every documented command offline
  python -m open_space.verify --report R.md --json L.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from . import report, runner
from .tree import RepoTree, TreeError


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m open_space.verify", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("root", nargs="?", default=".", type=Path, help="repository to verify (default: current directory)")
    source = p.add_mutually_exclusive_group()
    source.add_argument("--rev", help="Git revision to verify (default: HEAD)")
    source.add_argument("--worktree", action="store_true",
                        help="verify the working tree, staged through a throwaway index, instead of a commit")
    p.add_argument("--checks", default=",".join(runner.DEFAULT_CHECKS),
                   help=f"comma-separated checks to run, from: {', '.join(runner.ALL_CHECKS)}")
    p.add_argument("--replay", action="store_true", help="add the replay check (runs documented commands; slow)")
    p.add_argument("--result", action="append", metavar="NAME", help="limit to one result folder; may be repeated")
    p.add_argument("--json", type=Path, metavar="FILE", help="write the machine-readable ledger")
    p.add_argument("--full", action="store_true", help="include every match and every external receipt in the ledger")
    p.add_argument("--report", type=Path, metavar="FILE", help="write the Markdown report")
    p.add_argument("--strict", action="store_true", help="exit with status 1 on warnings as well as failures")
    p.add_argument("--quiet", action="store_true", help="print nothing but errors")
    p.add_argument("--python", metavar="EXE", help="interpreter for replayed commands (default: this one)")
    p.add_argument("--replay-timeout", type=float, default=120.0, metavar="SECONDS")
    p.add_argument("--replay-jobs", type=int, default=2, metavar="N")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")     # reports contain characters outside cp1252
    checks = [name.strip() for name in args.checks.split(",") if name.strip()]
    if args.replay and "replay" not in checks:
        checks.append("replay")
    try:
        tree = RepoTree.open(args.root, rev=args.rev, worktree=args.worktree)
    except TreeError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    only = set(args.result) if args.result else None
    if only:
        unknown = sorted(only - set(tree.result_dirs()))
        if unknown:
            print(f"error: no such result folder: {', '.join(unknown)}", file=sys.stderr)
            return 2

    def progress(message: str) -> None:
        if not args.quiet:
            print(f"… {message}", file=sys.stderr, flush=True)

    try:
        ledger = runner.run_checks(
            tree, checks, only, progress=progress,
            replay_options={"python": args.python, "timeout": args.replay_timeout, "jobs": args.replay_jobs})
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(runner.to_json(ledger, args.full), indent=1, ensure_ascii=False,
                                        sort_keys=False) + "\n", encoding="utf-8", newline="\n")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(report.markdown(ledger), encoding="utf-8", newline="\n")
    if not args.quiet:
        print(report.console(ledger))
    return runner.exit_code(ledger, args.strict)


if __name__ == "__main__":
    raise SystemExit(main())
