"""Advance an unattended TESS survey through deterministic, nonoverlapping sky samples."""

import argparse
import json
from pathlib import Path
import sys

from . import discovery, discovery_agents
from .archive import write_immutable
from .planet_catalog import snapshot
from .tess_archive import (
    SEED,
    CONTROL_TIC,
    inventory,
    atomic_bytes,
    fetch_plan,
    make_plan,
)


def encoded(value):
    return (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()


def run(
    root,
    search,
    validation,
    start,
    size,
    pages,
    workers,
    jobs,
    model_calls,
    models=False,
    seed=SEED,
    cache_root=None,
):
    """One survey driver per directory. Each page gets its own recoverable worker queue."""
    import fcntl

    if (
        min(search, validation, size, pages, workers, jobs) < 1
        or workers > 128
        or search >= validation
        or start < 0
        or model_calls < 0
    ):
        raise ValueError("Invalid survey limits or sector pair")
    root.mkdir(parents=True, exist_ok=True)
    cache = cache_root or root / "cache"
    paths = {k: cache / k for k in ("tess", "wise", "pixels", "catalog")}
    settings = discovery_agents.settings() if models else None
    project = {
        "schema_version": 1,
        "sectors": [search, validation],
        "start": start,
        "page_size": size,
        "seed": seed,
        "model": settings,
        "code_sha256": discovery.fingerprint(),
    }
    with (root / "driver.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(
                "Another survey driver is advancing this directory"
            ) from exc
        inventories = root / "inventories"
        inventories.mkdir(exist_ok=True)
        selections = []
        receipts = []
        for sector in (search, validation):
            frozen = inventories / f"s{sector}.sh"
            if not frozen.exists():
                inventory(sector, cache / "inventories")
                write_immutable(
                    frozen, (cache / "inventories" / f"s{sector}.sh").read_bytes()
                )
            products, receipt = inventory(sector, inventories)
            selections.append(products.keys())
            receipts.append(receipt)
        total = len((selections[0] & selections[1]) - {CONTROL_TIC})
        project.update(inventories=receipts, eligible_targets=total)
        write_immutable(root / "survey.json", encoded(project))
        cursor = root / "progress.json"
        state = (
            discovery.read(cursor)
            if cursor.exists()
            else {
                "next_offset": start,
                "completed_pages": 0,
                "completed_model_calls": 0,
                "explored_stars": 0,
                "exhausted": start >= total,
            }
        )
        count = 0
        while count < pages and not state["exhausted"]:
            offset = state["next_offset"]
            page = root / "pages" / f"{offset:08d}"
            page.mkdir(parents=True, exist_ok=True)
            manifest, catalog = page / "targets.json", page / "planet-catalog.json"
            if not manifest.exists():
                plan = make_plan(search, validation, size, offset, seed, inventories)
                write_immutable(
                    manifest, encoded(fetch_plan(plan, paths["tess"], min(workers, 16)))
                )
            targets = discovery.read(manifest)["targets"]
            if not catalog.exists():
                write_immutable(
                    catalog,
                    encoded(snapshot((t["tic_id"] for t in targets), paths["catalog"])),
                )
            discovery.initialize(page, [manifest], [], catalog, settings)
            ceiling = max(0, model_calls - state["completed_model_calls"])
            measured = discovery.work(page, paths, workers, jobs, 3, ceiling)
            if set(measured["tasks"]) - {"done"}:
                return {
                    **state,
                    "active_page": str(page),
                    "page_counts": measured,
                    "status": (
                        "needs_retry"
                        if measured["tasks"].get("failed")
                        else "in_progress"
                    ),
                }
            discovery.export(page, page / "export")
            state = {
                "next_offset": offset + len(targets),
                "explored_stars": state["explored_stars"] + len(targets),
                "completed_pages": state["completed_pages"] + 1,
                "exhausted": offset + len(targets) >= total,
                "completed_model_calls": state["completed_model_calls"]
                + measured["model_calls_reserved"],
            }
            atomic_bytes(cursor, encoded(state))
            count += 1
    return {
        **state,
        "status": "exhausted" if state["exhausted"] else "page_budget_reached",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--sectors", nargs=2, type=int, default=[1, 2])
    parser.add_argument(
        "--start", type=int, default=0, help="Initial rank in the deterministic sample"
    )
    parser.add_argument("--page-size", type=int, default=32)
    parser.add_argument(
        "--pages",
        type=int,
        default=1,
        help="Maximum completed new pages per invocation",
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--jobs-per-page", type=int, default=500)
    parser.add_argument(
        "--model-call-ceiling",
        type=int,
        default=0,
        help="Cumulative ceiling across all pages",
    )
    parser.add_argument("--models", action="store_true")
    parser.add_argument("--seed", default=SEED)
    parser.add_argument("--cache-root", type=Path)
    a = parser.parse_args()
    try:
        result = run(
            a.root,
            *a.sectors,
            a.start,
            a.page_size,
            a.pages,
            a.workers,
            a.jobs_per_page,
            a.model_call_ceiling,
            a.models,
            a.seed,
            a.cache_root,
        )
    except (OSError, ValueError) as exc:
        parser.exit(1, f"survey: {exc}\n")
    print(json.dumps(result, indent=2))
    return int(result["status"] == "needs_retry")


if __name__ == "__main__":
    sys.exit(main())
