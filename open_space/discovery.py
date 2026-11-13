"""Coordinate resumable, multi-mission searches and export their discovery register."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
import json
import platform
from importlib.metadata import version
from pathlib import Path
import sqlite3
import sys
import threading
import urllib.error

from . import discovery_agents, discovery_queue as queue
from .archive import write_immutable
from .dataio import canonical, digest
from .ledger import fingerprint, transaction
from .planet_catalog import assess
from .tess_archive import atomic_bytes

DEFAULT_TESS = Path("data/tess-plan.json")
DEFAULT_WISE = Path("data/wise-fields.json")
DEFAULT_CATALOG = Path("data/discovery-planet-catalog.json")
DB = "discovery.sqlite"
SIGNALS = {"two_sector_periodic_dip", "periodic_dip_with_veto"}


def read(path):
    return json.loads(Path(path).read_bytes())


def runtime_versions():
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "numpy": version("numpy"),
        "astropy": version("astropy"),
    }


def initialize(
    root,
    tess_plans,
    wise_plans,
    catalog_path,
    config=None,
    shards=1,
    shard=0,
    ffi_plans=(),
):
    if not 0 <= shard < shards:
        raise ValueError("Use a shard index from zero through shards minus one")
    catalog = read(catalog_path)
    all_targets = {}
    fields = {}
    for path in tess_plans:
        for t in read(path)["targets"]:
            ident = (
                str(t["tic_id"])
                + "/"
                + ",".join(str(p["sector"]) for p in t["products"])
            )
            if ident in all_targets and t != all_targets[ident]:
                raise ValueError("Conflicting TESS target products in input plans")
            all_targets[ident] = t
    for path in wise_plans:
        for f in read(path)["fields"]:
            if f["id"] in fields and f != fields[f["id"]]:
                raise ValueError("Conflicting WISE field plans")
            fields[f["id"]] = f
    image_fields = {}
    for path in ffi_plans:
        for f in read(path)["fields"]:
            if f != image_fields.setdefault(f["field_id"], f):
                raise ValueError("Conflicting FFI fields")
    if not (fields or all_targets or image_fields):
        raise ValueError("Empty campaign")

    def selected(key):
        return int(digest(key.encode()), 16) % shards == shard

    # Group both optical searches on the same shard for cache locality.
    jobs = [
        (
            kind,
            "TIC " + ident,
            {"target": t, "catalog_targets": catalog["targets"].get(str(t["tic_id"]))},
        )
        for ident, t in sorted(all_targets.items())
        if selected("TIC " + ident)
        for kind in ("transit", "variability")
    ]
    jobs.extend(
        ("motion", ident, {"field": f})
        for ident, f in sorted(fields.items())
        if selected(ident)
    )
    jobs.extend(
        ("ffi", ident, {"field": f})
        for ident, f in sorted(image_fields.items())
        if selected(ident)
    )
    project = {
        "schema_version": 1,
        "code_sha256": fingerprint(),
        "model": config,
        "runtime_versions": runtime_versions(),
        "shards": shards,
        "shard": shard,
        "tess_plan_sha256": [digest(Path(p).read_bytes()) for p in tess_plans],
        "ffi_plan_sha256": [digest(Path(p).read_bytes()) for p in ffi_plans],
        "wise_plan_sha256": [digest(Path(p).read_bytes()) for p in wise_plans],
        "catalog_sha256": digest(Path(catalog_path).read_bytes()),
        "seed_tasks": len(jobs),
    }
    path = root / DB
    queue.initialize(path, project)
    inputs = root / "inputs"
    inputs.mkdir(exist_ok=True)
    for source in [*tess_plans, *wise_plans, *ffi_plans, catalog_path]:
        raw = Path(source).read_bytes()
        write_immutable(inputs / (digest(raw) + ".json"), raw)
    with closing(queue.connect(path)) as db, transaction(db):
        for spec in jobs:
            queue.enqueue(db, *spec)
    return project


def compact(kind, data):
    if kind == "ffi":
        return {
            "field": data["field_id"],
            "role": data["role"],
            "candidates": [
                {k: s[k] for k in ("source_id", "ra_deg", "dec_deg", "disposition")}
                for s in data["sources"]
                if s["candidate"]
            ],
        }
    if kind == "transit":

        def event(s):
            return {k: v for k, v in s.items() if k != "events"}

        return {
            "tic_id": data["tic_id"],
            "role": data["role"],
            "disposition": data["disposition"],
            "ephemeris": data.get("ephemeris"),
            "discovery": event(data.get("discovery", {})),
            "validation": event(data.get("validation", {})),
            "flags": data.get("flags", []),
            "catalog": data["catalog"],
            "eclipsing_binaries": data.get("eclipsing_binaries"),
            "heldout_shapes": data.get("heldout_shapes"),
            "sectors": [s["sector"] for s in data["sources"]],
            "ra_deg": data["ra_deg"],
            "dec_deg": data["dec_deg"],
        }
    if kind == "variability":
        return {k: v for k, v in data.items() if k != "sources"}
    if kind == "motion":
        candidates = [
            {k: v for k, v in s.items() if k != "epoch_measurements"}
            for s in data["candidates"][:8]
        ]
        return {
            "field": data["field"]["id"],
            "counts": data["counts"],
            "candidates": candidates,
        }
    if kind in ("pixels", "variable_pixels"):
        loc = data.get("localization", {})
        return {
            "sector": data["sector"],
            "disposition": data["disposition"],
            "localization": {k: v for k, v in loc.items() if k != "events"},
            "source_sha256": data["source_sha256"],
            "archive_uri": data["archive_uri"],
        }
    if kind == "counterparts":
        return {
            "records": data["records"][:10],
            "total": len(data["records"]),
            "radius_arcsec": data["radius_arcsec"],
            "variable_catalog": {
                k: (
                    {
                        "status": v["status"],
                        "matches": v["matches"][:5],
                        "count": len(v["matches"]),
                    }
                    if isinstance(v, dict) and "matches" in v
                    else v
                )
                for k, v in data["variable_catalog"].items()
            },
        }
    if kind == "plan":
        return {
            k: data[k]
            for k in ("role", "assessment", "falsifier", "actions")
            if k in data
        }
    return data


def options(kind, result, payload):
    actions = {}
    if kind == "transit" and result.get("ephemeris"):
        ephem = {
            k: result["ephemeris"][k]
            for k in ("period_days", "epoch_btjd", "duration_days")
        }
        for p in payload["target"]["products"]:
            actions["pixels_sector_" + str(p["sector"])] = (
                "pixels",
                {
                    **ephem,
                    "tic_id": result["tic_id"],
                    "sector": p["sector"],
                    "product_name": p["filename"],
                    "role": result["role"],
                },
            )
    if kind == "variability" and "epoch_min_btjd" in result:
        for product in payload["target"]["products"]:
            actions["variable_pixels_sector_" + str(product["sector"])] = (
                "variable_pixels",
                {
                    "tic_id": result["tic_id"],
                    "role": result["role"],
                    "sector": product["sector"],
                    "product_name": product["filename"],
                    "epoch_min_btjd": result["epoch_min_btjd"],
                    "period_days": result["period_days"],
                },
            )
    if kind in ("transit", "variability"):
        actions["stellar_counterparts"] = (
            "counterparts",
            {
                **{k: result[k] for k in ("ra_deg", "dec_deg", "tic_id")},
                "period_days": result.get(
                    "period_days", result.get("ephemeris", {}).get("period_days")
                ),
            },
        )
    if kind == "ffi":
        for s in result["sources"]:
            if s["candidate"]:
                actions["counterpart_" + str(s["source_id"])] = (
                    "counterparts",
                    {
                        "ra_deg": s["ra_deg"],
                        "dec_deg": s["dec_deg"],
                        "period_days": s["discovery"]["period_days"],
                    },
                )
    if kind == "motion":
        for i, s in enumerate(result["candidates"][:8]):
            actions[f"counterpart_{i+1}"] = (
                "counterparts",
                {k: s[k] for k in ("ra_deg", "dec_deg")},
            )
    return actions


def signal(kind, r):
    if kind == "ffi":
        return any(s["candidate"] for s in r["sources"])
    if kind == "transit":
        return r["disposition"] in SIGNALS
    if kind == "variability":
        return r["status"] == "two_sector_coherent_variability"
    if kind == "motion":
        return bool(r["candidates"])
    return False


def execute(job, paths):
    kind, p = job["kind"], job["payload"]
    if kind == "ffi":
        from .ffi import run_field

        return run_field(p["field"], paths["ffi"])
    if kind == "transit":
        from .tess import run_target

        r = run_target(p["target"], paths["tess"])
        hosts = p["catalog_targets"]
        tic = str(r["tic_id"])
        r["catalog"] = assess(
            tic,
            r.get("ephemeris", {}).get("period_days"),
            {"targets": {} if hosts is None else {tic: hosts}},
        )
        if signal(kind, r):
            from .tess_eb_catalog import CATALOG_HASH, CATALOG_URL, load_catalog, lookup
            from .tess_morphology import run

            matches = lookup(
                load_catalog(paths["catalog"]), r["tic_id"], r["ra_deg"], r["dec_deg"]
            )
            r["eclipsing_binaries"] = {
                "catalog_url": CATALOG_URL,
                "sha256": CATALOG_HASH,
                "radius_arcsec": 90,
                "matches": matches,
            }
            r["heldout_shapes"] = run(
                {"targets": [p["target"]]}, {"targets": [r]}, paths["tess"]
            )[0]
        return r
    if kind == "variability":
        from .stellar import run_target

        return run_target(p["target"], paths["tess"])
    if kind == "motion":
        from .wise import run_field

        return run_field(p["field"], paths["wise"])
    if kind == "pixels":
        from .pixels import followup

        return followup(p, paths["pixels"])
    if kind == "variable_pixels":
        from .pixels import variable_followup

        return variable_followup(p, paths["pixels"])
    if kind == "counterparts":
        from .variable_catalog import search
        from .stellar import simbad

        match = simbad(p["ra_deg"], p["dec_deg"], paths["catalog"])
        match["variable_catalog"] = search(
            p["ra_deg"],
            p["dec_deg"],
            paths["catalog"],
            radius_arcsec=42,
            tic_id=p.get("tic_id"),
            period_days=p.get("period_days"),
        )
        match["query_position"] = p
        return match
    raise ValueError("Unknown measurement kind")


def step(db, job, paths, project, cap):
    kind, subject, p = job["kind"], job["subject"], job["payload"]
    if kind == "plan":
        actions = p["options"]
        # A zero-budget worker still performs the standard numerical follow-ups.
        use_model = bool(project["model"]) and queue.reserve_model(db, job, cap)
        decision = (
            discovery_agents.decide(project["model"], p["evidence"], list(actions))
            if use_model
            else {
                "role": "deterministic",
                "assessment": "Run available independent checks on the measured signal.",
                "actions": list(actions),
                "model_budget_fallback": bool(project["model"]),
            }
        )
        children = [
            (actions[a][0], subject, actions[a][1]) for a in decision["actions"]
        ]
        if use_model:
            child_ids = [digest(canonical(c).encode()) for c in children]
            children.append(
                (
                    "review",
                    subject,
                    {"parent": p["parent"], "plan": job["id"], "children": child_ids},
                )
            )
        return decision, children
    if kind == "review":
        # Claim only after dependencies are terminal (enforced by worker SQL).
        dependencies = [
            r[0]
            for r in db.execute(
                "SELECT required_id FROM dependency WHERE task_id=? ORDER BY required_id",
                (job["id"],),
            )
        ]
        rows = [
            db.execute(
                "SELECT kind,state,result,error FROM task WHERE id=?", (k,)
            ).fetchone()
            for k in [p["parent"], p["plan"], *dependencies]
        ]
        packet = [
            {
                "kind": r["kind"],
                "evidence": (
                    compact(r["kind"], json.loads(r["result"])) if r["result"] else None
                ),
                "error": r["error"],
            }
            for r in rows
        ]
        if not queue.reserve_model(db, job, cap):
            return {"role": "unreviewed", "reason": "campaign_call_ceiling"}, []
        return discovery_agents.decide(project["model"], packet, [], reviewer=True), []
    result = execute(job, paths)
    loc = result.get("localization", {})
    # Even when bootstrap scatter is large, an independently cataloged variable
    # at a bright displaced peak is a testable alternative to the aperture target.
    displaced = (loc.get("offset_pixels") or 0) >= 1 and (
        loc.get("peak_score") or 0
    ) >= 10
    if (
        kind in ("pixels", "variable_pixels")
        and "centroid_ra_deg" in loc
        and (result["disposition"] == "offset_source_followup" or displaced)
    ):
        source = {
            "ra_deg": loc["centroid_ra_deg"],
            "dec_deg": loc["centroid_dec_deg"],
            "position_kind": "difference_centroid",
            "localization_disposition": result["disposition"],
            "period_days": result["period_days"],
            "sector": result["sector"],
        }
        return result, [("counterparts", subject, source)]
    if not signal(kind, result):
        return result, []
    follow = options(kind, result, p)
    plan = {"parent": job["id"], "options": follow, "evidence": compact(kind, result)}
    return result, [("plan", subject, plan)] if follow else []


def work(root, paths, workers, max_jobs, max_attempts, cap):
    if (
        not 1 <= workers <= 128
        or max_jobs < 1
        or not 1 <= max_attempts <= 10
        or cap < 0
    ):
        raise ValueError("Invalid worker or budget settings")
    path = root / DB
    with closing(queue.connect(path)) as db:
        project = queue.config(db)
    if (
        "runtime_versions" in project
        and project["runtime_versions"] != runtime_versions()
    ):
        raise ValueError(
            "Restore the recorded Python, NumPy and Astropy versions before resuming"
        )
    if project["code_sha256"] != fingerprint():
        raise ValueError(
            "Worker code differs from the campaign. Restore its revision or initialize a new run."
        )
    allowance = max_jobs
    slot_lock = threading.Lock()

    def worker():
        nonlocal allowance
        with closing(queue.connect(path)) as db:
            while True:
                with slot_lock:
                    if allowance <= 0:
                        break
                    job = queue.claim(db, max_attempts)
                    if job is None:
                        break
                    allowance -= 1
                try:
                    with queue.Heartbeat(path, job):
                        result, children = step(db, job, paths, project, cap)
                        queue.complete(db, job, result, children)
                except queue.LostLease:
                    continue
                except Exception as exc:
                    transient = isinstance(
                        exc, (TimeoutError, ConnectionError, urllib.error.URLError)
                    ) and (
                        not isinstance(exc, urllib.error.HTTPError)
                        or exc.code in (408, 429, 500, 502, 503, 504)
                    )
                    try:
                        queue.fail(db, job, exc, max_attempts, transient)
                    except queue.LostLease:
                        pass

    n = min(max_jobs, workers)
    with ThreadPoolExecutor(max_workers=n) as pool:
        futures = [pool.submit(worker) for _ in range(n)]
        for f in futures:
            f.result()
    with closing(queue.connect(path)) as db:
        return queue.stats(db)


def individual_findings(kind, r, subject):
    base = {"measurement": kind, "subject": subject, "claim_status": "screening_only"}
    if kind == "ffi":
        for s in r["sources"]:
            yield {
                **base,
                "object": s["source_id"],
                "role": r["role"],
                "ra_deg": s["ra_deg"],
                "dec_deg": s["dec_deg"],
                "result": s["disposition"],
                "period_days": s["discovery"]["period_days"],
                "field": r["field_id"],
            }
    if kind == "motion":
        for s in r["sources"]:
            yield {
                **base,
                "object": s["source_id"],
                "ra_deg": s["ra_deg"],
                "dec_deg": s["dec_deg"],
                "result": s["disposition"],
                "role": r["field"]["role"],
                "motion_mas_yr": s["motion_mas_yr"],
                "catalog": s.get("gaia"),
                "field": r["field"]["id"],
            }
    elif kind in ("transit", "variability"):
        row = {
            **base,
            "object": "TIC " + str(r["tic_id"]),
            "tic_id": r["tic_id"],
            "ra_deg": r["ra_deg"],
            "dec_deg": r["dec_deg"],
            "role": r["role"],
            "result": r.get("status", r.get("disposition")),
        }
        if kind == "transit":
            row.update(
                period_days=r.get("ephemeris", {}).get("period_days"),
                flags=r.get("flags", []),
                catalog=r.get("catalog"),
                validation_snr=r.get("validation", {}).get("snr"),
            )
        else:
            row.update(
                period_days=r.get("period_days"),
                phase_offset_cycles=r.get("phase_offset_cycles"),
                semiamplitude_ppt=r.get("semiamplitude_ppt"),
                variance_fraction=r.get("variance_fraction"),
            )
        yield row


def export(root, output):
    """Keep separate per-task evidence and a static CSV registry readable on GitHub."""
    import fcntl

    output.mkdir(parents=True, exist_ok=True)
    with (output / ".export.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _export_snapshot(root, output)


def _export_snapshot(root, output):
    import csv
    import io

    output.mkdir(parents=True, exist_ok=True)
    for file in sorted((root / "inputs").glob("*.json")):
        raw = file.read_bytes()
        if digest(raw) != file.stem:
            raise ValueError("Frozen input bytes changed")
        atomic_bytes(output / "inputs" / file.name, raw)
    with closing(queue.connect(root / DB)) as db:
        db.execute("BEGIN")
        project = queue.config(db)
        totals = queue.stats(db)
        jobs = db.execute(
            "SELECT id,kind,subject,state,error,result,attempts,calls,payload FROM task ORDER BY subject,kind,id"
        )
        rows = []
        tasks = []
        for j in jobs:
            audit = {
                k: j[k]
                for k in (
                    "id",
                    "kind",
                    "subject",
                    "state",
                    "attempts",
                    "calls",
                    "error",
                )
            }
            audit.update(
                payload=json.loads(j["payload"]),
                payload_sha256=digest(j["payload"].encode()),
                parents=[
                    r[0]
                    for r in db.execute(
                        "SELECT parent_id FROM lineage WHERE child_id=? ORDER BY parent_id",
                        (j["id"],),
                    )
                ],
                prerequisites=[
                    r[0]
                    for r in db.execute(
                        "SELECT required_id FROM dependency WHERE task_id=? ORDER BY required_id",
                        (j["id"],),
                    )
                ],
            )
            if j["state"] != "done":
                tasks.append(audit)
                continue
            packet = json.loads(j["result"])
            raw = (json.dumps(packet, indent=2, allow_nan=False) + "\n").encode()
            rel = f"evidence/{j['id']}.json"
            atomic_bytes(output / rel, raw)
            tasks.append({**audit, "evidence": rel, "evidence_sha256": digest(raw)})
            for row in individual_findings(j["kind"], packet, j["subject"]):
                rows.append({**row, "task_id": j["id"], "evidence": rel})
        db.commit()
    registry = {
        "schema_version": 1,
        "campaign": project,
        "counts": totals,
        "findings": rows,
        "tasks": tasks,
    }
    atomic_bytes(
        output / "registry.json",
        (json.dumps(registry, indent=2, allow_nan=False) + "\n").encode(),
    )
    buf = io.StringIO(newline="")
    keys = (
        "measurement",
        "object",
        "subject",
        "role",
        "ra_deg",
        "dec_deg",
        "result",
        "period_days",
        "validation_snr",
        "motion_mas_yr",
        "claim_status",
        "task_id",
        "evidence",
    )
    table = csv.DictWriter(buf, keys, extrasaction="ignore", lineterminator="\n")
    table.writeheader()
    table.writerows(rows)
    atomic_bytes(output / "registry.csv", buf.getvalue().encode())
    return {
        "measurements": len(rows),
        "tasks": len(tasks),
        "path": str(output / "registry.json"),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("runs/discovery"))
    commands = p.add_subparsers(dest="command", required=True)
    init = commands.add_parser(
        "init", help="Seed an immutable campaign from one or more archive plans"
    )
    init.add_argument("--tess-plan", type=Path, action="append", default=[])
    init.add_argument("--wise-plan", type=Path, action="append", default=[])
    init.add_argument("--ffi-plan", type=Path, action="append", default=[])
    init.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    init.add_argument("--models", action="store_true")
    init.add_argument("--shards", type=int, default=1)
    init.add_argument("--shard", type=int, default=0)
    w = commands.add_parser(
        "work", help="Claim archive jobs, route measurements and drain a bounded queue"
    )
    w.add_argument("--workers", type=int, default=2)
    w.add_argument("--max-jobs", type=int, default=200)
    w.add_argument("--max-attempts", type=int, default=3)
    w.add_argument(
        "--model-call-ceiling",
        type=int,
        default=0,
        help="Total allowed reservations in this campaign",
    )
    for k in ("tess", "wise", "pixels", "catalog", "ffi"):
        w.add_argument("--" + k + "-cache", type=Path)
    commands.add_parser("status")
    retry = commands.add_parser(
        "retry", help="Explicitly reset the attempt count on a failed task"
    )
    retry.add_argument("task_id")
    show = commands.add_parser(
        "show", help="Read full evidence and execution state for a task"
    )
    show.add_argument("task_id")
    e = commands.add_parser("export")
    e.add_argument("output", type=Path)
    a = p.parse_args(argv)
    try:
        if a.command == "init":
            tess, wise = a.tess_plan, a.wise_plan
            if not tess and not wise and not a.ffi_plan:
                tess, wise = [DEFAULT_TESS], [DEFAULT_WISE]
            result = initialize(
                a.root,
                tess,
                wise,
                a.catalog,
                discovery_agents.settings() if a.models else None,
                a.shards,
                a.shard,
                a.ffi_plan,
            )
        elif a.command == "work":
            paths = {
                k: getattr(a, k + "_cache") or a.root / "cache" / k
                for k in ("tess", "wise", "pixels", "catalog", "ffi")
            }
            result = work(
                a.root,
                paths,
                a.workers,
                a.max_jobs,
                a.max_attempts,
                a.model_call_ceiling,
            )
        elif a.command == "export":
            result = export(a.root, a.output)
        else:
            with closing(queue.connect(a.root / DB)) as db:
                if a.command == "status":
                    result = queue.stats(db)
                elif a.command == "show":
                    row = db.execute(
                        "SELECT * FROM task WHERE id=?", (a.task_id,)
                    ).fetchone()
                    if not row:
                        raise ValueError("Unknown task")
                    result = {
                        **dict(row),
                        "payload": json.loads(row["payload"]),
                        "result": json.loads(row["result"]) if row["result"] else None,
                    }
                else:
                    result = queue.retry(db, a.task_id)
        print(json.dumps(result, indent=2, allow_nan=False))
        return int(a.command == "work" and result["tasks"].get("failed", 0) > 0)
    except (ValueError, KeyError, OSError, sqlite3.Error) as exc:
        print(f"discovery: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
