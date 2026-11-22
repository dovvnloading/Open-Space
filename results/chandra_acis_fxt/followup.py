"""Selectively repair dated-survey follow-up after the primary Actions deposit."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from unittest.mock import patch

import archive
import catalog_pack
import campaign
from migrate import BEFORE, INVENTORY, preserve_characterization, verify_scan

OLD = "76ddfd371cf8f02bd54c3653547f12c74a1a23d09fa7e0bbd2fa933ad5f79d8a"
NEW = "9f3615c635e42d60425be4738dd3ffd62dfe8aee74bc717736df951af2a996e8"
COMPANION = "7c945742d96a2f6918ce4028d56afffc967ab86af1fe1fa918a53ca9e9a0991b"
UNKNOWN = "No reported object types in a nonempty NED cone"
GUARD = ('    if label == "ned" and rows and not any(row["pretype"] for row in rows):\n'
         '        raise ValueError("' + UNKNOWN + '")\n').encode()
ANCHOR = b'    rows = [{k: scalar(row[columns[k]]) for k in required} for row in table]\n'
VETO = {"foreground_veto", "known_transient"}
TERMINAL = VETO | {"ambiguous", "no_counterpart", "host_candidate"}
SCOPED = Path("history/ned-untyped-v1")
FOLLOW = {"obsid", "catalog_cache", "vet", "flux", "diagnostics", "trigger_key"}


def sha(p):
    with p.open("rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def value_sha(item):
    return hashlib.sha256(json.dumps(item, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def safe(root, name):
    path = (root / name).resolve()
    if not path.is_relative_to(root) or path == root:
        raise ValueError(f"Path escapes evidence root: {name}")
    return path


def save(p, item):
    p.parent.mkdir(parents=True, exist_ok=True)
    archive.dump(p, item)


def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    temp = dst.with_name(dst.name + ".new")
    shutil.copyfile(src, temp)
    temp.replace(dst)


def old_path(root, obsid):
    return root / "observations" / f"{obsid}.json.gz"


def source_pair(bundle):
    stock = (campaign.HERE / "vet.py").read_bytes()
    current = hashlib.sha256(stock).hexdigest()
    if current == NEW and stock.count(ANCHOR) == 1:
        old = stock.replace(ANCHOR, ANCHOR + GUARD)
    elif current == OLD:
        old = stock
    else:
        raise ValueError("The original or reviewed vet.py is required")
    new = old.replace(ANCHOR + GUARD, ANCHOR)
    if (old.count(GUARD) != 1 or hashlib.sha256(old).hexdigest() != OLD
            or hashlib.sha256(new).hexdigest() != NEW):
        raise ValueError("The reviewed two-line deletion could not be verified")
    bundle.mkdir(parents=True, exist_ok=True)
    for label, data in (("before", old), ("after", new)):
        file = bundle / f"vet.{label}.py"
        if file.exists() and file.read_bytes() != data:
            raise ValueError("Conflicting stored source in repair bundle")
        file.write_bytes(data)


def checked(root, bundle):
    expected = {**BEFORE, "vet.py": OLD}
    if (sha(campaign.HERE / "vet.py") not in {OLD, NEW} or {**campaign.identity(), "vet.py": OLD} != expected
            or sha(campaign.HERE / "replay_untyped_ned.py") != COMPANION):
        raise ValueError("An unrelated scientific implementation changed")
    freeze = campaign.read(root / "inventory/campaign.json")
    if (freeze["implementation"] not in (expected, {**expected, "vet.py":NEW}) or freeze["inventory_sha256"] != INVENTORY
            or sha(root / "inventory/inventory.json") != INVENTORY
            or freeze["n_controls"] != 24 or freeze["n_injections_per_product"] != 12):
        raise ValueError("Unexpected primary inventory or science settings")
    if freeze["implementation"]["vet.py"] == NEW:
        history = safe(root, freeze.get("followup_revision", ""))
        if (not history.is_relative_to(root / SCOPED) or sha(history / "vet.before.py") != OLD
                or sha(history / "vet.after.py") != NEW or not (history / "result.json").is_file()
                or sha(history / "campaign.before.json") != freeze.get("prior_campaign_sha256")
                or campaign.read(history / "result.json")["plan_sha256"] != sha(history / "plan.json")):
            raise ValueError("Continuation lacks the prior completed revision and its source evidence")
    if sha(bundle / "vet.before.py") != OLD or sha(bundle / "vet.after.py") != NEW:
        raise ValueError("Repair source hash mismatch")
    return freeze


def load_vet(path):
    spec = importlib.util.spec_from_file_location("vet", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["vet"] = module
    spec.loader.exec_module(module)
    return module


def search_state(product, oid):
    scan = product.get("search")
    if not scan:
        return "absent"
    if (scan.get("obsid") is not None and int(scan["obsid"]) != oid
            or scan.get("seed") is not None and scan["seed"] != int(hashlib.sha256(
                ("chandra-acis-fxt-v1:" + product["input"]["sha256"]).encode()).hexdigest()[:16], 16) >> 1):
        return "conflicting"
    try:
        verify_scan(scan, product["input"], oid)
        for trial in (scan, *scan["controls"]):
            if not isinstance(trial["triggers"], list):
                return "partial"
    except (ValueError, KeyError, TypeError):
        return "partial"
    return "complete"


def triggers(product):
    scan = product.get("search")
    if scan:
        yield from (("real", i, t) for i, t in enumerate(scan.get("triggers", [])))
        for c in scan.get("controls", []):
            yield from ((f"control-{c['replica']}", i, t) for i, t in enumerate(c.get("triggers", [])))


def affected(v):
    return (v.get("errors", {}).get("ned") == UNKNOWN or
            any(q.get("query_error") == UNKNOWN for q in v.get("queries", {}).values()))


def certify(v):
    if v.get("status") not in TERMINAL or v.get("errors") != {}:
        return False
    # A successful NED decode in the old pinned program proves the removed
    # guard was never reached. Check the primary and mirror dispositions too.
    if affected(v):
        return False
    q = v.get("queries", {})
    routes = (("gaia", "gaia_aip"), ("simbad", "simbad_idia"), ("ned", "ned_vo"),
              ("quirola2023", "quirola2023_idia"), ("dillmann2025", "dillmann2025_idia"))
    return (all(any(q.get(k, {}).get("ok") and isinstance(q[k].get("row_count"), int)
                    and q[k]["row_count"] >= 0 and q[k].get("sha256") for k in alternatives)
                for alternatives in routes)
            and q.get("cudidi2026", {}).get("ok") and bool(q["cudidi2026"].get("sha256")))


def inspect(oid, row, original_root=None, worker=None):
    reasons = set()
    work = {"real": 0, "control": 0, "ned_real": 0, "ned_control": 0}
    if not row:
        return ["unattempted"], work
    if row["obsid"] != oid:
        raise ValueError(f"Ledger contents do not match ObsID {oid}")
    if row["status"] != "complete":
        reasons.add(row["status"])
    expected = [(p["filename"], p["sha256"]) for p in row.get("acquired", [])]
    found = [(p["input"]["filename"], p["input"]["sha256"]) for p in row.get("products", [])]
    if not expected:
        reasons.add("no_complete_acquisition_manifest")
    if len(expected) != len(set(expected)) or len(found) != len(set(found)):
        raise ValueError(f"Duplicate product in ObsID {oid}")
    if set(expected) != set(found):
        reasons.add("partial_product_list")
    for p in row.get("products", []):
        if not p.get("complete") or not p.get("search"):
            reasons.add("partial_search")
        state = search_state(p, oid)
        if state != "complete":
            reasons.add(state + "_search")
            continue
        for stream, i, event in triggers(p):
            kind = "real" if stream == "real" else "control"
            work[kind] += 1
            if not event.get("trigger_key"):
                reasons.add("partial_trigger_annotation")
            if event.get("trigger_key") not in {None, f"{stream}/{i}"}:
                raise ValueError(f"Wrong trigger identity in {oid} {stream}/{i}")
            v = event.get("vet", {})
            if original_root is not None and row["status"] == "complete":
                try:
                    cache_root = original_root
                    if worker:
                        base, namespace_name = worker
                        path = Path(event.get("catalog_cache") or row["catalog_cache"])
                        if path.parts[:2] == ("catalogs", namespace_name):
                            cache_root = base
                    folder = old_cache(cache_root, event, row)
                    verify_original(v, folder, cache_root)
                    if not memo_matches(folder, v, cache_root):
                        reasons.add("catalog_memo_requires_reconstruction")
                except (KeyError, OSError, ValueError):
                    reasons.add("missing_or_altered_original_catalog_inputs")
            if affected(v):
                work["ned_" + kind] += 1
                reasons.add("all_blank_ned_type")
            if not certify(v):
                reasons.add("catalog_requires_review")
            if (kind == "real" and event["instrument_pass"] and v.get("status") not in VETO
                    and ("peak_flux_erg_cm2_s" not in event.get("flux", {})
                         or not event.get("diagnostics", {}).get("files"))):
                reasons.add("physical_characterization")
    return sorted(reasons), work


def prepare(root, bundle, groups):
    if bundle.exists() and any(bundle.iterdir()):
        raise ValueError("Use an empty directory for a new plan")
    source_pair(bundle)
    freeze = checked(root, bundle)
    plan, tasklist = {}, []
    inventory = campaign.read(root / "inventory/inventory.json")
    for entry in inventory:
        oid = entry["obsid"]
        path = old_path(root, oid)
        row = campaign.read_obs(path) if path.is_file() else None
        reasons, counts = inspect(oid, row, root)
        plan[str(oid)] = {"primary_sha256": sha(path) if path.is_file() else None, "reasons": reasons,
                          "primary_status": row["status"] if row else "unattempted", "rows_checked": counts}
        if reasons:
            tasklist.append((oid, entry.get("event_count") or 0))
    if len(plan) != len(inventory):
        raise ValueError("Duplicate ObsID in original census")
    if not 1 <= groups <= 192:
        raise ValueError("Use 1 through 192 groups")
    chunks = [[] for _ in range(min(groups, len(tasklist)))]
    weights = [0] * len(chunks)
    for oid, count in sorted(tasklist, key=lambda pair: (-pair[1], pair[0])):
        i = min(range(len(chunks)), key=weights.__getitem__)
        chunks[i].append(oid)
        weights[i] += max(count, 1000)
    copy(root / "inventory/campaign.json", bundle / "campaign.before.json")
    copy(campaign.HERE / "replay_untyped_ned.py", bundle / "companion.py")
    save(bundle / "plan.json", {"obsids": plan, "groups": chunks, "inventory_sha256": INVENTORY,
         "before": freeze["implementation"], "after": {**freeze["implementation"], "vet.py": NEW},
         "companion_sha256": COMPANION, "repair_sha256": sha(Path(__file__)),
         "catalog_reader_sha256": sha(Path(catalog_pack.__file__)),
         "primary_git_commit": os.environ["GITHUB_SHA"], "run_id": os.environ["GITHUB_RUN_ID"],
         "run_attempt": os.environ["GITHUB_RUN_ATTEMPT"]})
    if os.getenv("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as out:
            print("matrix=" + json.dumps(list(range(len(chunks)))), file=out)
    print(f"{len(plan)} census ObsIDs audited; {len(tasklist)} require work, in {len(chunks)} groups.")


def primary(root, plan, oid):
    path = old_path(root, oid)
    if (sha(path) if path.is_file() else None) != plan["obsids"][str(oid)]["primary_sha256"]:
        raise ValueError(f"Primary ObsID {oid} changed after audit")
    return campaign.read_obs(path) if path.is_file() else None


def old_cache(root, event, ledger):
    return safe(root, event.get("catalog_cache") or ledger["catalog_cache"])


def memo_matches(folder, vetting, root=None):
    for route, q in vetting.get("queries", {}).items():
        if not q.get("url"):
            continue
        key = hashlib.sha256(q["url"].encode()).hexdigest()[:24]
        path = folder / f"{route}-{key}.json"
        if not path.is_file() or (root is not None and not catalog_pack.memo_unchanged(root, path)):
            return False
        try:
            current = campaign.read(path)
        except (OSError, ValueError):
            return False
        if any(current.get(k) != q.get(k) for k in ("url", "ok", "sha256", "body_file", "size_bytes")):
            return False
    return True


def request_files(v, folder):
    """The original raw query receipt and checksum-qualified body selections."""
    for route, q in v.get("queries", {}).items():
        if not q.get("url"):
            continue
        key = hashlib.sha256(q["url"].encode()).hexdigest()[:24]
        memo = safe(folder, f"{route}-{key}.json")
        if q.get("ok") and not q.get("body_file"):
            raise ValueError(f"Successful {route} has no original response body")
        bodies = {}
        for item in [q, *q.get("attempts", [])]:
            name = item.get("body_file")
            if not name:
                continue
            if Path(name).name != name or not name.endswith(".bin"):
                raise ValueError("Unexpected original response name")
            body = safe(folder, name)
            if body in bodies and any(bodies[body][k] != item[k] for k in ("sha256", "size_bytes")):
                raise ValueError("Inconsistent original HTTP attempt byte identities")
            bodies[body] = item
        yield route, q, memo, bodies


def verify_original(v, folder, root):
    for _, _, _, bodies in request_files(v, folder):
        for body, receipt in bodies.items():
            catalog_pack.recorded_reply(root, body, receipt)


def companion_run(root, ledger, directory, old_vet):
    load_vet(old_vet)
    import replay_untyped_ned as companion
    path = old_path(root, ledger["obsid"])
    view = directory / "view"
    copy(path, old_path(view, ledger["obsid"]))
    for p in ledger.get("products", []):
        for _, _, t in triggers(p):
            v = t.get("vet", {})
            if v.get("errors", {}).get("ned") != UNKNOWN:
                continue
            folder = old_cache(root, t, ledger)
            for _, _, memo, bodies in request_files(v, folder):
                if memo.is_file():
                    copy(memo, view / memo.relative_to(root))
                for file, receipt in bodies.items():
                    catalog_pack.recorded_reply(root, file, receipt, view / file.relative_to(root))
    report = directory / "offline-companion.json"
    # The function reports individual failed/missing offline inputs in its JSON.
    companion.replay(view, report)
    return campaign.read(report)


def scientific(item):
    """The companion explicitly uses None for either representation of blank."""
    if isinstance(item, dict):
        return {k: None if k == "pretype" and v == "" else scientific(v) for k, v in item.items()}
    if isinstance(item, list):
        return [scientific(v) for v in item]
    return item


def seed_cache(root, row, event, target, *, live):
    original = event.get("vet", {})
    if not original.get("queries"):
        return {}
    folder = old_cache(root, event, row)
    target.mkdir(parents=True, exist_ok=True)
    import vet
    failed = {route for label, error in original.get("errors", {}).items() if error != UNKNOWN
              for route in (label, vet.MIRRORS.get(label, (label,))[0])}
    pins = {}
    for route, q, file, bodies in request_files(original, folder):
        if live and route in failed and (not q.get("ok") or q.get("query_error") != UNKNOWN):
            continue
        for body, receipt in bodies.items():
            catalog_pack.recorded_reply(root, body, receipt, target / body.name)
        # Reconstruct this exact historical response selection if a shared cache
        # manifest subsequently advanced to another successful response.
        memo = {k: v for k, v in q.items() if k not in {"row_count", "query_error"}}
        save(target / file.name, memo)
        if q.get("ok"):
            pins[route] = (q["url"], q["sha256"])
    return pins


def require_pins(result, pins):
    for route, (url, checksum) in pins.items():
        receipt = result.get("queries", {}).get(route)
        if receipt is None:
            # If the primary now decodes, the original fallback is never opened.
            if (route == "ned_vo" and "ned" in pins and
                    result.get("queries", {}).get("ned", {}).get("row_count") is not None):
                continue
            raise ValueError(f"Previously successful {route} receipt not accounted for")
        if not receipt.get("ok") or receipt.get("sha256") != checksum or receipt.get("url") != url:
            raise ValueError(f"Previously successful {route} bytes changed")


def event_names(filename, oid):
    if Path(filename).name != filename or not filename.startswith("acis"):
        raise ValueError("Invalid archived product filename")
    return [f"{base}/{oid % 10}/{oid}/primary/{filename}" for base in archive.ROOTS]


def event_bytes(p, oid, directory, retained, attempts):
    """Verify a retained copy or obtain the same named revision and bytes."""
    name = p["filename"]
    urls = event_names(name, oid)
    pinned = p.get("sha256")
    if p.get("url"):
        archive.check_archive_url(p["url"])
        urls.insert(0, p["url"])
    target = directory / name
    sources = []
    if pinned and retained:
        for base in retained:
            sources.extend(("file", q) for q in (base / str(oid) / name,
                                                  *base.glob(f"*/{oid}/{name}")) if q.is_file())
    sources.extend(("url", url) for url in dict.fromkeys(urls))
    for kind, src in sources:
        record = {"obsid": oid, "filename": name, "source": str(src), "kind": kind}
        attempts.append(record)
        try:
            if kind == "file":
                copy(src, target)
                rec = {"sha256": sha(target), "bytes": target.stat().st_size, "url": p["url"]}
            else:
                rec = archive.download(src, target)
            record["response"] = rec
            if pinned and (rec["sha256"] != pinned or rec["bytes"] != p["bytes"]):
                raise ValueError("Archived input differs from the original bytes")
            fits = archive.inspect_product(target, oid)
            if p.get("fits_sha256") and fits["fits_sha256"] != p["fits_sha256"]:
                raise ValueError("Decompressed FITS hash changed")
            record["fits_sha256"] = fits["fits_sha256"]
            return {**rec, **fits, **p, "filename": name, "obsid": oid, "path": str(target)}
        except (OSError, ValueError, EOFError, KeyError, archive.http.client.HTTPException) as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
            target.unlink(missing_ok=True)
    raise ValueError(f"No verified archive or retained copy of {name}")


def partial_getter(prior, obsid, retained, attempts):
    old_indexes, pinned = set(), {}
    def add(p):
        name = p["filename"]
        event_names(name, obsid)
        if name in pinned and (pinned[name]["sha256"], pinned[name]["bytes"]) != (p["sha256"], p["bytes"]):
            raise ValueError(f"Conflicting primary byte receipts for {name}")
        pinned[name] = {**pinned.get(name, {}), **p}
    for a in prior.get("archive_attempts", []):
        if a.get("selected_files"):
            old_indexes.add(tuple(sorted(a["selected_files"])))
        for p in a.get("downloaded", []):
            add(p)
    for p in prior.get("products", []):
        add(p["input"])
    if not pinned:
        return archive.fetch_events
    if len(old_indexes) != 1:
        raise ValueError("Partial acquisition has no uniquely witnessed original product list")
    names, = old_indexes
    if not set(pinned) <= set(names) or len(set(names)) != len(names):
        raise ValueError("Incomplete or contradictory original product membership")
    def fetch(oid, directory):
        if oid != obsid:
            raise ValueError("Unexpected ObsID in partial retrieval")
        return [event_bytes(pinned.get(name, {"filename": name}), oid, directory, retained, attempts)
                for name in names]
    return fetch


def raw_scan(product):
    if "search" not in product:
        return None
    from copy import deepcopy
    cloned = deepcopy(product)
    for _, _, event in triggers(cloned):
        for key in FOLLOW:
            event.pop(key, None)
    return value_sha(cloned["search"])


def prior_preserved(old, new):
    if not old:
        return
    if old.get("acquired") and new.get("acquired") != old["acquired"]:
        raise ValueError("Previously acquired event manifest changed")
    after = {p["input"]["sha256"]: p for p in new.get("products", [])}
    for p in old.get("products", []):
        if search_state(p, old["obsid"]) == "conflicting":
            raise ValueError("The original search identity conflicts with its product")
        if search_state(p, old["obsid"]) != "complete":
            continue
        new_product = after.get(p["input"]["sha256"])
        if not new_product or raw_scan(p) != raw_scan(new_product):
            raise ValueError("Existing event search, seeds, control rows, or injection outcomes changed")
        for (_, _, before), (_, _, now) in zip(triggers(p), triggers(new_product), strict=True):
            if certify(before.get("vet", {})) and before["vet"] != now.get("vet"):
                raise ValueError("An unaffected, previously complete classification changed")


def prepare_event(row, p, stream, i, event, root, out, corrected, checked_results, changes):
    old = event.get("vet", {})
    affected_before = affected(old)
    dubious = not certify(old)
    rebuild_memo = bool(old.get("queries")) and not memo_matches(old_cache(root, event, row), old, root)
    if not (affected_before or dubious or rebuild_memo):
        return
    name = f"{stream}/{i}"
    rel = Path("catalogs") / namespace() / str(row["obsid"]) / p["input"]["sha256"] / stream / str(i)
    cache = out / rel
    pins = seed_cache(root, row, event, cache, live=False)
    result = corrected.vet_trigger(event, cache, replay=True) if old.get("queries") else None
    if result is not None:
        require_pins(result, pins)
    key = (p["input"]["sha256"], name)
    check = checked_results.get(key, {})
    if check.get("result") and not check["result"].get("errors"):
        fields = ("status", "host", "matches", "position", "errors")
        if not result or any(scientific(result[k]) != scientific(check["result"][k]) for k in fields):
            raise ValueError(f"Two offline implementations disagree for {name}")
    if result and not result["errors"] and result["status"] in TERMINAL:
        if rebuild_memo and not affected_before and not dubious:
            if any(result[k] != old[k] for k in ("status", "host", "matches", "position", "errors")):
                raise ValueError("Reconstructed historical cache does not reproduce the complete classification")
            event["vet"] = old
        else:
            event["vet"] = result
        event["catalog_cache"] = str(rel)
    else:
        # A live retry may refresh only historically failed logical services.
        # Successful response bytes are seeded separately for every trigger.
        if cache.exists():
            shutil.rmtree(cache)
        pins = seed_cache(root, row, event, cache, live=True)
        event["vet"] = {"status": "catalog_incomplete", "reason": "Residual query requires a recorded live retry"}
    verified_failures = [route for route, receipt in old.get("queries", {}).items()
                         if receipt.get("url") and not receipt.get("ok") and
                         any(route == logical and old["queries"].get(mirror[0], {}).get("ok")
                             and logical not in old.get("errors", {}) for logical, mirror in corrected.MIRRORS.items())]
    changes[key] = {"cache": str(rel), "replay_only_failed_primaries": verified_failures, "pins": pins, "previous_status": old.get("status"),
                    "offline_status": result["status"] if result else None,
                    "offline_errors": result["errors"] if result else None}
    if stream == "real" and event["vet"]["status"] in VETO:
        for k in ("flux", "diagnostics"):
            if event.get(k, {}).get("error"):
                del event[k]


def namespace():
    return "followup-" + os.environ["GITHUB_RUN_ID"] + "-" + os.environ["GITHUB_RUN_ATTEMPT"]


def run_one(root, bundle, plan, oid, out, retained):
    checked(root, bundle)
    meta = plan["obsids"][str(oid)]
    if not meta["reasons"]:
        raise ValueError(f"ObsID {oid} is outside the residual plan")
    prior = primary(root, plan, oid)
    out.mkdir(parents=True, exist_ok=True)
    status_path = out / "repair_records" / f"{oid}.json"
    if status_path.is_file():
        receipt = campaign.read(status_path)
        hashed = receipt.get("revised_sha256")
        if (not hashed or sha(old_path(out, oid)) == hashed) and receipt["plan_sha256"] == sha(bundle / "plan.json"):
            return receipt
        raise ValueError("Existing worker result differs from the registered plan")
    audit = {"obsid": oid, "plan_sha256": sha(bundle / "plan.json"),
             "primary_sha256": meta["primary_sha256"], "reasons": meta["reasons"],
             "run_id": os.environ["GITHUB_RUN_ID"], "run_attempt": os.environ["GITHUB_RUN_ATTEMPT"],
             "event_byte_attempts": [], "outcome": "unresolved", "revised_sha256": None}
    changes = {}
    new_path = old_path(out, oid)
    leftovers = [p for p in (new_path, out / "catalogs" / namespace() / str(oid),
                             out / "photometry" / str(oid)) if p.exists()]
    if leftovers:
        quarantine = out / "repair_records" / str(oid) / "interrupted"
        quarantine.mkdir(parents=True, exist_ok=True)
        attempt = Path(tempfile.mkdtemp(prefix="worker-", dir=quarantine))
        for previous in leftovers:
            moved = attempt / previous.relative_to(out)
            moved.parent.mkdir(parents=True, exist_ok=True)
            previous.rename(moved)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            temp = Path(tmp)
            stage = temp / "prior"
            raw = temp / "events"
            raw.mkdir()
            # A separate old-ledger view includes historical bytes; no original
            # observation, catalog response or photometry file is edited.
            original = campaign.read_obs(old_path(root, oid)) if prior else {}
            for product in original.get("products", []):
                state = search_state(product, oid)
                if state == "conflicting":
                    raise ValueError("The original search identity conflicts with its product")
                if state == "partial":
                    product.pop("search")
                    product["complete"] = False
            relevant = any(e.get("vet", {}).get("errors", {}).get("ned") == UNKNOWN
                           for p in original.get("products", []) for _, _, e in triggers(p))
            comparison = {}
            if relevant:
                report = companion_run(root, original, temp, bundle / "vet.before.py")
                copy(temp / "offline-companion.json", out / "repair_records" / str(oid) / "offline-companion.json")
                if (report["vet_sha256"] != OLD or report["replay_implementation_sha256"] != COMPANION
                        or any(t["ledger_sha256"] != meta["primary_sha256"] for t in report["triggers"])):
                    raise ValueError("Companion replay used different original source or ledger bytes")
                comparison = {(t["input_sha256"], t["trigger_key"]):t for t in report["triggers"]}
                if len(comparison) != len(report["triggers"]):
                    raise ValueError("Ambiguous offline companion trigger identity")
            corrected = load_vet(bundle / "vet.after.py")
            for p in original.get("products", []):
                for stream, i, event in triggers(p):
                    event["obsid"] = oid
                    if event.get("vet", {}).get("queries"):
                        verify_original(event["vet"], old_cache(root, event, original), root)
                    prepare_event(original, p, stream, i, event, root, out, corrected, comparison, changes)
                    if stream == "real" and event["instrument_pass"] and event.get("vet", {}).get("status") not in VETO:
                        preserve_characterization(event, root, out, p["input"], stream, i)
            if prior:
                ppath = old_path(stage, oid)
                ppath.parent.mkdir(parents=True, exist_ok=True)
                campaign.save_obs(ppath, original)
            cache = out / "catalogs" / namespace() / str(oid) / "fresh"
            cache.mkdir(parents=True, exist_ok=True)
            capture = audit["event_byte_attempts"]
            def replay(manifest, event_dir):
                for p in campaign.read(manifest):
                    event_bytes(p, oid, event_dir, retained, capture)
            fetch = archive.fetch_events if original.get("acquired") else partial_getter(original, oid, retained, capture)
            by_filename = {p["input"]["filename"]:p["input"]["sha256"] for p in original.get("products", [])}
            annotate = campaign.annotate
            validation_errors = []
            def annotate_pinned(trigger, directory, event_file, stream, index, default_cache, output):
                item = changes.get((by_filename.get(event_file.name), f"{stream}/{index}"))
                selected = out / item["cache"] if item else default_cache
                selected.mkdir(parents=True, exist_ok=True)
                read_input = corrected.read_input
                def pinned_primary(label, url, cache, replay, records, params=None):
                    locked = item and label in item["replay_only_failed_primaries"]
                    return read_input(label, url, cache, bool(replay or locked), records, params)
                with patch.object(corrected, "read_input", pinned_primary):
                    annotate(trigger, directory, event_file, stream, index, selected, output)
                if item:
                    try:
                        require_pins(trigger["vet"], item["pins"])
                    except ValueError as exc:
                        validation_errors.append(str(exc))
                        raise
                if stream == "real" and trigger.get("vet", {}).get("status") in VETO:
                    for k in ("flux", "diagnostics"):
                        if trigger.get(k, {}).get("error"):
                            del trigger[k]
            with (patch.object(campaign, "EVIDENCE", stage), patch.object(campaign, "annotate", annotate_pinned),
                  patch.object(campaign, "replay", replay), patch.object(campaign, "fetch_events", fetch)):
                campaign.process(oid, out, raw, cache)
            if validation_errors:
                raise ValueError("Catalog provenance failed: " + " | ".join(validation_errors))
            if not new_path.exists():
                raise ValueError("Worker produced no observation ledger")
            revised = campaign.read_obs(new_path)
            prior_preserved(prior, revised)
            reasons, counts = inspect(oid, revised, root, (out, namespace()))
            # Errors during input replay can leave campaign.process holding the
            # prior reference. Promotion requires a successful byte verification.
            if prior and prior.get("acquired"):
                verified = {a["filename"] for a in capture if a.get("fits_sha256") and not a.get("error")}
                if verified != {p["filename"] for p in prior["acquired"]}:
                    raise ValueError("Not all originally acquired EVT2 bytes could be verified")
            audit.update(revised_sha256=sha(new_path), revised_status=revised["status"],
                         remaining=reasons, revised_rows_checked=counts,
                         outcome="complete" if not reasons else "unresolved")
    except Exception as exc:
        audit["error"] = f"{type(exc).__name__}: {exc}"
        # Quarantine an unsuccessfully validated output for diagnosis. Its
        # contents must not replace the primary observation in final reduction.
        if new_path.exists():
            copy(new_path, out / "repair_records" / str(oid) / "quarantined.json.gz")
            new_path.unlink()
    audit["catalog_changes"] = [{"product_sha256": h, "trigger_key": key, **v} for (h,key),v in changes.items()]
    save(status_path, audit)
    print(f"ObsID {oid}: {audit['outcome']}, {audit.get('error', audit.get('remaining'))}", flush=True)
    return audit


def deferred(old, why, root=None):
    """Make unverified prior classifications visible to the corrected reducer."""
    if not old:
        return None
    expected = {(p["filename"], p["sha256"]) for p in old.get("acquired", [])}
    scanned = set()
    for p in old["products"]:
        if search_state(p, old["obsid"]) != "complete":
            p["complete"] = False
            continue
        if p.get("complete"):
            scanned.add((p["input"]["filename"], p["input"]["sha256"]))
        for _, _, t in triggers(p):
            valid = certify(t.get("vet", {}))
            if valid and root:
                try:
                    folder = old_cache(root, t, old)
                    verify_original(t["vet"], folder, root)
                    valid = memo_matches(folder, t["vet"], root)
                except (KeyError, OSError, ValueError):
                    valid = False
            if not valid:
                previous = t.get("vet", {})
                t["vet"] = {**previous, "status": "catalog_incomplete", "reason": why,
                            "previous_status": previous.get("previous_status", previous.get("status"))}
    if old["status"] in ("complete", "followup_incomplete"):
        old["status"] = "followup_incomplete" if expected and expected == scanned else "error"
    old.setdefault("errors", []).append(why)
    return old


def archive_replace(src, dest, before):
    if dest.is_file():
        if sha(src) == sha(dest):
            return
        if before.is_file() and sha(before) != sha(dest):
            raise ValueError(f"Conflicting existing historical snapshot {before}")
        if not before.exists():
            copy(dest, before)
    copy(src, dest)


def finish(root, bundle, staged, install):
    plan = campaign.read(bundle / "plan.json")
    revision = root / SCOPED / f"{plan['run_id']}-{plan['run_attempt']}"
    plan_hash = sha(bundle / "plan.json")
    final = revision / "result.json"
    if final.is_file():
        result = campaign.read(final)
        if (result["plan_sha256"] != plan_hash or sha(campaign.HERE / "vet.py") != NEW
                or campaign.read(root / "inventory/campaign.json")["implementation"] != plan["after"]):
            raise ValueError("A changed final repair receipt already exists")
        return result
    freeze_file = root / "inventory/campaign.json"
    # A previous interrupted finalization may already have installed the exact
    # replacement source; the freeze is written last after the results ledger.
    freeze = campaign.read(bundle / "campaign.before.json")
    if freeze["implementation"] != plan["before"] or freeze["inventory_sha256"] != INVENTORY:
        raise ValueError("Frozen pre-repair receipt does not match plan")
    if (sha(bundle / "companion.py") != COMPANION or sha(bundle / "vet.before.py") != OLD
            or sha(bundle / "vet.after.py") != NEW or sha(root / "inventory/inventory.json") != INVENTORY
            or plan["repair_sha256"] != sha(Path(__file__)) or plan["catalog_reader_sha256"] != sha(Path(catalog_pack.__file__))
            or {**campaign.identity(), "vet.py":OLD} != {**BEFORE, "vet.py": OLD}):
        raise ValueError("Repair sources or inventory changed")
    if campaign.read(freeze_file)["implementation"] not in (plan["before"], plan["after"]):
        raise ValueError("Unknown current campaign identity")
    if sha(campaign.HERE / "vet.py") != NEW and not install:
        raise ValueError("Finish requires the reviewed corrected vet.py or the --install-reviewed-vet flag")
    revision.mkdir(parents=True, exist_ok=True)
    for name in ("plan.json", "vet.before.py", "vet.after.py", "companion.py", "campaign.before.json"):
        dest = revision / name
        if dest.is_file() and sha(dest) != sha(bundle / name):
            raise ValueError(f"Different prior repair artifact {name}")
        copy(bundle / name, dest)
    expected = {int(oid) for oid, status in plan["obsids"].items() if status["reasons"]}
    records = {}
    for p in staged.rglob("repair_records/*.json"):
        r = campaign.read(p)
        oid = r["obsid"]
        if oid not in expected or p.name != f"{oid}.json" or oid in records or r["plan_sha256"] != plan_hash:
            raise ValueError(f"Unrecognized or duplicate worker result {p}")
        records[oid] = (r, p.parent.parent)
    results = {}
    for key, m in plan["obsids"].items():
        oid = int(key)
        dest = old_path(root, oid)
        before = old_path(revision / "before", oid)
        original_path = before if before.is_file() else dest
        if m["primary_sha256"] is None:
            original = None
            if dest.exists() and not (revision / "new-ledgers" / (key + ".json")).is_file():
                raise ValueError(f"Unexpected primary ledger appeared: {oid}")
        else:
            if not original_path.is_file() or sha(original_path) != m["primary_sha256"]:
                raise ValueError(f"Primary ledger changed: {oid}")
            original = campaign.read_obs(original_path) if oid in expected else None
        if oid not in expected:
            discrepancies, _ = inspect(oid, campaign.read_obs(original_path), root)
            if discrepancies:
                raise ValueError(f"Carried observation changed after planning: {oid} {discrepancies}")
            continue
        if original is not None and not before.exists():
            copy(original_path, before)
        r, job = records.get(oid, ({}, None))
        if r.get("primary_sha256", m["primary_sha256"]) != m["primary_sha256"]:
            raise ValueError(f"Worker used another primary ledger for {oid}")
        reason = r.get("error") or "Residual worker missing or its output could not be validated"
        revised = None
        revised_path = old_path(job, oid) if job else None
        approved = bool(r.get("revised_sha256"))
        if approved:
            if not revised_path.is_file() or sha(revised_path) != r["revised_sha256"]:
                raise ValueError(f"Changed revised ledger for {oid}")
            revised = campaign.read_obs(revised_path)
            prior_preserved(original, revised)
            worker_namespace = "followup-" + r["run_id"] + "-" + r["run_attempt"]
            problems, _ = inspect(oid, revised, root, (job, worker_namespace))
            if r["outcome"] == "complete" and problems:
                raise ValueError(f"Worker claimed completion with residual issues: {oid} {problems}")
            if problems and revised["status"] == "complete":
                revised = deferred(revised, "Repair eligibility checks incomplete: " + ", ".join(problems))
        else:
            revised = deferred(original, reason + "; " + ", ".join(m["reasons"]), root)
        if job:
            worker_namespace = "followup-" + r["run_id"] + "-" + r["run_attempt"]
            for tree in (Path("catalogs") / worker_namespace / key, Path("photometry") / key,
                         Path("repair_records") / key):
                for file in (job / tree).rglob("*"):
                    if not file.is_file():
                        continue
                    rel = file.relative_to(job)
                    if rel.parts[0] == "repair_records":
                        copy(file, revision / rel)
                    elif not approved:
                        copy(file, revision / "quarantined" / rel)
                    else:
                        target = root / rel
                        historical = revision / "before" / rel
                        # Repeating the same promotion after interruption is safe:
                        # destination may already hold the worker's exact bytes.
                        archive_replace(file, target, historical)
            copy(job / "repair_records" / f"{oid}.json", revision / "repair_records" / f"{oid}.json")
        if revised is not None:
            dest.parent.mkdir(parents=True, exist_ok=True)
            if original is None:
                save(revision / "new-ledgers" / (key + ".json"), {"worker_ledger_sha256":sha(revised_path)})
            campaign.save_obs(dest, revised)
        results[key] = {"status": revised["status"] if revised else "unattempted",
                        "validated_worker": approved, "worker_outcome": r.get("outcome", "unattempted"),
                        "errors": revised["errors"] if revised else [reason],
                        "before_sha256": m["primary_sha256"],
                        "after_sha256": sha(dest) if dest.is_file() else None}
    if sha(campaign.HERE / "vet.py") != NEW:
        copy(bundle / "vet.after.py", campaign.HERE / "vet.py")
    if campaign.identity() != plan["after"]:
        raise ValueError("Current source cannot reproduce the repaired results")
    # Small historic aggregate tables remain available alongside changed rows.
    # Unchanged catalogs, observations and photometry stay at their original paths.
    for name in ("summary.json", "leads.json", "control_trials.csv", "completeness.csv", "obsid_coverage.csv",
                 "event_receipts.csv", "followup_execution.json"):
        dest = revision / "before" / name
        if (root / name).is_file() and not dest.exists():
            copy(root / name, dest)
    new_freeze = {**freeze, "implementation": plan["after"], "git_commit": os.environ["GITHUB_SHA"],
                  "followup_revision": str(revision.relative_to(root)), "prior_campaign_sha256": sha(bundle / "campaign.before.json")}
    save(freeze_file, new_freeze)
    outcome = {"plan_sha256": plan_hash, "changed_or_attempted": results, "inventory_sha256": INVENTORY,
               "carried_complete_without_rescan": len(plan["obsids"]) - len(expected),
               "remaining_obsids": [int(k) for k, v in results.items() if v["status"] != "complete"],
               "corrected_vet_sha256": NEW, "primary_vet_sha256": OLD,
               "mode": "The historical two-line all-blank NED guard was removed; all other primary analysis code is unchanged."}
    save(final, outcome)
    print(f"{len(records)} residual workers received; {len(outcome['remaining_obsids'])} ObsIDs remain incomplete.")
    return outcome


def main():
    import subprocess
    import time
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REPOSITORY") != "dovvnloading/Open-Space"
            or os.getenv("GITHUB_REF") != "refs/heads/research/chandra-acis-fast-transients"
            or os.getenv("GITHUB_EVENT_NAME") not in {"workflow_dispatch", "push"}):
        raise RuntimeError("Use the Chandra research branch workflow")
    cli = argparse.ArgumentParser(description=__doc__)
    cmds = cli.add_subparsers(dest="command", required=True)
    for command in ("prepare", "one", "group", "finish"):
        cmd = cmds.add_parser(command)
        cmd.add_argument("evidence", type=Path)
        cmd.add_argument("bundle", type=Path)
        if command == "prepare":
            cmd.add_argument("--groups", type=int, default=32)
        if command in {"one", "group"}:
            cmd.add_argument("index", type=int)
            cmd.add_argument("out", type=Path)
            cmd.add_argument("--retained", action="append", type=Path, default=[])
        if command == "finish":
            cmd.add_argument("staged", type=Path)
            cmd.add_argument("--install-reviewed-vet", action="store_true")
    args = cli.parse_args()
    root, bundle = args.evidence.resolve(), args.bundle.resolve()
    if root == bundle or root.is_relative_to(bundle) or bundle.is_relative_to(root):
        raise ValueError("Keep the repair bundle separate from primary evidence")
    if args.command == "prepare":
        prepare(root, bundle, args.groups)
        return
    if args.command == "finish":
        staged = args.staged.resolve()
        if root == staged or root.is_relative_to(staged) or staged.is_relative_to(root):
            raise ValueError("Use a separate worker staging directory")
        finish(root, bundle, staged, args.install_reviewed_vet)
        return
    checked(root, bundle)
    plan = campaign.read(bundle / "plan.json")
    if (plan["repair_sha256"] != sha(Path(__file__)) or plan["inventory_sha256"] != INVENTORY
            or plan["catalog_reader_sha256"] != sha(Path(catalog_pack.__file__))):
        raise ValueError("Worker source differs from the original repair plan")
    out = args.out.resolve()
    if any(out == p or out.is_relative_to(p) or p.is_relative_to(out) for p in (root, bundle)):
        raise ValueError("Use a separate worker output directory")
    if not os.getenv("RETAIN_EVENTS") or not os.getenv("CIAO_PREFIX"):
        raise ValueError("Residual workers require RETAIN_EVENTS and CIAO_PREFIX")
    if args.command == "one":
        run_one(root, bundle, plan, args.index, out, [p.resolve() for p in args.retained])
    else:
        if args.index not in range(len(plan["groups"])):
            raise ValueError("Unknown residual group")
        deadline = time.monotonic() + 295 * 60
        for oid in plan["groups"][args.index]:
            remaining = int(deadline - time.monotonic())
            if remaining < 120:
                break
            cmd = [sys.executable, __file__, "one", str(root), str(bundle), str(oid), str(out)]
            for d in args.retained:
                cmd.extend(["--retained", str(d.resolve())])
            try:
                subprocess.run(cmd, check=True, timeout=remaining)
            except (subprocess.SubprocessError, OSError) as exc:
                receipt = out / "repair_records" / f"{oid}.json"
                if not receipt.is_file():
                    interrupted = old_path(out, oid)
                    if interrupted.exists():
                        copy(interrupted, out / "repair_records" / str(oid) / "quarantined.json.gz")
                    save(receipt, {"obsid": oid, "plan_sha256": sha(bundle / "plan.json"),
                                  "primary_sha256": plan["obsids"][str(oid)]["primary_sha256"],
                                  "run_id": os.environ["GITHUB_RUN_ID"], "run_attempt": os.environ["GITHUB_RUN_ATTEMPT"],
                                  "outcome": "unresolved", "revised_sha256": None, "error": str(exc)})
                print(f"ObsID {oid} worker interrupted: {exc}", flush=True)


if __name__ == "__main__":
    main()
