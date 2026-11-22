"""Branch-scoped archive orchestration. Each completed observation is an atomic checkpoint."""
import argparse
import gzip
import hashlib
import json
import importlib.metadata
import platform
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from archive import dump, fetch_events, inventory, replay

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"
N_CONTROL = 24
N_INJECT = 12


def identity():
    # Changing a search cut, vetter, or dependency starts a separate experiment.
    names = ["archive.py", "campaign.py", "search.py", "vet.py", "flux.py", "requirements.txt",
             "calibrate.sh", "report.py", "diagnostics.py"]
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in names}


def read(path):
    return json.loads(path.read_text())


def read_obs(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def save_obs(path, obj):
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    tmp = path.with_suffix(".tmp")
    tmp.write_bytes(gzip.compress(payload, compresslevel=6, mtime=0))
    tmp.replace(path)


def plan():
    folder = EVIDENCE / "inventory"
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / "campaign.json").exists():
        freeze = read(folder / "campaign.json")
        if freeze["implementation"] != identity():
            raise RuntimeError("Search revision differs from frozen evidence; archive this campaign before restarting.")
        observations = read(folder / "inventory.json")
        if hashlib.sha256((folder / "inventory.json").read_bytes()).hexdigest() != freeze["inventory_sha256"]:
            raise ValueError("Frozen inventory differs from its receipt")
    else:
        observations = inventory(folder, "2021-01-01", "2026-11-21")
        dump(folder / "campaign.json", {
            "implementation": identity(), "git_commit": os.environ["GITHUB_SHA"],
            "inventory_sha256": hashlib.sha256((folder / "inventory.json").read_bytes()).hexdigest(),
            "n_controls": N_CONTROL, "n_injections_per_product": N_INJECT,
        })
    missing = []
    for obs in observations:
        path = EVIDENCE / "observations" / f"{obs['obsid']}.json.gz"
        if not path.exists() or read_obs(path)["status"] != "complete":
            missing.append(obs)
    groups = [[] for _ in range(max(1, min(192, (len(missing) + 15) // 16)))]
    load = [0] * len(groups)
    for obs in sorted(missing, key=lambda r: (-(r["event_count"] or 0), r["obsid"])):
        group = min(range(len(groups)), key=lambda i: (load[i], i))
        groups[group].append(obs["obsid"])
        load[group] += max(1000, obs["event_count"] or 0)
    dump(folder / "shards.json", groups)
    with open(os.environ["GITHUB_OUTPUT"], "a") as target:
        print("matrix=" + json.dumps(list(range(len(groups)))), file=target)


def response(trigger, evt, directory):
    from flux import peak_flux
    prefix = os.environ.get("CIAO_PREFIX")
    if prefix and not (Path(prefix) / "bin/eff2evt").is_file():
        subprocess.run(["bash", str(HERE / "calibrate.sh")], check=True, timeout=1800)
    return peak_flux(evt, trigger, directory)


def annotate(trigger, directory, event_file, stream, index, cache, output):
    from vet import vet_trigger
    trigger["trigger_key"] = f"{stream}/{index}"
    previous_vet = trigger.get("vet", {})
    if not previous_vet or previous_vet["status"] == "catalog_incomplete":
        trigger["catalog_cache"] = str(cache.relative_to(output))
        try:
            if trigger.get("ra_deg") is None or trigger.get("dec_deg") is None:
                raise ValueError("Missing WCS localization")
            trigger["vet"] = vet_trigger(trigger, cache)
        except Exception as exc:
            trigger["vet"] = {"status": "catalog_incomplete", "reason": f"{type(exc).__name__}: {exc}"}
    # Photometry characterizes real photons after all selection decisions.
    # A scrambled time assignment has no physical peak-flux measurement.
    if stream == "real" and trigger["instrument_pass"] and trigger["vet"]["status"] not in (
        "foreground_veto", "known_transient"
    ):
        flux_dir = directory / stream / str(index)
        flux_dir.mkdir(parents=True, exist_ok=True)
        if "diagnostics" not in trigger or "error" in trigger["diagnostics"]:
            from diagnostics import save_diagnostics
            try:
                trigger["diagnostics"] = save_diagnostics(event_file, trigger, flux_dir)
            except Exception as exc:
                trigger["diagnostics"] = {"error": f"{type(exc).__name__}: {exc}"}
        if "peak_flux_erg_cm2_s" not in trigger.get("flux", {}):
            try:
                trigger["flux"] = response(trigger, event_file, flux_dir)
            except Exception as exc:
                trigger["flux"] = {"error": f"{type(exc).__name__}: {exc}"}


def process(obsid, output, event_dir, cache):
    from search import search_file
    ledger = {"obsid": obsid, "status": "error", "products": [], "errors": [],
              "catalog_cache": str(cache.relative_to(output))}
    path = output / "observations" / f"{obsid}.json.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    old = EVIDENCE / "observations" / path.name
    prior = read_obs(old) if old.exists() else {}
    old_products = {r["input"]["sha256"]: r for r in prior.get("products", [])}
    try:
        if pinned := prior.get("acquired"):
            with tempfile.TemporaryDirectory() as temp:
                manifest = Path(temp) / "products.json"
                dump(manifest, pinned)
                replay(manifest, event_dir)
            products = [dict(p, path=str(event_dir / p["filename"])) for p in pinned]
        else:
            products = fetch_events(obsid, event_dir)
        ledger["acquired"] = [{k: v for k, v in p.items() if k != "path"} for p in products]
        ledger["products"] = [old_products.get(p["sha256"], {"input": receipt})
                              for p, receipt in zip(products, ledger["acquired"], strict=True)]
        save_obs(path, ledger)
        for p, record in zip(products, ledger["products"], strict=True):
            try:
                record["complete"] = False
                seed = int(hashlib.sha256(("chandra-acis-fxt-v1:" + p["sha256"]).encode()).hexdigest()[:16], 16) >> 1
                event_file = Path(p["path"])
                result = record.get("search")
                if result is None:
                    result = search_file(event_file, seed=seed, n_controls=N_CONTROL, n_injections=N_INJECT)
                    record["search"] = result
                    save_obs(path, ledger)
                directory = output / "photometry" / str(obsid) / p["filename"]
                streams = [("real", result["triggers"])] + [
                    (f"control-{i}", trial["triggers"]) for i, trial in enumerate(result["controls"])
                ]
                if len(result["controls"]) != N_CONTROL:
                    raise ValueError("Incomplete time-scramble ensemble")
                keep = False
                for stream, triggers in streams:
                    for i, trigger in enumerate(triggers):
                        trigger["obsid"] = obsid
                        annotate(trigger, directory, event_file, stream, i, cache, output)
                        if stream == "real" and trigger["instrument_pass"] and trigger["vet"]["status"] not in (
                            "foreground_veto", "known_transient"
                        ):
                            keep = True
                if keep:
                    depot = Path(os.environ["RETAIN_EVENTS"]) / str(obsid)
                    depot.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(event_file, depot / p["filename"])
                record["complete"] = True
            except Exception as exc:
                ledger["errors"].append(f"{p['filename']} {type(exc).__name__}: {exc}")
            save_obs(path, ledger)
        if ledger["products"] and all(r.get("complete") for r in ledger["products"]):
            # A worker can finish scanning even if an external service failed.
            # A later dispatch retries the frozen input in those observations.
            ledger["status"] = "complete"
            for record in ledger["products"]:
                scan = record["search"]
                for trial in [scan, *scan["controls"]]:
                    for event in trial["triggers"]:
                        if (event["vet"]["status"] == "catalog_incomplete" or event.get("flux", {}).get("error")
                                or event.get("diagnostics", {}).get("error")):
                            ledger["status"] = "followup_incomplete"
    except Exception as exc:
        if prior and not ledger["products"]:
            ledger = prior
        ledger["errors"].append(f"{type(exc).__name__}: {exc}")
    if (event_dir / "acquisition.json").is_file():
        ledger["archive_attempts"] = read(event_dir / "acquisition.json")
    elif "archive_attempts" in prior:
        ledger["archive_attempts"] = prior["archive_attempts"]
    save_obs(path, ledger)


def shard(group, output):
    deadline = time.monotonic() + 295 * 60
    groups = read(EVIDENCE / "inventory/shards.json")
    output.mkdir(parents=True, exist_ok=True)
    dump(output / f"environment-{group}.json", {
        "git_commit": os.environ["GITHUB_SHA"], "run_id": os.environ["GITHUB_RUN_ID"],
        "python": sys.version, "platform": platform.platform(), "runner_image": os.getenv("ImageVersion"),
        "packages": {p.metadata["Name"]: p.version for p in importlib.metadata.distributions() if p.metadata["Name"]},
    })
    cache = output / "catalogs" / f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}-{group}"
    cache.mkdir(parents=True, exist_ok=True)
    print(f"Shard {group} has {len(groups[group])} assigned observations", flush=True)
    for obsid in groups[group]:
        remaining = int(deadline - time.monotonic())
        if remaining < 120:
            break
        print(f"ObsID {obsid} started", flush=True)
        with tempfile.TemporaryDirectory() as tmp:
            try:
                run = subprocess.run([sys.executable, __file__, "one", str(obsid), str(output),
                                      tmp, str(cache)], timeout=remaining)
                if run.returncode:
                    raise RuntimeError(f"worker exited {run.returncode}")
            except (subprocess.TimeoutExpired, RuntimeError) as exc:
                path = output / "observations" / f"{obsid}.json.gz"
                path.parent.mkdir(parents=True, exist_ok=True)
                row = read_obs(path) if path.exists() else {"obsid": obsid, "products": [], "errors": []}
                row["status"] = "error"
                row["errors"].append(str(exc))
                save_obs(path, row)
        current = read_obs(output / "observations" / f"{obsid}.json.gz")
        print(f"ObsID {obsid} {current['status']}: {str(current['errors'])[:500]}", flush=True)
    dump(output / f"shard-{group}.json", {"group": group, "assigned": groups[group],
         "completed_worker": True, "commit": os.environ["GITHUB_SHA"]})


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    cmd = args.add_subparsers(dest="mode", required=True)
    cmd.add_parser("plan")
    sh = cmd.add_parser("shard")
    sh.add_argument("group", type=int)
    sh.add_argument("output", type=Path)
    one = cmd.add_parser("one")
    one.add_argument("obsid", type=int)
    for field in ("output", "event_dir", "cache"):
        one.add_argument(field, type=Path)
    a = args.parse_args()
    if a.mode == "plan":
        plan()
    elif a.mode == "one":
        process(a.obsid, a.output, a.event_dir, a.cache)
    else:
        shard(a.group, a.output)
