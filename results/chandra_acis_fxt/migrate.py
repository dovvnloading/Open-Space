"""Reclassify a pinned campaign after the CHANMASTER-independent catalog fix.

Run before campaign.py plan. STAGING holds separately extracted chandra-evidence-N
artifacts from run 44981778618. Supply after-hashes literally from reviewed code
in the same git revision's workflow, including any newly tracked policy files.
The archived run and its original catalog decisions remain available in history.
"""
import argparse
import copy
import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile

from archive import dump
from campaign import EVIDENCE, N_CONTROL, N_INJECT, identity, read, read_obs, save_obs

RUN = "44981778618"
COMMIT = "75870e846a6139ad52d756b4409db768ecd901ad"
INVENTORY = "172b673a39b7c2bca8b82294563a62a7cb5807fdcc393a8da27d31e56f3c518e"
BEFORE = {
    "archive.py": "4623fb84527aa4952a2beab6d2405b6b6e0f5ddd0ffdfd69e1c96bfeacca0a7a",
    "calibrate.sh": "fde33cbc008ab369e4539c3c6abd2fba3ab686a60031f32239fd3ba70a113ec8",
    "campaign.py": "654b5d6aba7592f4c4b385b75a275b34d48c4e1d428ef555b55939ea818f94b5",
    "diagnostics.py": "c0844ef26db821503267349ded5dcd6a872e84b2ab6caaeb6769c239e578829e",
    "flux.py": "e4a1906e24c987dcc5f2d8320d9347a69ae96c3da2022283d3fcdd67def0eaad",
    "report.py": "2b8adc1133502baa7b7fcbb9f16b06bb0953df1427a3143cb5dac0664021688d",
    "requirements.txt": "f3ee4d1c5cbf58c313adf4c3926484e4e3f5736834225d3dbf5976f518660790",
    "search.py": "17bdbb0ec1fa5c86bf1a3394b8df29d3b245e436dc0bb8ce2376a5c9e33a4f99",
    "vet.py": "b1ba7998a4ededd38c9bcbf2d803ecb579eb9eedc09f2b7c02f4d582dc3ee503",
}


def file_sha(path):
    with path.open("rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def pin_transition(literals):
    reviewed = {}
    allowed = {"vet.py", "campaign.py", "report.py", "alerts.py", "visits.py", "exceptions.py"}
    for spec in literals:
        name, sep, value = spec.partition("=")
        if not sep or name not in allowed or name in reviewed or not re.fullmatch("[a-f0-9]{64}", value):
            raise ValueError("Invalid reviewed implementation entry: " + spec)
        reviewed[name] = value
    if reviewed.get("vet.py", BEFORE["vet.py"]) == BEFORE["vet.py"]:
        raise ValueError("The correction requires an explicitly reviewed new vet.py hash")
    after = {**BEFORE, **reviewed}
    if identity() != after:
        raise ValueError(f"Current implementation differs from the reviewed transition: {identity()}")
    return after


def preserve_characterization(event, root, out, receipt, stream, index):
    # Retain only input-bound products that replay can read, independently of
    # the new catalog disposition. Error objects must be recomputed by workers.
    original_flux = event.pop("flux", {})
    original_diagnostics = event.pop("diagnostics", {})
    if stream != "real":
        return
    rel = Path("photometry") / str(receipt["obsid"]) / receipt["filename"] / "real" / str(index)
    directory = root / rel
    files = set()
    if (original_flux.get("event_sha256") == receipt["sha256"]
            and "peak_flux_erg_cm2_s" in original_flux
            and (directory / "flux.json").is_file()
            and read(directory / "flux.json") == original_flux
            and (directory / "peak_weights.fits.gz").is_file()
            and file_sha(directory / "peak_weights.fits.gz") == original_flux.get("calibrated_fits_sha256")
            and (directory / "peak_selection.json").is_file()):
        event["flux"] = original_flux
        files.update(name for name in ("flux.json", "peak_weights.fits.gz", "peak_selection.json",
                                       "commands.json", "environment.json", "ciaover.log",
                                       "eff2evt.log", "dmcopy.log") if (directory / name).is_file())
    if original_diagnostics.get("input_sha256") == receipt["sha256"]:
        artifacts = original_diagnostics.get("files", {})
        if (artifacts and (directory / "diagnostics.json").is_file()
                and read(directory / "diagnostics.json") == original_diagnostics
                and all(Path(n).name == n and (directory / n).is_file()
                        and file_sha(directory / n) == digest for n, digest in artifacts.items())):
            event["diagnostics"] = original_diagnostics
            files.update([*artifacts, "diagnostics.json"])
    if files:
        (out / rel).mkdir(parents=True, exist_ok=True)
        for name in files:
            shutil.copyfile(directory / name, out / rel / name)


def verify_scan(result, receipt, oid):
    expected_seed = int(hashlib.sha256(("chandra-acis-fxt-v1:" + receipt["sha256"]).encode()).hexdigest()[:16], 16) >> 1
    if (result["seed"] != expected_seed or int(result["obsid"]) != oid
            or [c["replica"] for c in result["controls"]] != list(range(N_CONTROL))
            or [j["trial"] for j in result["injections"]] != list(range(N_INJECT))):
        raise ValueError(f"ObsID {oid}: scan identity, controls or injections differ from frozen campaign")


def transition_ledger(path, historical_root, output, inventory):
    row = read_obs(path)
    oid = int(path.name.split(".")[0])
    if oid not in inventory or row["obsid"] != oid:
        raise ValueError(f"Out-of-census ledger or ObsID mismatch: {path.name}")
    original_status = row["status"]
    receipts = {}
    for source in row.get("acquired", []):
        filename = source["filename"]
        if (not re.fullmatch(r"acisf[A-Za-z0-9_]+_evt2\.fits(?:\.gz)?", filename) or
                filename in receipts or int(source["obsid"]) != oid or source["bytes"] <= 0 or
                not re.fullmatch("[0-9a-f]{64}", source["sha256"])):
            raise ValueError(f"Invalid or duplicated archived input receipt, ObsID {oid}")
        receipts[filename] = source
    scanned = set()
    triggers = 0
    for product in row["products"]:
        source = product["input"]
        filename = source["filename"]
        if source != receipts.get(filename) or filename in scanned:
            raise ValueError(f"ObsID {oid}: product does not match acquired bytes")
        result = product.get("search")
        if result is None:
            product["complete"] = False
            continue
        verify_scan(result, source, oid)
        scanned.add(filename)
        for stream, events in [("real", result["triggers"]), *[
                (f"control-{c['replica']}", c["triggers"]) for c in result["controls"]]]:
            for index, event in enumerate(events):
                triggers += 1
                event.pop("vet", None)
                event.pop("catalog_cache", None)
                preserve_characterization(event, historical_root, output, source, stream, index)
        if triggers:
            product["complete"] = False
    can_promote = (original_status == "complete" and receipts and scanned == receipts.keys()
                   and all(r.get("complete") for r in row["products"]) and triggers == 0)
    row["status"] = "complete" if can_promote else "pending_corrected_followup"
    # Missing catalog annotations are never presented to the reducer as finished.
    if not can_promote:
        for product in row["products"]:
            product["complete"] = False
    save_obs(output / "observations" / path.name, row)
    return {"obsid": oid, "prior_status": original_status, "status": row["status"], "scanned_products": len(scanned),
            "triggers_to_revet": triggers, "historical_ledger": str(path.relative_to(output)),
            "prior_sha256": file_sha(path)}


def migrate(staging, reviewed):
    after = pin_transition(reviewed)
    freeze_path = EVIDENCE / "inventory/campaign.json"
    freeze = read(freeze_path)
    if file_sha(EVIDENCE / "inventory/inventory.json") != INVENTORY or freeze["inventory_sha256"] != INVENTORY:
        raise ValueError("Frozen population differs from the reviewed source campaign")
    if (N_CONTROL, N_INJECT, freeze["n_controls"], freeze["n_injections_per_product"]) != (24, 12, 24, 12):
        raise ValueError("Control or injection design changed")
    if freeze["implementation"] == after:
        record = read(EVIDENCE / "history" / RUN / "migration.json")
        if record["before"] != BEFORE or record["after"] != after:
            raise ValueError("Current receipt lacks the authorized campaign history")
        print("The reviewed catalog migration is already recorded.")
        return
    if freeze["implementation"] != BEFORE or freeze["git_commit"] != COMMIT:
        raise ValueError("Only the original 44981778618 implementation can be migrated")
    population = {r["obsid"] for r in read(EVIDENCE / "inventory/inventory.json")}
    assignments = read(EVIDENCE / "inventory/shards.json")
    artifacts = sorted(staging.glob("chandra-evidence-*"))
    # Assemble and validate a complete replacement offline before touching the
    # checked-out evidence. A failed validation leaves the prior tree intact.
    with tempfile.TemporaryDirectory(dir=EVIDENCE.parent) as work:
        new = Path(work) / "evidence"
        shutil.copytree(EVIDENCE, new)
        history = new / "history" / RUN
        history.mkdir(parents=True)
        branch = history / "checkpoint"
        branch.mkdir()
        for name in ("observations", "catalogs", "photometry"):
            path = new / name
            if path.exists():
                path.rename(branch / name)
        (new / "observations").mkdir()
        for name in ("inventory/campaign.json", "inventory/shards.json", "pilot.json"):
            archived = branch / name
            archived.parent.mkdir(parents=True, exist_ok=True)
            if (new / name).is_file():
                shutil.copyfile(new / name, archived)
        origins = [(branch, population)]
        for artifact in artifacts:
            match = re.fullmatch(r"chandra-evidence-(\d+)", artifact.name)
            if not match or not artifact.is_dir() or any(p.is_symlink() for p in artifact.rglob("*")):
                raise ValueError(f"Invalid source artifact {artifact.name}")
            group = int(match[1])
            environment = read(artifact / f"environment-{group}.json")
            if (group >= len(assignments) or str(environment["run_id"]) != RUN
                    or environment["git_commit"] != COMMIT):
                raise ValueError(f"Wrong run, commit or assignment for {artifact.name}")
            if (artifact / f"shard-{group}.json").exists():
                shard = read(artifact / f"shard-{group}.json")
                if shard["commit"] != COMMIT or shard["assigned"] != assignments[group]:
                    raise ValueError(f"ObsID assignment disagrees for {artifact.name}")
            archived = history / artifact.name
            shutil.copytree(artifact, archived)
            origins.append((archived, set(assignments[group])))
        index = []
        prior_files = {}
        for origin, assigned in origins:
            for path in sorted((origin / "observations").glob("*.json.gz")):
                oid = int(path.name.split(".")[0])
                if oid not in assigned:
                    raise ValueError(f"Unassigned ObsID {oid} in {origin.name}")
                if oid in prior_files:
                    if file_sha(path) != file_sha(prior_files[oid]):
                        raise ValueError(f"Conflicting source ledgers for ObsID {oid}")
                    continue
                prior_files[oid] = path
                index.append(transition_ledger(path, origin, new, population))
        receipt = {"source_run": RUN, "source_git_commit": COMMIT,
                   "migration_run": os.environ["GITHUB_RUN_ID"], "migration_git_commit": os.environ["GITHUB_SHA"],
                   "before": BEFORE, "after": after, "frozen_inventory_sha256": INVENTORY,
                   "artifacts": [a.name for a in artifacts], "ledgers": index,
                   "promoted_no_trigger": sum(r["status"] == "complete" for r in index),
                   "triggers_to_revet": sum(r["triggers_to_revet"] for r in index)}
        dump(history / "migration.json", receipt)
        revised = copy.deepcopy(freeze)
        revised["implementation_history"] = [str((history / "migration.json").relative_to(new))]
        revised["implementation"] = after
        revised["git_commit"] = os.environ["GITHUB_SHA"]
        revised["original_git_commit"] = COMMIT
        dump(new / "inventory/campaign.json", revised)
        backup = Path(work) / "before"
        EVIDENCE.rename(backup)
        try:
            new.rename(EVIDENCE)
        except BaseException:
            backup.rename(EVIDENCE)
            raise
    print(f"Migrated {len(index)} ledgers from {len(artifacts)} artifacts and the branch checkpoint. "
          f"{receipt['promoted_no_trigger']} require no catalog operations. "
          f"{receipt['triggers_to_revet']} trigger classifications are queued for replacement.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("staging", type=Path)
    parser.add_argument("--reviewed", action="append", required=True)
    args = parser.parse_args()
    migrate(args.staging, args.reviewed)
