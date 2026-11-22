"""Repeat a completed primary group's missing upload, keeping the two producers distinct."""
import argparse
from contextlib import contextmanager
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import sys
import tempfile
from urllib.parse import unquote, urlsplit
from unittest.mock import patch

import archive
import campaign
import followup
import ingest_primary as ingest

need, sha, save = ingest.require, ingest.sha, ingest.save
HISTORY = Path("upload_recovery") / f"{ingest.RUN}-{ingest.ATTEMPT}"
PART = 20 * 1024**2
# A retained artifact identity witnessed before this repair was written.
WITNESSED = {72: (137101488955, 14747907629, 4863648,
                 "sha256:1d73a35f0639700215c53232855715e18fd267d9956d8486c207a593733af2cb")}
FIELDS = ("id", "run_attempt", "head_sha", "head_branch", "path", "status", "conclusion", "html_url")


def run_info(run_id, attempt, commit):
    row = ingest.api(f"runs/{run_id}")
    need((row["id"], row["run_attempt"], row["head_sha"], row["head_branch"], row["repository"]["full_name"]) ==
         (run_id, attempt, commit, ingest.BRANCH, ingest.REPO) and row["path"].split("@")[0] == ingest.FLOW,
         "Unexpected Actions producer")
    return {k: row[k] for k in FIELDS} | {"repository_id": row["repository"]["id"]}


def get_item(name, source):
    response = ingest.api(f"runs/{ingest.RUN}/artifacts?name={name}&per_page=100")
    items = response["artifacts"]
    need(response["total_count"] == len(items) <= 1, "Nonunique artifact lookup")
    if not items:
        return None
    entry = items[0]
    parent = entry["workflow_run"]
    need(entry["name"] == name and parent["id"] == ingest.RUN and parent["head_sha"] == ingest.CODE and
         parent["head_branch"] == ingest.BRANCH and parent["repository_id"] == source["repository_id"] ==
         parent["head_repository_id"] and re.fullmatch("sha256:[0-9a-f]{64}", entry.get("digest") or ""),
         "Foreign artifact or missing provider SHA-256")
    return {k: entry[k] for k in ("id", "name", "size_in_bytes", "digest", "created_at", "expired", "expires_at")}


def failed_job(job, group, job_id, count, log):
    need(job["id"] == job_id and job.get("run_id") == ingest.RUN and job["run_attempt"] == ingest.ATTEMPT and
         job["head_sha"] == ingest.CODE and job["status"] == "completed" and job["conclusion"] == "failure" and
         re.fullmatch(rf"observation \((?:group: )?{group}\)", job["name"]), "Wrong failed original worker")
    steps = sorted(job["steps"], key=lambda s: s["number"])
    search = [i for i, s in enumerate(steps) if s["name"] == "Search assigned observations and record all outcomes"]
    need(len(search) == 1 and search[0] + 1 < len(steps), "Original worker science/upload boundary missing")
    i = search[0]
    failed = steps[i+1]
    need(steps[i]["conclusion"] == "success" and failed["conclusion"] == "failure" and
         "upload-artifact" in failed["name"] and f"Shard {group} has {count} assigned observations" in log and
         "CreateArtifact" in log and re.search(r"(?:timed? ?out|timeout|ECONNRESET|503|502|service unavailable)", log, re.I),
         "No authenticated successful search and failed evidence upload")
    status_lines = re.findall(r"\bObsID (\d+) (complete|followup_incomplete|error):", log)
    return {"id": job_id, "url": job["html_url"], "name": job["name"], "conclusion": job["conclusion"],
            "original_worker_status_lines": status_lines,
            "search_step": steps[i], "failed_upload_step": failed}


def save_blob(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        need(path.read_bytes() == content, "Different bytes at " + path.name)
    else:
        tmp = path.with_name(path.name + ".tmp")
        tmp.write_bytes(content)
        tmp.replace(path)


def stream_key(name):
    return re.sub(r"N\d+", "N#", name.removesuffix(".gz"), count=1, flags=re.I)


def source_products(stage, temp, artifact, ids):
    if not artifact:
        return {}, []
    need(not artifact["expired"], "Original retained EVT2 artifact expired")
    zipped = temp / "original.zip"
    ingest.api(f"artifacts/{artifact['id']}/zip", zipped)
    need(sha(zipped) == artifact["digest"][7:] and zipped.stat().st_size == artifact["size_in_bytes"],
         "Original retained EVT2 ZIP differs from GitHub receipt")
    parts = []
    with zipped.open("rb") as src:
        while block := src.read(PART):
            part = stage / "original" / "event_zip" / f"part-{len(parts):04d}.bin"
            save_blob(part, block)
            parts.append({"path": str(part.relative_to(stage)), "sha256": sha(part), "bytes": len(block)})
    unpacked = temp / "retained"
    unpacked.mkdir()
    ingest.unzip(zipped, unpacked)
    products = {}
    for name, file in sorted(ingest.logical_files(unpacked).items()):
        p = Path(name)
        need(len(p.parts) == 2 and p.parts[0].isdigit() and int(p.parts[0]) in ids and
             str(int(p.parts[0])) == p.parts[0] and
             archive.product_names(f'<a href="{p.name}">event</a>', int(p.parts[0])) == [p.name],
             "Foreign/malformed product in original retained artifact: " + name)
        oid = int(p.parts[0])
        item = {"member": name, "filename": p.name, "obsid": oid, "bytes": file.stat().st_size,
                "sha256": sha(file), **archive.inspect_product(file, oid)}
        k = (oid, stream_key(p.name))
        need(k not in products, "Duplicate original retained EVT2 stream")
        products[k] = (file, item)
    need(products, "Published retained artifact contains no EVT2 product")
    return products, parts


@contextmanager
def exact_retained(oid, products, artifact):
    selected = {key: v for (i, key), v in products.items() if i == oid}
    names, download = archive.product_names, archive.download

    def choose(html, obsid):
        need(obsid == oid, "Recovery selector called for another ObsID")
        streams = {stream_key(n): n for n in names(html, obsid)}
        need(set(selected) <= streams.keys(), "Original retained stream absent from new archive inventory")
        return sorted(selected[k][1]["filename"] if k in selected else n for k, n in streams.items())

    def fetch(url, dest, data=None):
        name = unquote(urlsplit(url).path.rsplit("/", 1)[-1])
        key = stream_key(name)
        if key not in selected:
            return download(url, dest, data)
        archive.check_archive_url(url)
        src, item = selected[key]
        need(item["filename"] == name and data is None and
             f"/{oid}/primary/" in urlsplit(url).path, "Pinned retained stream or archive path differs")
        shutil.copyfile(src, dest)
        need(sha(dest) == item["sha256"], "Restored EVT2 changed")
        return {"url": url, "resolved_url": None, "sha256": item["sha256"], "bytes": item["bytes"],
                "source_kind": "original_retained_actions_member", "artifact_id": artifact["id"],
                "artifact_digest": artifact["digest"], "source_member": item["member"]}

    with patch.object(archive, "product_names", choose), patch.object(archive, "download", fetch):
        yield


def witness_match(oid, row, products):
    actual = {(e["filename"], e["sha256"], e["fits_sha256"], e["bytes"]) for e in row.get("acquired", [])}
    needed = {(r["filename"], r["sha256"], r["fits_sha256"], r["bytes"])
              for (i, _), (_, r) in products.items() if i == oid}
    return needed <= actual


def attest_environment():
    return {"python": sys.version, "platform": platform.platform(), "image": os.getenv("ImageVersion"),
            "packages": {p.metadata["Name"]: p.version for p in importlib.metadata.distributions() if p.metadata["Name"]}}


def filemap(folder):
    need(not any(p.is_symlink() for p in folder.rglob("*")), "Unexpected symbolic link in recovery")
    return {name: {"bytes": p.stat().st_size, "sha256": sha(p)} for name, p in sorted(ingest.logical_files(folder).items())
            if name != "seal.json"}


def verify_stage(stage):
    seal = campaign.read(stage / "seal.json")
    need(seal.get("schema") == 1 and seal["files"] == filemap(stage), "Recovery file manifest mismatch")
    return seal


def recover(root, stage, group, job_id):
    need(0 <= group < 192 and int(os.environ["GITHUB_RUN_ID"]) != ingest.RUN, "Require a distinct repeat run")
    frozen = ingest.Frozen(root)
    ids = frozen.groups[group]
    stage.mkdir(parents=True, exist_ok=True)
    source = run_info(ingest.RUN, ingest.ATTEMPT, ingest.CODE)
    repeat = run_info(int(os.environ["GITHUB_RUN_ID"]), int(os.environ["GITHUB_RUN_ATTEMPT"]), os.environ["GITHUB_SHA"])
    checked = ingest.git("rev-parse", "HEAD").decode().strip()
    need(repeat["head_sha"] == checked, "Recovery checkout must match its own triggering SHA")
    need(get_item(f"chandra-evidence-{group}", source) is None, "The original evidence upload exists: use primary ingest")
    artifact = get_item(f"chandra-event-inputs-{group}", source)
    if group in WITNESSED:
        need(artifact and (job_id, artifact["id"], artifact["size_in_bytes"], artifact["digest"]) == WITNESSED[group],
             "Unexpected witnessed original job or EVT2 artifact")
    stable = {"group": group, "assigned": ids, "source_run": ingest.RUN, "source_attempt": ingest.ATTEMPT,
              "source_code": ingest.CODE, "assignment_commit": ingest.PLAN, "assignment": frozen.sources["inventory/shards.json"],
              "job_id": job_id, "repeat_run": repeat["id"], "repeat_attempt": repeat["run_attempt"],
              "repeat_head": checked, "science": frozen.science, "adapter_sha256": sha(Path(__file__))}
    if (stage / "seal.json").is_file():
        need(verify_stage(stage)["identity"] == stable, "Different recovery filled this stage")
        return
    identity = stage / "identity.json"
    if identity.exists():
        need(campaign.read(identity) == stable, "Conflicting partial recovery identity")
    else:
        need(not any(stage.iterdir()), "Refuse unowned recovery destination")
        save(identity, stable)
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as tmp:
        temp = Path(tmp)
        job = ingest.api(f"jobs/{job_id}")
        original_log = temp / "job.log"
        ingest.api(f"jobs/{job_id}/logs", original_log)
        raw = original_log.read_bytes()
        failure = failed_job(job, group, job_id, len(ids), raw.decode("utf-8", errors="replace"))
        save_blob(stage / "original/job.log.gz", gzip.compress(raw, mtime=0))
        save(stage / "original/failure.json", {"job": failure, "run": source, "index_checked_missing": f"chandra-evidence-{group}",
             "retained_artifact": artifact, "log_sha256": hashlib.sha256(raw).hexdigest()})
        products, parts = source_products(stage, temp, artifact, ids)
        save(stage / "original/retained.json", {"provider": artifact, "parts": parts,
             "products": [record for _, record in products.values()]})
        need({**campaign.identity(), "vet.py": followup.OLD} == frozen.science and frozen.science["vet.py"] == followup.OLD,
             "Unreviewed change to original search")
        followup.source_pair(temp / "code")
        followup.load_vet(temp / "code/vet.before.py")
        save(stage / "repeat.json", {"producer": repeat, "checkout": checked, "environment": attest_environment()})
        prior = temp / "prior"
        data = stage / "science"
        cache = data / "catalogs" / f"upload-recovery-{repeat['id']}-{repeat['run_attempt']}-{group}"
        cache.mkdir(parents=True, exist_ok=True)
        outcomes = {}
        for oid in ids:
            rel = f"observations/{oid}.json.gz"
            if rel in frozen.blobs:
                save_blob(prior / rel, frozen.baseline(rel))
            checkpoint = stage / "receipts" / f"{oid}.json"
            if checkpoint.is_file():
                receipt = campaign.read(checkpoint)
                need(sha(data / rel) == receipt["ledger_sha256"], "A finished recovery ledger changed")
            else:
                print(f"Recovery group {group} ObsID {oid} started", flush=True)
                with tempfile.TemporaryDirectory(dir=temp) as events:
                    event_dir = Path(events)
                    with patch.object(campaign, "EVIDENCE", prior), exact_retained(oid, products, artifact):
                        campaign.process(oid, data, event_dir, cache)
                    for src in [*event_dir.glob("primary_*.html"), event_dir / "acquisition.json"]:
                        if src.is_file():
                            ingest.install(src, stage / "acquisition" / str(oid) / src.name)
                row = ingest.ledger(data / rel, oid)
                witnessed = witness_match(oid, row, products)
                # A failed intake remains visible to the ordinary repair plan; never count a complete mismatch.
                need(witnessed or row["status"] not in {"complete", "followup_incomplete"},
                     "Recovery omitted or changed a surviving original EVT2 product")
                receipt = {"obsid": oid, "status": row["status"], "ledger_sha256": sha(data / rel),
                           "original_retained_products_accounted": witnessed,
                           "checkpoint_reuse": frozen.blobs.get(rel),
                           "products": [{"filename": e["filename"], "sha256": e["sha256"],
                                        "fits_sha256": e["fits_sha256"], "url": e["url"], "bytes": e["bytes"]}
                                        for e in row.get("acquired", [])]}
                save(checkpoint, receipt)
            outcomes[str(oid)] = receipt
        save(stage / "seal.json", {"schema": 1, "identity": stable, "observations": outcomes, "files": filemap(stage),
             "scope": "Independent repeat. Original archive bytes are attested only for listed retained and checkpoint products. "
                      "Remaining archive and catalog receipts identify the repeat. The original producer upload is missing."})
        verify_stage(stage)


def import_repeat(root, group, rid, attempt, commit):
    need(rid != ingest.RUN and re.fullmatch("[0-9a-f]{40}", commit) and 0 <= group < 192,
         "Invalid repeat producer identity")
    repeat = run_info(rid, attempt, commit)
    jobs = ingest.pages(f"runs/{rid}/attempts/{attempt}/jobs", "jobs")
    producers = [j for j in jobs if j["name"] == "recover-upload"]
    need(len(producers) == 1, "Unexpected repeat producer job")
    producer = producers[0]
    need(producer["run_id"] == rid and producer["head_sha"] == commit and producer["run_attempt"] == attempt
         and producer["status"] == "completed" and producer["conclusion"] == "success",
         "The repeat producer job must finish successfully before import")
    name = f"chandra-missing-upload-{group}"
    result = ingest.api(f"runs/{rid}/artifacts?name={name}&per_page=100")
    need(result["total_count"] == len(result["artifacts"]) == 1, "Expected one archived repeat artifact")
    item = result["artifacts"][0]
    parent = item["workflow_run"]
    need(not item["expired"] and item["name"] == name and parent["id"] == rid and
         parent["repository_id"] == repeat["repository_id"] == parent["head_repository_id"] and
         parent["head_sha"] == commit and parent["head_branch"] == ingest.BRANCH and
         re.fullmatch("sha256:[0-9a-f]{64}", item.get("digest") or ""), "Wrong repeat artifact provenance")
    stage = root / HISTORY / str(group) / f"{rid}-{attempt}" / "stage"
    transport = stage.parent / "transport.json"
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as directory:
        temp = Path(directory)
        target = temp / "repeat.zip"
        ingest.api(f"artifacts/{item['id']}/zip", target)
        need(sha(target) == item["digest"][7:] and item["size_in_bytes"] == target.stat().st_size,
             "Repeat artifact differs from provider bytes")
        unpacked = temp / "unpacked"
        unpacked.mkdir()
        ingest.unzip(target, unpacked)
        seal = verify_stage(unpacked)
        identity = seal["identity"]
        need((identity["repeat_run"], identity["repeat_attempt"], identity["repeat_head"], identity["group"]) ==
             (rid, attempt, commit, group) and identity["source_run"] == ingest.RUN and identity["source_code"] == ingest.CODE
             and identity["source_attempt"] == ingest.ATTEMPT and identity["assignment_commit"] == ingest.PLAN,
             "Repeat artifact metadata does not match the two producers")
        for rel, src in ingest.logical_files(unpacked).items():
            dst = stage / rel
            if dst.exists():
                need(sha(dst) == sha(src), "A deposited repeat artifact member changed")
            else:
                ingest.install(src, dst)
        verify_stage(stage)
        stamp = {"artifact": {k: item[k] for k in ("id", "name", "digest", "size_in_bytes", "created_at")},
                 "producer_job": {k: producer[k] for k in ("id","name","html_url","started_at","completed_at","conclusion")},
                 "repeat_run": rid, "repeat_attempt": attempt, "repeat_head": commit,
                 "group": group, "seal_sha256": sha(stage / "seal.json")}
        if transport.exists():
            need(campaign.read(transport) == stamp, "A different transport is already deposited")
        else:
            save(transport, stamp)
    print(stage, flush=True)
    return stage


def adoption_state(root, group, job, ids, source, coverage):
    """Only a terminal missing original and a separately sealed original held-out evaluation may advance."""
    import heldout_original
    name = f"chandra-evidence-{group}"
    row = coverage["groups"][str(group)]
    need(source["run"]["status"] == "completed" and row["artifact_status"] == "missing_artifact" and
         row["artifact_id"] is None and name not in source["artifacts"] and not ingest.receipt_file(root, str(group)).exists()
         and source["jobs"][str(group)]["id"] == row["job"]["id"] == job and row["job"]["conclusion"] == "failure"
         and set(ids) == set(row["assigned_without_verified_worker_ledger"]), "Original missing-upload status changed")
    seal = root / "held_out" / heldout_original.RECEIPT
    need(seal.is_file() and campaign.read(seal)["source_run"] == ingest.RUN,
         "Seal the original held-out evidence before adopting an independent repeat")
    return seal


def adopt(root, stage):
    import heldout_original
    seal = verify_stage(stage)
    i = seal["identity"]
    g, job = i["group"], i["job_id"]
    rel = HISTORY / str(g) / f"{i['repeat_run']}-{i['repeat_attempt']}" / "stage"
    need(stage == (root / rel).resolve(), "Use the GitHub-digest-verified imported stage")
    transport = campaign.read(stage.parent / "transport.json")
    need((transport["group"], transport["repeat_run"], transport["repeat_attempt"], transport["repeat_head"],
          transport["seal_sha256"]) == (g, i["repeat_run"], i["repeat_attempt"], i["repeat_head"], sha(stage/"seal.json")),
         "Repeat transport not authenticated")
    marker = root / HISTORY / str(g) / "adoption.json"
    if marker.exists():
        old = campaign.read(marker)
        need(old["seal_sha256"] == sha(stage / "seal.json") and old["transport_sha256"] == sha(stage.parent / "transport.json"),
             "A different recovery run already supplied this group")
        return
    source, coverage = heldout_original.registry(root)
    frozen = ingest.Frozen(root)
    ids = frozen.groups[g]
    need(i["assigned"] == ids and i["assignment"] == frozen.sources["inventory/shards.json"] and
         i["science"] == frozen.science and set(seal["observations"]) == set(map(str, ids)), "Frozen repeat scope mismatch")
    original_heldout = adoption_state(root, g, job, ids, source, coverage)
    heldout_original.run(root, "verify")
    need(not (root/"history/ned-untyped-v1").exists() or not any((root/"history/ned-untyped-v1").rglob("result.json")),
         "Adopt the independent repeat before any post-primary correction")
    old_art = source["retained_event_artifacts"].get(f"chandra-event-inputs-{g}")
    used = campaign.read(stage / "original/retained.json")["provider"]
    keys = ("id", "name", "size_in_bytes", "digest")
    need((old_art is None and used is None) or (old_art and used and all(old_art[k] == used[k] for k in keys)),
         "Terminal index and retained original EVT2 artifact disagree")
    scopes = {}
    files = {}
    for path, file in sorted(ingest.logical_files(stage / "science").items()):
        p = Path(path)
        match = re.fullmatch(r"observations/(\d+)\.json\.gz", path)
        if match:
            oid = int(match[1])
            need(oid in ids and str(oid) == match[1], "Unexpected repeated ObsID")
            observed = ingest.ledger(file, oid)
            record = seal["observations"][str(oid)]
            need(record["ledger_sha256"] == sha(file) and record["status"] == observed["status"] and
                 (record["original_retained_products_accounted"] or observed["status"] not in {"complete","followup_incomplete"}),
                 "Repeated ObsID receipt mismatch")
            old = json.loads(gzip.decompress(frozen.baseline(path))) if path in frozen.blobs else None
            chosen = ingest.revisit(old, observed)
            scopes[str(oid)] = {"selection": chosen, "repeat_sha256": sha(file), "repeat_status": observed["status"],
                                "repeat_original_retained_matches": record["original_retained_products_accounted"]}
            if chosen != "worker":
                need((root/path).is_file() and sha(root/path) == hashlib.sha256(frozen.baseline(path)).hexdigest(),
                     "Preserved prior checkpoint differs")
                continue
        elif not (len(p.parts) >= 3 and
                  (p.parts[:2] == ("catalogs", f"upload-recovery-{i['repeat_run']}-{i['repeat_attempt']}-{g}") or
                   (p.parts[0] == "photometry" and p.parts[1] in {str(oid) for oid in ids}))):
            raise ValueError("Unowned recovery path: " + path)
        files[path] = file
    need(set(scopes) == set(map(str, ids)), "Missing repeated ObsID receipt")
    # Validate all collisions up front, then atomically install each file. A retry may see a partial install.
    for path, file in files.items():
        dest = root / path
        if dest.exists() and sha(dest) != sha(file):
            need(path in frozen.blobs and sha(dest) == hashlib.sha256(frozen.baseline(path)).hexdigest(),
                 "Refusing to overwrite unrelated evidence at " + path)
    for path, file in files.items():
        dest = root / path
        if dest.exists() and sha(dest) == sha(file):
            continue
        if path in frozen.blobs:
            save_blob(marker.parent / "before" / path, frozen.baseline(path))
        ingest.install(file, dest)
    save(marker, {"group": g, "assignment_commit": ingest.PLAN, "original_run": ingest.RUN, "original_job": job,
         "original_upload_status": "missing_artifact", "repeat_run": i["repeat_run"], "repeat_code": i["repeat_head"],
         "seal_sha256": sha(stage / "seal.json"), "transport_sha256": sha(stage.parent / "transport.json"),
         "terminal_coverage_sha256": sha(root/ingest.HISTORY/"coverage.json"),
         "original_heldout_sha256": sha(original_heldout), "observations": scopes,
         "scope": "Independent repeat supplied dated primary inputs to the corrective plan. Original worker transport and held-out outcomes remain registered separately."})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("recover", "import", "adopt"))
    p.add_argument("--evidence", type=Path, default=campaign.EVIDENCE)
    p.add_argument("--stage", "--output", dest="stage", type=Path)
    p.add_argument("--group", type=int)
    p.add_argument("--job", type=int)
    p.add_argument("--run", type=int)
    p.add_argument("--attempt", type=int, default=1)
    p.add_argument("--sha")
    a = p.parse_args()
    need(os.getenv("GITHUB_ACTIONS") == "true" and os.getenv("GITHUB_REPOSITORY") == ingest.REPO and
         os.getenv("GITHUB_REF") == "refs/heads/" + ingest.BRANCH and
         os.environ.get("GITHUB_WORKFLOW_REF", "").split("@", 1)[0] == ingest.REPO + "/" + ingest.FLOW,
         "Use the Chandra research workflow")
    root = a.evidence.resolve()
    if a.mode == "recover":
        need(a.stage and a.group is not None and a.job and os.getenv("RETAIN_EVENTS") and os.getenv("CIAO_PREFIX"),
             "Supply group, original job ID, output, RETAIN_EVENTS and CIAO_PREFIX")
        recover(root, a.stage.resolve(), a.group, a.job)
    elif a.mode == "import":
        need(a.group is not None and a.run and a.sha, "Supply group and repeat producer run/attempt/SHA")
        import_repeat(root, a.group, a.run, a.attempt, a.sha)
    else:
        need(a.stage is not None, "Supply a previously authenticated stage")
        adopt(root, a.stage.resolve())


if __name__ == "__main__":
    main()
