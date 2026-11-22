"""Recover and authenticate frozen primary Actions artifacts, one producer at a time.

In the Chandra workflow, run index, then batch --groups 0-15 (or supplement),
then audit. Repeat batches as required. --zip-dir accepts original, unopened
<artifact-name>.zip or <artifact-id>.zip files. The GitHub SHA-256 applies to
those exact ZIP bytes; pre-extracted directories cannot establish that digest.
"""
import argparse
from contextlib import contextmanager
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tempfile
import time
from zipfile import ZipFile

import campaign
import catalog_pack
from followup import raw_scan
from migrate import INVENTORY, verify_scan

REPO = "dovvnloading/Open-Space"
BRANCH = "research/chandra-acis-fast-transients"
RUN = 44985879443
ATTEMPT = 1
CODE = "4663ab98666b5d9b5b116478f61d6cae5d49e104"
PLAN = "0655e5b54215a3d7bbac1f6b2561379bb50c81cf"
FLOW = ".github/workflows/chandra-acis-fxt.yml"
PREFIX = "results/chandra_acis_fxt/evidence"
HISTORY = Path("primary_ingest") / f"{RUN}-{ATTEMPT}"
SUP = "chandra-supplemental-search"
HEX = re.compile("[0-9a-f]{64}")
EXPECTED = {f"chandra-evidence-{g}":str(g) for g in range(192)} | {SUP:"supplement"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return catalog_pack.sha(path)


def save(path, obj, *, zipped=False):
    data = (json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    if zipped:
        data = gzip.compress(data, compresslevel=6, mtime=0)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".new")
    tmp.write_bytes(data)
    tmp.replace(path)


def git(*args, **kwargs):
    return subprocess.check_output(["git", "-C", str(campaign.HERE.parents[1]), *args], **kwargs)


def immutable(rev):
    fast = subprocess.run(["git","-C",str(campaign.HERE.parents[1]),"cat-file","-e",rev+"^{commit}"],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if fast.returncode:
        git("fetch", "--no-tags", "--depth=1", "--filter=blob:none", "origin", rev)
    require(git("rev-parse", rev+"^{commit}").decode().strip() == rev, "Incorrect immutable commit")


def blob(commit, path):
    return git("show", commit + ":" + path)


class Frozen:
    def __init__(self, root):
        immutable(PLAN)
        immutable(CODE)
        tree = git("ls-tree", "-r", "-z", PLAN + ":" + PREFIX)
        self.blobs = {}
        for row in tree.split(b"\0"):
            if row:
                meta, name = row.split(b"\t", 1)
                mode, kind, oid = meta.split()
                if kind == b"blob":
                    self.blobs[name.decode()] = oid.decode()
        sources = {}
        for name in ("inventory/shards.json", "inventory/inventory.json", "inventory/campaign.json",
                     "release_supplement/supplement.json"):
            raw = git("cat-file", "blob", self.blobs[name])
            sources[name] = {"git_blob":self.blobs[name], "sha256":hashlib.sha256(raw).hexdigest()}
            if name == "inventory/shards.json":
                self.groups = json.loads(raw)
                # A later campaign may have replanned the working copy.
            else:
                require(sha(root/name) == sources[name]["sha256"], f"Conflicting frozen input {name}")
        require(sources["inventory/inventory.json"]["sha256"] == INVENTORY, "Unexpected parent census")
        self.population = {r["obsid"] for r in campaign.read(root/"inventory/inventory.json")}
        self.supplement = {int(r["obsid"]) for r in campaign.read(root/"release_supplement/supplement.json")["unknown"]}
        assigned = [n for group in self.groups for n in group]
        require(len(self.groups) == 192 and len(set(assigned)) == len(assigned) and set(assigned) <= self.population
                and self.population.isdisjoint(self.supplement) and len(self.supplement) == 14,
                "Invalid dated or supplemental assignments")
        self.assigned = set(assigned)
        self.science = campaign.read(root/"inventory/campaign.json")["implementation"]
        base = PREFIX.rsplit("/", 1)[0]
        for name, expected in self.science.items():
            require(hashlib.sha256(blob(PLAN, base+"/"+name)).hexdigest() == expected
                    == hashlib.sha256(blob(CODE, base+"/"+name)).hexdigest(), "Scientific code differs: " + name)
        self.ambiguous_sha = hashlib.sha256(blob(CODE, base+"/ambiguous.py")).hexdigest()
        self.sources = sources

    def baseline(self, path):
        return git("cat-file", "blob", self.blobs[path])


def api(endpoint, destination=None):
    route = f"repos/{REPO}/actions/{endpoint}"
    for i in range(3):
        args = ["gh", "api", "--method", "GET", "-H", "X-GitHub-Api-Version: 2022-11-28", route]
        try:
            if destination is None:
                return json.loads(subprocess.check_output(args))
            with destination.open("wb") as out:
                subprocess.run(args, stdout=out, check=True)
            return
        except (subprocess.CalledProcessError, json.JSONDecodeError):
            if i == 2:
                raise
            time.sleep(2 ** i)


def listed_page_set(endpoint, field):
    output = []
    page = 1
    # An active run can add artifacts during pagination. Retry any incoherent snapshot.
    while True:
        data = api(f"{endpoint}?per_page=100&page={page}")
        output.extend(data[field])
        if len(data[field]) < 100:
            break
        page += 1
    require(len(output) == data["total_count"] and len({r["id"] for r in output}) == len(output),
            "GitHub pagination changed while reading; repeat index")
    return output


def pages(endpoint, field):
    for attempt in range(3):
        try:
            return listed_page_set(endpoint,field)
        except ValueError:
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def index(root, frozen):
    run = api(f"runs/{RUN}")
    require(run["id"] == RUN and run["run_attempt"] == ATTEMPT and run["head_sha"] == CODE
            and run["head_branch"] == BRANCH and run["path"].split("@")[0] == FLOW
            and run["repository"]["full_name"] == REPO, "Unexpected source Actions run")
    listings = pages(f"runs/{RUN}/artifacts", "artifacts")
    jobs = pages(f"runs/{RUN}/attempts/{ATTEMPT}/jobs", "jobs")
    selected, inputs, workers = {}, {}, {}
    for entry in listings:
        name = entry["name"]
        if name not in EXPECTED and not re.fullmatch(r"chandra-(?:event-inputs-\d+|supplemental-event-inputs)", name):
            continue
        parent = entry["workflow_run"]
        require(parent["id"] == RUN and parent["head_sha"] == CODE and parent["head_branch"] == BRANCH
                and parent["repository_id"] == run["repository"]["id"] == parent["head_repository_id"],
                "Cross-run artifact "+name)
        dest = selected if name in EXPECTED else inputs
        require(name not in dest, "Duplicate source artifact: "+name)
        require(re.fullmatch(r"sha256:[0-9a-f]{64}", entry.get("digest") or ""), "Absent GitHub artifact digest: "+name)
        dest[name] = {k:entry[k] for k in ("id","name","size_in_bytes","digest","expired","created_at","expires_at")}
    for job in jobs:
        match = re.fullmatch(r"observation \((?:group: )?(\d+)\)", job["name"])
        group = match[1] if match else "supplement" if job["name"] == "supplement" else None
        if group is not None:
            require(group not in workers and (group == "supplement" or int(group) in range(192)), "Duplicate source job")
            workers[group] = {k:job.get(k) for k in ("id","name","status","conclusion","started_at","completed_at","html_url")}
    data = {"run":{k:run[k] for k in ("id","run_attempt","head_sha","head_branch","path","status","conclusion","html_url")},
            "plan_commit":PLAN, "assignment_sources":frozen.sources, "artifacts":selected,
            "retained_event_artifacts":inputs, "jobs":workers, "missing_artifacts":sorted(set(EXPECTED)-set(selected))}
    path = root/HISTORY/"source.json"
    if path.exists():
        previous = campaign.read(path)
        for name, prior in previous["artifacts"].items():
            require(name in selected and all(prior[k] == selected[name][k] for k in ("id","size_in_bytes","digest")),
                    "A previously indexed artifact disappeared or changed: "+name)
    save(path, data)
    return data


def selection(value):
    out = []
    for term in value.split(","):
        if term == "supplement":
            out.append(term)
            continue
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", term)
        require(match is not None, "Expected group numbers, ranges, or supplement")
        low, high = int(match[1]), int(match[2] or match[1])
        require(0 <= low <= high < 192, "Group outside frozen assignment")
        out.extend(str(i) for i in range(low, high+1))
    require(len(set(out)) == len(out), "Duplicate requested group")
    return out


def unzip(archive, destination):
    with ZipFile(archive) as bundle:
        files = [p for p in bundle.infolist() if not p.is_dir()]
        require(len({p.filename for p in files}) == len(files), "Duplicate ZIP member")
        require(sum(p.file_size for p in files) < shutil.disk_usage(destination).free, "Insufficient scratch for this artifact")
        for member in bundle.infolist():
            name = member.filename.removesuffix("/")
            p = PurePosixPath(name)
            require(p.parts and not p.is_absolute() and str(p) == name and ".." not in p.parts and "\\" not in name
                    and stat.S_IFMT(member.external_attr >> 16) in (0, stat.S_IFREG, stat.S_IFDIR),
                    "Unsafe original ZIP member: "+name)
            if member.is_dir():
                continue
            out = destination / str(p)
            out.parent.mkdir(parents=True, exist_ok=True)
            with bundle.open(member) as src, out.open("xb") as dst:
                shutil.copyfileobj(src, dst, length=1024*1024)
            require(out.stat().st_size == member.file_size, "Incomplete ZIP extraction")


def owned(path, group, allowed):
    p = PurePosixPath(path)
    if path in {f"environment-{group}.json", f"shard-{group}.json"}:
        return
    if len(p.parts) >= 3 and p.parts[:2] == ("catalogs", f"{RUN}-{ATTEMPT}-{group}"):
        return
    if len(p.parts) >= 3 and p.parts[0] == "photometry" and p.parts[1] in {str(n) for n in allowed}:
        return
    obs = re.fullmatch(r"observations/(\d+)\.json\.gz", path)
    require(obs and int(obs[1]) in allowed and str(int(obs[1])) == obs[1], f"Foreign artifact path in group {group}: {path}")


def ledger(path, oid):
    row = campaign.read_obs(path)
    require(row.get("obsid") == oid, f"ObsID mismatch in {path.name}")
    acquired = {p["sha256"]:p for p in row.get("acquired", [])}
    require(len(acquired) == len(row.get("acquired", [])), "Duplicated acquired product")
    scans = {}
    for product in row.get("products", []):
        rec = product["input"]
        h = rec["sha256"]
        require(h in acquired and h not in scans and rec["filename"] == acquired[h]["filename"]
                and int(rec["obsid"]) == oid and HEX.fullmatch(h), "Mismatch between acquired and processed products")
        scan = product.get("search")
        if scan:
            verify_scan(scan, rec, oid)
            durations = {"60","180","600","1800","3600"}
            support = scan.get("start_support_by_duration",{})
            exposure = scan.get("exposure_by_duration_s",{})
            gti = scan.get("gti_livetime_s")
            require(set(support) == set(exposure) == durations and isinstance(gti,(int,float)) and math.isfinite(gti)
                    and gti>=0 and all(isinstance(support[k],list) and math.isfinite(exposure[k]) and exposure[k]>=0 for k in durations)
                    and all(isinstance(s.get("triggers"),list) for s in [scan,*scan["controls"]])
                    and all(type(j.get("attempted")) is bool and (not j["attempted"] or type(j.get("recovered")) is bool)
                            for j in scan["injections"]), "Incomplete GTI, duration, trigger or injection receipt")
        scans[h] = product
    if row["status"] in {"complete","followup_incomplete"}:
        require(acquired and scans.keys() == acquired.keys() and all(p.get("complete") and p.get("search") for p in scans.values()),
                "Source claimed a complete scan without every acquired product")
    return row


def revisit(old, updated):
    """Prefer a nonregressing worker ledger and preserve both original byte streams."""
    if not old:
        return "worker"
    if old.get("acquired") and not updated.get("acquired"):
        return "checkpoint_after_interrupted_retry"
    if old.get("acquired"):
        require(old["acquired"] == updated["acquired"], "Retry changed the archived input bytes")
    after = {p["input"]["sha256"]:p for p in updated.get("products",[])}
    for p in old.get("products", []):
        fresh = after.get(p["input"]["sha256"])
        if p.get("search") and not (fresh and fresh.get("search")):
            return "checkpoint_after_interrupted_retry"
        if p.get("search"):
            require(raw_scan(p) == raw_scan(fresh), "Retry changed a frozen real scan, control, injection or exposure")
            for a,b in zip([p["search"],*p["search"]["controls"]],[fresh["search"],*fresh["search"]["controls"]],strict=True):
                for before,now in zip(a["triggers"],b["triggers"],strict=True):
                    known = before.get("vet",{})
                    if known and known.get("status") != "catalog_incomplete":
                        require(known == now.get("vet"), "Retry changed a completed original catalog classification")
                    if before.get("flux",{}).get("peak_flux_erg_cm2_s") is not None:
                        require(before["flux"] == now.get("flux"), "Retry changed an already calibrated peak")
                    if before.get("diagnostics",{}).get("files"):
                        require(before["diagnostics"] == now.get("diagnostics"), "Retry changed original detector evidence")
    if old.get("status") == "complete" and updated.get("status") != "complete":
        return "checkpoint_after_interrupted_retry"
    return "worker"


def logical_files(root):
    return {p.relative_to(root).as_posix():p for p in root.rglob("*") if p.is_file()}


def install(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as stream:
        temp = Path(stream.name)
    try:
        shutil.copyfile(source, temp)
        require(sha(source) == sha(temp), "Transit error installing "+str(target))
        temp.replace(target)
    finally:
        temp.unlink(missing_ok=True)


def publish(root, frozen, source, logical, group):
    target = catalog_pack.safe(root/logical, root)
    incoming = sha(source)
    existing = sha(target) if target.is_file() else None
    if logical in frozen.blobs:
        previous = frozen.baseline(logical)
        old_hash = hashlib.sha256(previous).hexdigest()
        require(existing in {None,incoming,old_hash}, "Destination differs from both the checkpoint and the artifact: "+logical)
        if old_hash != incoming:
            prior = root/HISTORY/"checkpoint_before"/group/logical
            if not prior.exists():
                prior.parent.mkdir(parents=True,exist_ok=True)
                temp=prior.with_name(prior.name+".new")
                temp.write_bytes(previous)
                temp.replace(prior)
            require(sha(prior) == old_hash, "Conflict in preserved checkpoint file")
    else:
        require(existing in {None,incoming}, "Refusing to overwrite different, non-checkpoint evidence: "+logical)
    if existing == incoming:
        return
    install(source, target)


def receipt_file(root, group):
    return root/HISTORY/"groups"/(group+".json")


def verify_import(root, frozen, group, entry):
    rec = campaign.read(receipt_file(root,group))
    require(rec["group"] == group and rec["artifact_id"] == entry["id"] and rec["zip_digest"] == entry["digest"]
            and rec["plan_commit"] == PLAN and rec["code_commit"] == CODE, "Conflicting stored import receipt")
    manifest = root/HISTORY/"members"/(group+".json.gz")
    require(sha(manifest) == rec["member_manifest_sha256"], "Changed import member manifest")
    members = campaign.read_obs(manifest)
    for line in members:
        target = root/line["destination"]
        if not target.exists() and line["destination"].startswith("catalogs/") and target.suffix == ".bin":
            catalog_pack.recorded_reply(root,target,{"size_bytes":line["bytes"],"sha256":line["sha256"]})
        else:
            require(target.is_file() and target.stat().st_size == line["bytes"] and sha(target) == line["sha256"],
                    "Original artifact member missing or changed: "+str(target))
    for name, meta in rec["checkpoint_before"].items():
        prior = root/HISTORY/"checkpoint_before"/group/name
        require(name in frozen.blobs and frozen.blobs[name] == meta["git_blob"] and sha(prior) == meta["sha256"],
                "Prior checkpoint snapshot changed: "+name)
    require(rec["source_member_count"] == len(members), "Original artifact member count changed")
    return rec


def validated_stage(frozen, stage, group, files):
    optional = []
    choices, statuses = {}, {}
    if group == "supplement":
        ids = frozen.supplement
        if "release_bounds.json" in files:
            header = campaign.read(files["release_bounds.json"])
            require(str(header["run_id"]) == str(RUN) and header["git_commit"] == CODE
                    and header["implementation"] == frozen.science and header["supplemental_program_sha256"] == frozen.ambiguous_sha
                    and header["supplement_input_sha256"] == frozen.sources["release_supplement/supplement.json"]["sha256"]
                    and sorted(r["obsid"] for r in header["obsids"]) == sorted(ids), "Supplement uses different scope or code")
        else:
            optional.append("missing_release_bounds")
        if "summary.json" in files:
            summary = campaign.read(files["summary.json"])
            require(sorted(r["obsid"] for r in summary["status"]) == sorted(ids), "Supplemental summary assignment changed")
        else:
            optional.append("missing_supplement_summary")
        prefix = Path("supplemental_search") / f"{RUN}-{ATTEMPT}"
        for name in files:
            p = PurePosixPath(name)
            if p.parts[0] == "observations":
                require(re.fullmatch(r"observations/\d+\.json\.gz",name) and int(p.name.split(".")[0]) in ids,
                        "Unknown supplemental ObsID")
            if p.parts[0] == "photometry":
                require(len(p.parts)>=3 and p.parts[1] in {str(n) for n in ids}, "Foreign supplemental photometry")
            if p.parts[0] == "catalogs":
                require(len(p.parts)>=3 and p.parts[1] == f"{RUN}-{ATTEMPT}", "Foreign supplemental cache")
    else:
        ids = set(frozen.groups[int(group)])
        for name in files:
            owned(name, group, ids)
        env = f"environment-{group}.json"
        marker = f"shard-{group}.json"
        if env in files:
            row = campaign.read(files[env])
            require(str(row["run_id"]) == str(RUN) and row["git_commit"] == CODE, "Wrong worker environment")
        else:
            optional.append("missing_environment")
        if marker in files:
            row = campaign.read(files[marker])
            require(row["group"] == int(group) and row["assigned"] == frozen.groups[int(group)]
                    and row["commit"] == CODE and row["completed_worker"] is True, "Changed worker assignments")
        else:
            optional.append("missing_worker_end_receipt")
        prefix = Path(".")
    for oid in sorted(ids):
        name = f"observations/{oid}.json.gz"
        prior = None
        if group != "supplement" and name in frozen.blobs:
            prior = json.loads(gzip.decompress(frozen.baseline(name)))
            require(prior["obsid"] == oid, "Incorrect checkpoint observation identity")
        if name not in files:
            statuses[str(oid)] = {"worker_status":"unattempted","checkpoint_status":prior["status"] if prior else None}
            continue
        incoming = ledger(files[name], oid)
        choice = revisit(prior, incoming)
        if choice != "worker":
            choices[str(oid)] = choice
        statuses[str(oid)] = {"worker_status":incoming["status"], "checkpoint_status":prior["status"] if prior else None,
                              "selected_source":choice, "selected_status":(prior if choice != "worker" else incoming)["status"]}
    return prefix, choices, statuses, optional


def import_one(root, frozen, source, group, packed, zip_dir):
    name = SUP if group == "supplement" else "chandra-evidence-"+group
    require(name in source["artifacts"], "No indexed artifact for "+name+"; refresh index or record a gap")
    entry = source["artifacts"][name]
    if receipt_file(root,group).exists():
        return verify_import(root,frozen,group,entry)
    with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as tmp:
        temp = Path(tmp)
        archived = None
        if zip_dir:
            for file in (zip_dir/(name+".zip"), zip_dir/(str(entry["id"])+".zip")):
                if file.is_file():
                    archived = file
                    break
        if archived is None:
            require(not entry["expired"], "Source artifact expired and no checked local ZIP was supplied")
            archived = temp/"download.zip"
            api(f"artifacts/{entry['id']}/zip",archived)
        require(archived.stat().st_size == entry["size_in_bytes"] and "sha256:"+sha(archived) == entry["digest"],
                "Artifact ZIP size or GitHub SHA-256 mismatch: "+name)
        stage = temp/"stage"
        stage.mkdir()
        unzip(archived,stage)
        files = logical_files(stage)
        prefix, choices, statuses, issues = validated_stage(frozen,stage,group,files)
        for oid in choices:
            name=f"observations/{oid}.json.gz"
            require((root/name).is_file() and sha(root/name)==hashlib.sha256(frozen.baseline(name)).hexdigest(),
                    "Cannot retain a missing or modified checkpoint ledger: "+name)
        def target_name(rel):
            logical = prefix/rel
            p = Path(rel)
            oid = p.parts[1].split(".")[0] if len(p.parts)>1 and p.parts[0] in {"observations","photometry"} else None
            if oid in choices:
                logical = HISTORY/"interrupted_worker"/group/rel
            return logical.as_posix()
        members = [{"source":n,"destination":target_name(n),"bytes":p.stat().st_size,"sha256":sha(p)}
                   for n,p in sorted(files.items())]
        # Keep the original source membership, even where the transferred representation is losslessly packed.
        if packed and group != "supplement":
            ns = f"{RUN}-{ATTEMPT}-{group}"
            if (stage/"catalogs"/ns).exists():
                catalog_pack.pack(stage,ns,catalog_pack.protected(stage),False,True)
        for rel,p in sorted(logical_files(stage).items()):
            publish(root,frozen,p,target_name(rel),group)
        if packed and group != "supplement" and (stage/"catalog_store").exists():
            readme = root/"catalog_store/README.md"
            if readme.exists():
                require(readme.read_text() == catalog_pack.README, "Conflicting catalog-store introduction")
            else:
                readme.write_text(catalog_pack.README)
        snapshots = {}
        base = root/HISTORY/"checkpoint_before"/group
        for rel,p in logical_files(base).items():
            require(rel in frozen.blobs and hashlib.sha256(frozen.baseline(rel)).hexdigest() == sha(p), "Invalid historical file")
            snapshots[rel] = {"git_blob":frozen.blobs[rel],"sha256":sha(p)}
        manifest = root/HISTORY/"members"/(group+".json.gz")
        save(manifest,members,zipped=True)
        receipt = {"group":group,"artifact_name":name,"artifact_id":entry["id"],"zip_bytes":entry["size_in_bytes"],
                   "zip_digest":entry["digest"],"code_commit":CODE,"plan_commit":PLAN,"assignment_sources":frozen.sources,
                   "source_member_count":len(members),"member_manifest_sha256":sha(manifest),"checkpoint_before":snapshots,
                   "worker_ledger_status":statuses,"missing_worker_receipts":issues,"packed_catalogs":bool(packed and group!="supplement")}
        # Verify the installed logical files before marking this producer durably imported.
        for line in members:
            target=root/line["destination"]
            if not target.exists() and line["destination"].startswith("catalogs/") and target.suffix==".bin":
                catalog_pack.recorded_reply(root,target,{"size_bytes":line["bytes"],"sha256":line["sha256"]})
            else:
                require(target.is_file() and target.stat().st_size==line["bytes"] and sha(target)==line["sha256"],
                        "Failed final member transfer: "+line["destination"])
        save(receipt_file(root,group),receipt)
        return receipt


def audit(root,frozen,source):
    rows, groups = {}, {}
    transport_gaps, worker_gaps, residual = [], [], []
    for name,g in EXPECTED.items():
        job = source["jobs"].get(g,{})
        entry = source["artifacts"].get(name)
        rec = verify_import(root,frozen,g,entry) if entry and receipt_file(root,g).is_file() else None
        state = "verified" if rec else "missing_artifact" if entry is None else "available_unimported"
        if state != "verified":
            transport_gaps.append(g)
        assigned = sorted(frozen.supplement) if g=="supplement" else frozen.groups[int(g)]
        markers = rec["missing_worker_receipts"] if rec else ["no_verified_artifact"]
        if not job or job.get("conclusion") != "success" or markers:
            worker_gaps.append(g)
        statuses = rec["worker_ledger_status"] if rec else {}
        missing = []
        for oid in assigned:
            s = statuses.get(str(oid),{"worker_status":"unverified"})
            if s["worker_status"] in {"unattempted","unverified"}:
                missing.append(oid)
            if g!="supplement":
                rows[str(oid)] = {"group":int(g),**s}
                if s.get("selected_status") != "complete":
                    residual.append(oid)
        groups[g] = {"artifact_status":state,"artifact_id":entry["id"] if entry else None,"job":job,
                     "missing_worker_receipts":markers,"assigned_without_verified_worker_ledger":missing}
        if missing and g not in worker_gaps:
            worker_gaps.append(g)
    for oid in sorted(frozen.population - frozen.assigned):
        rel=f"observations/{oid}.json.gz"
        require(rel in frozen.blobs, "Unassigned census ObsID lacks checkpoint ledger")
        old = frozen.baseline(rel)
        require((root/rel).is_file() and sha(root/rel) == hashlib.sha256(old).hexdigest(), "Unassigned checkpoint ledger changed")
        row = ledger(root/rel,oid)
        require(row["status"] == "complete", "Unassigned census ObsID lacked completed checkpoint status")
        rows[str(oid)] = {"group":None,"selected_source":"immutable_checkpoint","selected_status":row["status"]}
    require(len(rows)==len(frozen.population), "Incomplete census assignment accounting")
    final = {"source_run":RUN,"attempt":ATTEMPT,"source_code":CODE,"assignment_commit":PLAN,
             "source_status":source["run"],"transport_gaps":transport_gaps,"incomplete_worker_groups":worker_gaps,
             "dated_residual_obsids":sorted(residual),"dated_obsids":rows,"groups":groups,
             "complete_artifact_transport":not transport_gaps,
             "all_source_workers_finished_assignment":not worker_gaps and source["run"]["status"]=="completed",
             "scope":"Worker statuses describe original products. The scientific reducer and repair plan verify all exposure, controls, injections and follow-up separately."}
    save(root/HISTORY/"coverage.json",final)
    return final


@contextmanager
def lock(root):
    lockfile=Path(os.getenv("RUNNER_TEMP",tempfile.gettempdir()))/("chandra-import-"+hashlib.sha256(str(root).encode()).hexdigest()+".lock")
    with lockfile.open("a") as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        yield


def main():
    require(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==REPO
            and os.getenv("GITHUB_REF")=="refs/heads/"+BRANCH and os.getenv("GITHUB_EVENT_NAME") in {"push","workflow_dispatch"}
            and os.getenv("GITHUB_WORKFLOW_REF","").split("@")[0]==REPO+"/"+FLOW, "Use the single Chandra research workflow")
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--evidence",type=Path,default=campaign.EVIDENCE)
    sub=cli.add_subparsers(dest="mode",required=True)
    sub.add_parser("index")
    b=sub.add_parser("batch")
    b.add_argument("--groups",required=True)
    b.add_argument("--pack-catalogs",action="store_true")
    b.add_argument("--zip-dir",type=Path)
    b.add_argument("--available-only",action="store_true",help="skip source artifacts still missing; audit records every gap")
    a=sub.add_parser("audit")
    a.add_argument("--allow-partial",action="store_true",help="record explicit missing workers/artifacts in a bounded checkpoint")
    args=cli.parse_args()
    root=args.evidence.resolve()
    with lock(root):
        frozen=Frozen(root)
        if args.mode=="index":
            data=index(root,frozen)
            print(f"Indexed {len(data['artifacts'])}/{len(EXPECTED)} source evidence artifacts.")
            return
        data=campaign.read(root/HISTORY/"source.json")
        require(data["run"]["id"]==RUN and data["run"]["head_sha"]==CODE and data["plan_commit"]==PLAN
                and data["assignment_sources"]==frozen.sources,"Import plan or source metadata changed")
        if args.mode=="batch":
            for group in selection(args.groups):
                if args.available_only and (SUP if group=="supplement" else "chandra-evidence-"+group) not in data["artifacts"]:
                    print(f"Group {group} has no source artifact in the pinned run index.",flush=True)
                    continue
                rec=import_one(root,frozen,data,group,args.pack_catalogs,args.zip_dir)
                print(f"Verified group {group}, {rec['source_member_count']} original artifact members.",flush=True)
        else:
            result=audit(root,frozen,data)
            print(f"{len(result['transport_gaps'])} artifact gaps, {len(result['incomplete_worker_groups'])} unfinished source groups, "
                  f"{len(result['dated_residual_obsids'])} dated ObsIDs requiring follow-up.")
            if not args.allow_partial and not (result["complete_artifact_transport"] and result["all_source_workers_finished_assignment"]):
                raise SystemExit("Incomplete source campaign. Exact deficits are recorded in coverage.json.")


if __name__=="__main__":
    main()
