"""Checkpoint individual frozen _scan calls in large, unfinished ACIS products."""
import argparse
import gzip
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import tempfile
from unittest.mock import patch

import numpy as np

import campaign
import search

FROZEN_SEARCH = "17bdbb0ec1fa5c86bf1a3394b8df29d3b245e436dc0bb8ce2376a5c9e33a4f99"
FROZEN_REQUIREMENTS = "f3ee4d1c5cbf58c313adf4c3926484e4e3f5736834225d3dbf5976f518660790"
FROZEN_CAMPAIGN = "654b5d6aba7592f4c4b385b75a275b34d48c4e1d428ef555b55939ea818f94b5"
SUFFIX = ".scan-unit.json.gz"
ORIGINAL_SCAN = search._scan


def encode(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def atomic(path, data):
    """Never replace an existing committed unit, even with concurrent writers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() == data:
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        temp = Path(f.name)
        try:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            temp.unlink(missing_ok=True)
            raise
    try:
        try:
            os.link(temp, path)
        except FileExistsError:
            if path.read_bytes() != data:
                raise ValueError(f"Conflicting completed result at {path.name}")
    finally:
        temp.unlink(missing_ok=True)


def environment():
    identity = campaign.identity()
    if (identity["search.py"] != FROZEN_SEARCH or identity["campaign.py"] != FROZEN_CAMPAIGN
            or identity["requirements.txt"] != FROZEN_REQUIREMENTS):
        raise ValueError("Unknown original search implementation or environment prescription")
    required = dict(line.split("==") for line in (campaign.HERE / "requirements.txt").read_text().splitlines() if line)
    actual = {name:version(name) for name in required}
    if actual != required:
        raise ValueError("Installed package versions differ from frozen requirements")
    return {"search_sha256": FROZEN_SEARCH, "requirements_sha256": FROZEN_REQUIREMENTS,
            "campaign_sha256": FROZEN_CAMPAIGN, "wrapper_sha256": file_sha(Path(__file__)),
            "packages": actual, "python": platform.python_version(), "architecture": platform.machine(),
            "blas_threads": os.getenv("OPENBLAS_NUM_THREADS"), "omp_threads": os.getenv("OMP_NUM_THREADS")}


def repository_guard():
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REPOSITORY") != "dovvnloading/Open-Space"
        or os.getenv("GITHUB_REF") != "refs/heads/research/chandra-acis-fast-transients"
        or os.getenv("GITHUB_EVENT_NAME") not in {"workflow_dispatch", "push"}):
        raise RuntimeError("Use the Chandra research branch workflow")


def unpack(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def verify_evt(path, receipt):
    if (path.name != receipt["filename"] or path.stat().st_size != receipt["bytes"]
            or file_sha(path) != receipt["sha256"]):
        raise ValueError("EVT2 filename, byte count or hash differs from the source receipt")


def product(ledger, filename):
    acquired = [r for r in ledger.get("acquired", []) if r["filename"] == filename]
    products = [r for r in ledger.get("products", []) if r["input"]["filename"] == filename]
    if len(acquired) != 1 or len(products) != 1 or products[0]["input"] != acquired[0]:
        raise ValueError("A unique, fully acquired product with a matching ledger receipt is required")
    if products[0].get("search"):
        raise ValueError("This product already has a search checkpoint")
    return acquired[0], products[0]


class Calls:
    """Call identity comes from the unchanged stock loop order and generated arrays."""
    def __init__(self, handler):
        self.handler = handler
        self.descriptors = []
        self.segments = []
        self.stream_count = 0
        self.injection_count = 0

    def __call__(self, data, seg, meta, wcs, grids=None):
        geom = [seg.ccd, seg.start, seg.stop, seg.duty, seg.frame,
                None if seg.orientation is None else seg.orientation.tolist()]
        if grids is None:
            stream = f"injection-{self.injection_count:03d}"
            self.injection_count += 1
        else:
            if not self.segments or self.segments[-1] != geom:
                if self.segments and self.stream_count != 1 + campaign.N_CONTROL:
                    raise ValueError("Segment ended before all control calls")
                if geom in self.segments or self.injection_count:
                    raise ValueError("Segment ordering changed")
                self.segments.append(geom)
                self.stream_count = 0
            trial = self.stream_count
            self.stream_count += 1
            stream = f"segment-{len(self.segments)-1:04d}-" + (f"control-{trial-1:02d}" if trial else "real")
        if not data.flags.c_contiguous or data.dtype != search.EVENT_DTYPE:
            raise ValueError("Unexpected event representation")
        desc = {"id":stream, "ordinal":len(self.descriptors), "segment":geom, "rows":len(data),
                "event_rows_sha256":digest(memoryview(data).cast("B") if len(data) else b"")}
        self.descriptors.append(desc)
        return self.handler(desc, data, seg, meta, wcs, grids)


def walk(evt, seed, handler):
    calls = Calls(handler)
    with patch.object(search, "_scan", calls):
        result = search.search_file(evt, seed=seed, n_controls=campaign.N_CONTROL, n_injections=campaign.N_INJECT)
    if calls.segments and calls.stream_count != 1 + campaign.N_CONTROL:
        raise ValueError("Final real/control segment was incomplete")
    actual = [j["trial"] for j in result["injections"] if j["attempted"]]
    if len(actual) != calls.injection_count:
        raise ValueError("Injection call/attempt mismatch")
    return result, calls.descriptors, actual


def skeleton(result):
    # A substituted empty scan result influences only these output fields in
    # the frozen search.py; it never controls a future _scan call or RNG draw.
    return {**result, "triggers":[], "controls":[{**c, "triggers":[]} for c in result["controls"]],
            "injections":[{**j, "recovered":None} for j in result["injections"]]}


def plan(ledger_path, evt, destination, n_groups):
    ledger = campaign.read_obs(ledger_path)
    receipt, _ = product(ledger, evt.name)
    if ledger["obsid"] != int(receipt["obsid"]):
        raise ValueError("Wrong ObsID")
    verify_evt(evt, receipt)
    env = environment()
    seed = int(hashlib.sha256(("chandra-acis-fxt-v1:" + receipt["sha256"]).encode()).hexdigest()[:16], 16) >> 1
    empty, units, injections = walk(evt, seed, lambda *args: [])
    if int(empty["obsid"]) != ledger["obsid"]:
        raise ValueError("EVT2 header and observation ledger disagree")
    if not 1 <= n_groups <= 192:
        raise ValueError("Use one through 192 work groups")
    groups = [[] for _ in range(min(len(units), n_groups))]
    cost = [0] * len(groups)
    for u in sorted(units, key=lambda r:(-r["rows"],r["ordinal"])):
        k = min(range(len(groups)), key=cost.__getitem__)
        groups[k].append(u["id"])
        cost[k] += max(1,u["rows"])
    data = {"environment":env, "obsid":ledger["obsid"], "evt2":receipt, "seed":seed,
            "original_ledger_sha256":file_sha(ledger_path), "units":units, "groups":groups,
            "attempted_injection_trials_in_call_order":injections, "empty_result_sha256":digest(encode(skeleton(empty))),
            "controls":campaign.N_CONTROL, "injections":campaign.N_INJECT}
    atomic(destination, encode(data))
    return data


def validate_plan(manifest, evt):
    data = json.loads(manifest.read_bytes())
    if (data["environment"] != environment() or data["controls"] != campaign.N_CONTROL
            or data["injections"] != campaign.N_INJECT):
        raise ValueError("Plan, numerical environment or source code changed")
    verify_evt(evt, data["evt2"])
    ids = [d["id"] for d in data["units"]]
    grouped = [name for group in data["groups"] for name in group]
    if len(set(ids)) != len(ids) or sorted(grouped) != sorted(ids):
        raise ValueError("Repeated or missing work assignment")
    return data, file_sha(manifest)


def payload(path, plan_id, desc):
    capsule = unpack(path)
    raw = capsule["result"]
    if (capsule["plan_sha256"] != plan_id or capsule["unit"] != desc
            or digest(encode(raw)) != capsule["result_sha256"]):
        raise ValueError(f"Inconsistent or corrupt scan unit {desc['id']}")
    restored = []
    for t, indices in raw:
        if (not isinstance(indices, list) or len(indices) != t["on_counts"]
                or any(type(i) is not int or not 0 <= i < desc["rows"] for i in indices)):
            raise ValueError("Malformed photon indices in unit")
        restored.append((t, np.asarray(indices, dtype=int)))
    return restored


def existing_files(staged, plan_id, expected):
    files = {}
    for path in staged.rglob("*" + SUFFIX):
        if path.parent.name != plan_id:
            continue
        ident = path.name.removesuffix(SUFFIX)
        if ident not in expected:
            raise ValueError("Unexpected unit file in this plan")
        if ident in files and file_sha(files[ident]) != file_sha(path):
            raise ValueError(f"Workers disagree on duplicate unit {ident}")
        files[ident] = path
    return files


def attest(output, plan_id, task, checksum):
    run = os.getenv("GITHUB_RUN_ID", "fixture")
    attempt = os.getenv("GITHUB_RUN_ATTEMPT", "1")
    receipt = {"plan_sha256":plan_id, "task":task, "output_sha256":checksum, "run_id":run,
               "run_attempt":attempt, "commit":os.getenv("GITHUB_SHA"), "job":os.getenv("GITHUB_JOB")}
    atomic(output / "executions" / plan_id / f"{run}-{attempt}" / (task + ".json"), encode(receipt))


def execute(manifest, evt, output, *, group=None, unit=None, staged=None, ledger=None):
    data, plan_id = validate_plan(manifest, evt)
    expected = {u["id"]:u for u in data["units"]}
    assembling = staged is not None
    if assembling:
        files = existing_files(staged, plan_id, expected)
        missing = sorted(set(expected) - files.keys())
        if missing:
            raise ValueError(f"Missing {len(missing)} units. First missing: {missing[:3]}")
        selected = set()
    else:
        if (group is None) == (unit is None):
            raise ValueError("Select one work group or one unit")
        if group is not None:
            if group not in range(len(data["groups"])):
                raise ValueError("Unknown work group")
            selected = set(data["groups"][group])
        else:
            if unit not in expected:
                raise ValueError("Unknown work unit")
            selected = {unit}
    def dispatch(desc, events, seg, meta, wcs, grids):
        ident = desc["id"]
        if desc != expected.get(ident):
            raise ValueError(f"Original scan input or segment changed for {ident}")
        if assembling:
            return payload(files[ident], plan_id, desc)
        if ident in selected:
            path = output / plan_id / (ident + SUFFIX)
            if path.exists():
                payload(path, plan_id, desc)
            else:
                result = ORIGINAL_SCAN(events, seg, meta, wcs, grids)
                rows = [[event, ids.tolist()] for event, ids in result]
                capsule = {"plan_sha256":plan_id, "unit":desc, "result_sha256":digest(encode(rows)), "result":rows}
                atomic(path, gzip.compress(encode(capsule), compresslevel=6, mtime=0))
                print(f"Saved {ident}", flush=True)
            attest(output, plan_id, ident, file_sha(path))
        return []
    result, seen, attempted = walk(evt, data["seed"], dispatch)
    if (seen != data["units"] or attempted != data["attempted_injection_trials_in_call_order"]
            or digest(encode(skeleton(result))) != data["empty_result_sha256"]):
        raise ValueError("Original scan call sequence, support or injection plan changed")
    if not assembling:
        return
    raw = encode(result)
    packed = gzip.compress(raw, compresslevel=6, mtime=0)
    atomic(output / "search.json.gz", packed)
    receipt = {"plan_sha256":plan_id, "event_sha256":data["evt2"]["sha256"], "search_json_sha256":digest(raw),
               "compressed_sha256":digest(packed), "unit_files_sha256":{i:file_sha(p) for i,p in files.items()},
               "environment":data["environment"]}
    if ledger is not None:
        if file_sha(ledger) != data["original_ledger_sha256"]:
            raise ValueError("Original ledger changed; cannot construct a seeded worker ledger")
        obs = campaign.read_obs(ledger)
        p, record = product(obs, evt.name)
        if p != data["evt2"]:
            raise ValueError("Ledger product differs from the plan")
        record.update(search=result, complete=False)
        obs["status"] = "error"
        obs["errors"].append("Exact raw scan assembled from completed units. Catalog and photometry follow-up pending.")
        name = Path("observations") / f"{obs['obsid']}.json.gz"
        raw_ledger = gzip.compress(encode(obs), compresslevel=6, mtime=0)
        atomic(output / name, raw_ledger)
        receipt.update(seeded_ledger=str(name), seeded_ledger_sha256=digest(raw_ledger),
                       original_ledger_sha256=data["original_ledger_sha256"])
    atomic(output / "assembly.json", encode(receipt))
    attest(output, plan_id, "assembly", digest(encode(receipt)))
    return result


def main():
    repository_guard()
    args = argparse.ArgumentParser(description=__doc__)
    sub = args.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("plan")
    p.add_argument("ledger", type=Path)
    p.add_argument("evt2", type=Path)
    p.add_argument("manifest", type=Path)
    p.add_argument("--groups", type=int, default=8)
    w = sub.add_parser("work")
    a = sub.add_parser("assemble")
    for cmd in (w,a):
        cmd.add_argument("manifest", type=Path)
        cmd.add_argument("evt2", type=Path)
    w.add_argument("output", type=Path)
    choice = w.add_mutually_exclusive_group(required=True)
    choice.add_argument("--group", type=int)
    choice.add_argument("--unit")
    a.add_argument("staged", type=Path)
    a.add_argument("output", type=Path)
    a.add_argument("--ledger", type=Path)
    opt = args.parse_args()
    if opt.mode == "plan":
        p = plan(opt.ledger, opt.evt2, opt.manifest, opt.groups)
        if os.getenv("GITHUB_OUTPUT"):
            with open(os.environ["GITHUB_OUTPUT"], "a") as target:
                print("matrix=" + json.dumps(list(range(len(p["groups"])))), file=target)
        print(f"{len(p['units'])} calls in {len(p['groups'])} worker groups.")
    elif opt.mode == "work":
        execute(opt.manifest, opt.evt2, opt.output, group=opt.group, unit=opt.unit)
    else:
        execute(opt.manifest, opt.evt2, opt.output, staged=opt.staged, ledger=opt.ledger)


if __name__ == "__main__":
    main()
