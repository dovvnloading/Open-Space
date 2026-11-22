"""Obtain versioned timing releases and reduce their published models with PINT.

fetch writes inputs.json under --out and extracted, relative paths under --work.
reduce starts from those inputs; no reduction is inferred from a missing file.
"""
import argparse
import contextlib
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile
import tempfile
import traceback
import zipfile
from collections import defaultdict

RELEASES = {
    "ng": ("16051178", "NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz",
           "557d42dd8486a5f8272d90dec9b228a8"),
    "epta": ("8300645", "EPTA-DR2.zip", "307b2b99cd352876b07409899448af8b"),
}
ALIASES = {"B1855+09": "J1857+0943", "B1937+21": "J1939+2134", "B1953+29": "J1955+2908"}
CLOCK_REV = "99a698887521a3000c1f37d027284eac15d92d9e"  # IPTA, NG15 v2.1.0 release day
PSR = re.compile(r"(?<![A-Za-z0-9])([JB]\d{4}[+-]\d{2,4})(?!\d)", re.I)


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n")


def digest(path, alg="sha256"):
    with path.open("rb") as stream:
        h = hashlib.new(alg)
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, path):
    subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location",
                    "--connect-timeout", "30", "--retry", "4", "--retry-all-errors",
                    "--speed-limit", "1024", "--speed-time", "120", "--max-time", "1800",
                    "--output", str(path), url], check=True)


def retained(path):
    return (path.suffix.lower() in {".tim", ".par", ".toa"}
            or bool({"clock", "tims", "tim", "toas"} & {p.lower() for p in path.parts[:-1]})
            or path.name.lower().startswith("readme"))


def emit_member(member, stream, root):
    path = PurePosixPath(member)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe archive path {member}")
    if not retained(path) or any(p.startswith(".") or p == "__MACOSX" for p in path.parts):
        return
    dst = root.joinpath(*path.parts)
    dst.parent.mkdir(parents=True, exist_ok=True)
    with stream() as source, dst.open("wb") as target:
        shutil.copyfileobj(source, target)


def unpack(archive, root):
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as z:
            for item in z.infolist():
                if not item.is_dir():
                    emit_member(item.filename, lambda: z.open(item), root)
    else:
        with tarfile.open(archive, mode="r|gz") as z:
            for item in z:
                if item.isfile():
                    emit_member(item.name, lambda: z.extractfile(item), root)


def dataset_for(path, source):
    parts = [p.lower() for p in path.parts]
    if source == "epta":
        # DR2new is nested data selection. DR2full+ adds InPTA measurements.
        if any("new" in p or "inpta" in p or "+" == p[-1:] or "fullplus" in p
               for p in parts if not PSR.search(p)):
            return None
        return "epta"
    for tag, dataset in (("narrowband", "nb"), ("wideband", "wb")):
        if tag in parts:
            return dataset
    return None


def rank(path, dataset, suffix):
    full = str(path).lower()
    parts = {p.lower() for p in path.parts}
    name = path.name.lower()
    if suffix == ".par":
        if any(tag in full for tag in ("predictive", "norednoise", "no_red_noise", "no-red-noise", "template")):
            return None
        return ("full" in full if dataset == "epta" else True,
                bool(PSR.fullmatch(path.stem)) if dataset == "epta" else True,
                "pint" in full, "tempo2" not in full and ".t2." not in name,
                "par" in parts, -len(path.parts))
    if suffix == ".tim":
        return ("full" in full if dataset == "epta" else True,
                "excise" in full, "all" in name or "combined" in name,
                bool(PSR.fullmatch(path.stem)) if dataset == "epta" else True,
                "tim" in parts, -len(path.parts))
    return None


def unique_best(files, dataset, suffix):
    candidates = [(rank(p, dataset, suffix), p) for p in files]
    candidates = [(s, p) for s, p in candidates if s is not None]
    if not candidates:
        raise ValueError(f"missing {suffix}")
    highest = max(s for s, _ in candidates)
    best = sorted(p for s, p in candidates if s == highest)
    if len(best) > 1 and len({digest(p) for p in best}) != 1:
        raise ValueError("ambiguous published top-level input: " + ", ".join(map(str, best)))
    return best[0]


def timing_includes(path, root, active=()):
    """Hash the exact PINT INCLUDE closure before handing any file to PINT."""
    if path in active:
        raise ValueError("cyclic INCLUDE in " + str(path))
    files = {}
    with path.open("rb") as stream:
        for raw in stream:
            line = raw.lstrip()
            if not line[:7].upper().startswith(b"INCLUDE"):
                continue
            fields = line.split()
            if len(fields) < 2 or fields[0].upper() != b"INCLUDE":
                continue
            child = (path.parent / os.fsdecode(fields[1])).resolve()
            if not child.is_relative_to(root) or not child.is_file():
                raise ValueError(f"missing or out-of-release INCLUDE {child} in {path}")
            key = str(child.relative_to(root))
            files[key] = digest(child)
            files.update(timing_includes(child, root, (*active, path)))
    return files


def inventory(work):
    by_name = defaultdict(lambda: defaultdict(list))
    for source in RELEASES:
        for path in sorted((work / source).rglob("*")):
            if path.suffix.lower() not in {".par", ".tim"}:
                continue
            dataset = dataset_for(path.relative_to(work), source)
            matches = PSR.findall(path.name) or PSR.findall(str(path.parent.relative_to(work)))
            if not dataset or not matches:
                continue
            name = matches[-1].upper()
            by_name[dataset, ALIASES.get(name, name)][path.suffix.lower()].append(path)
    result = []
    for (dataset, name), candidates in sorted(by_name.items()):
        # Ignore backend-only directories if there is no released model.
        if ".par" not in candidates:
            continue
        row = {"dataset": dataset, "pulsar": name}
        try:
            for suffix in (".par", ".tim"):
                chosen = unique_best(candidates[suffix], dataset, suffix)
                row[suffix[1:]] = str(chosen.relative_to(work))
                row[suffix[1:] + "_sha256"] = digest(chosen)
                if suffix == ".tim":
                    row["include_sha256"] = timing_includes(chosen, work)
        except ValueError as e:
            row["error"] = str(e).replace(str(work) + "/", "")
        result.append(row)
    return result


def fetch(work, out):
    work.mkdir(parents=True, exist_ok=True)
    manifest = {"archives": {}, "sources": [], "status": "in_progress"}
    try:
        fetch_all(work, out, manifest)
    except Exception as e:
        manifest.update(status="inconclusive", error=f"{type(e).__name__}: {e}")
        raise
    else:
        manifest.update(status="complete", stage="complete")
    finally:
        save(work / "inputs.json", manifest)
        save(out / "inputs.json", manifest)


def fetch_all(work, out, manifest):
    def checkpoint(stage):
        manifest["stage"] = stage
        save(work / "inputs.json", manifest)
        save(out / "inputs.json", manifest)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for source, (record, filename, md5) in RELEASES.items():
            archive = tmp / filename
            url = f"https://zenodo.org/records/{record}/files/{filename}?download=1"
            entry = {"doi": "10.5281/zenodo." + record, "url": url, "md5": md5}
            manifest["archives"][source] = entry
            checkpoint(f"downloading {source}")
            download(url, archive)
            checkpoint(f"checking {source}")
            if digest(archive, "md5") != md5:
                raise ValueError(f"{source} release MD5 mismatch")
            entry.update(sha256=digest(archive), bytes=archive.stat().st_size)
            checkpoint(f"extracting {source}")
            unpack(archive, work / source)
            archive.unlink()
            manifest["sources"] = inventory(work)
            checkpoint(f"inventoried {source}")
    for dataset in ("nb", "wb", "epta"):
        rows = [r for r in manifest["sources"] if r["dataset"] == dataset]
        print(f"{dataset}: {len(rows)} pulsars, {sum('error' not in r for r in rows)} unambiguous pairs")
    # Archive release clocks take precedence for observatories represented there.
    # PINT also needs clocks for EPTA sites and BIPM.
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        url = f"https://codeload.github.com/ipta/pulsar-clock-corrections/tar.gz/{CLOCK_REV}"
        manifest["clocks"] = {"commit": CLOCK_REV, "url": url}
        checkpoint("downloading clock snapshot")
        clock_tar = tmp / "clocks.tar.gz"
        download(url, clock_tar)
        unpack(clock_tar, tmp / "clocksrc")
        clock = work / "standards" / "clock"
        clock.mkdir(parents=True, exist_ok=True)
        for base in [tmp / "clocksrc", work / "ng"]:
            for src in sorted(base.rglob("*")):
                if "clock" in {p.lower() for p in src.parts} and src.suffix.lower() in {".clk", ".dat"}:
                    shutil.copyfile(src, clock / src.name)
        # PINT 1.1.8 opens a global index even when using overrides. Mirror
        # the supported index format locally so its second lookup is pinned too.
        files = sorted(p.name for p in clock.iterdir() if p.name != "index.txt")
        (clock / "index.txt").write_text("".join(f"{f} 365000 ---\n" for f in files))
        manifest["clocks"] = {"commit": CLOCK_REV, "url": url, "snapshot_sha256": digest(clock_tar),
                              "files_sha256": {p.name: digest(p) for p in sorted(clock.iterdir())}}
    ephnames = set()
    for r in manifest["sources"]:
        if "par" in r:
            header = (work / r["par"]).read_text(errors="replace")
            token = re.search(r"^\s*EPHEM\s+(DE\d+)\b", header, re.M | re.I)
            if token:
                ephnames.add(token.group(1).lower())
    manifest["ephemerides"] = {}
    for name in sorted(ephnames | {"de440"}):
        dst = work / "standards" / (name + ".bsp")
        mirrors = ["https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/",
                   "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/",
                   "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/"]
        for base in mirrors:
            url = base + dst.name
            manifest["ephemerides"][name] = {"url": url}
            checkpoint(f"downloading {name}")
            try:
                download(url, dst)
                with dst.open("rb") as stream:
                    if stream.read(8).startswith(b"DAF/SPK"):
                        break
            except subprocess.CalledProcessError:
                continue
        else:
            raise ValueError(f"could not obtain a valid JPL kernel for {name}")
        manifest["ephemerides"][name] = {"url": url, "sha256": digest(dst)}
    checkpoint("checking release coverage")
    for dataset in ("nb", "wb"):
        sources = [r for r in manifest["sources"] if r["dataset"] == dataset]
        if len(sources) != 68:
            raise ValueError(f"expected 68 {dataset} pulsars, found {len(sources)}; check inputs.json")
    if not any(r["dataset"] == "epta" for r in manifest["sources"]):
        raise ValueError("EPTA archive contained no identified timing models")


def ecorr_blocks(model, toas):
    """PINT ECORR has disjoint selection masks and one-second anchor epochs."""
    import numpy as np
    import astropy.units as u
    from pint.models.noise_model import get_ecorr_epochs
    values = np.zeros(len(toas))
    groups = np.full(len(toas), -1)
    seconds = toas.get_tdb_seconds(dtype=np.float64)
    index = 0
    for name in model.params:
        if not re.fullmatch(r"ECORR\d+", name):
            continue
        ecorr = getattr(model, name)
        mask = np.asarray(ecorr.select_toa_mask(toas))
        if mask.dtype == bool:
            mask = np.flatnonzero(mask)
        for epoch in get_ecorr_epochs(seconds[mask]):
            members = mask[epoch]
            if np.any(groups[members] >= 0):
                raise ValueError("overlapping ECORR components require a multi-block covariance")
            groups[members] = index
            values[members] = ecorr.quantity.to_value(u.us)
            index += 1
    return groups, values


def compress(t, y, sigma, freq, backend, design, group, ecorr, dataset):
    """Two inverse-frequency subbands per site/system/epoch, no ECORR mixing."""
    import numpy as np
    if dataset == "wb":
        members = [[i] for i in range(len(t))]
    else:
        keys = defaultdict(list)
        for i, (date, be, g) in enumerate(zip(t, backend, group)):
            keys[int(np.floor(date)), be, g].append(i)
        members = []
        for ids in keys.values():
            ids = np.asarray(ids)
            frequencies = freq[ids]
            median = np.median(frequencies)
            for mask in (frequencies < median, frequencies >= median):
                if mask.any():
                    members.append(ids[mask])
    result = defaultdict(list)
    for ids in members:
        ids = np.asarray(ids)
        w = 1 / sigma[ids]**2
        variance = 1 / w.sum()
        w *= variance
        result["t"].append(w @ t[ids])
        result["y"].append(w @ y[ids])
        result["err"].append(np.sqrt(variance))
        result["freq"].append((w @ (freq[ids] ** -2)) ** -0.5)
        result["backend"].append(backend[ids[0]])
        result["design"].append(w @ design[ids])
        result["ecorr_group"].append(group[ids[0]])
        result["ecorr_us"].append(ecorr[ids[0]])
    order = np.argsort(result["t"], kind="stable")
    return {k: np.asarray(v)[order] for k, v in result.items()}, order


def extract(work, row, out):
    import numpy as np
    import astropy.units as u
    import pint
    from astropy.utils import iers
    from pint.models import get_model_and_toas
    from pint.observatory import global_clock_corrections
    from pint.residuals import Residuals
    from pint.solar_system_ephemerides import load_kernel

    iers.conf.auto_download = False
    clock = work / "standards" / "clock"
    os.environ["PINT_CLOCK_OVERRIDE"] = str(clock)
    local_uri = clock.as_uri() + "/"
    global_clock_corrections.global_clock_correction_url_base = local_uri
    global_clock_corrections.global_clock_correction_url_mirrors = [local_uri]
    par, tim = (work / row[k] for k in ("par", "tim"))
    for source, expected in [(par, row["par_sha256"]), (tim, row["tim_sha256"])] + [
            (work / name, h) for name, h in row["include_sha256"].items()]:
        if digest(source) != expected:
            raise ValueError(f"input changed after acquisition: {source.relative_to(work)}")
    text = par.read_text(errors="replace")
    token = re.search(r"^\s*EPHEM\s+(DE\d+)\b", text, re.M | re.I)
    eph = token.group(1).lower() if token else "de440"
    bipm = "BIPM2021" if row["dataset"] == "epta" else "BIPM2019"
    explicit = re.search(r"^\s*(?:CLK|CLOCK)\s+TT\((BIPM\d{4})\)", text, re.M | re.I)
    if explicit:
        bipm = explicit.group(1).upper()
    kernel = work / "standards" / (eph + ".bsp")
    if not kernel.is_file():
        raise ValueError(f"missing hashed kernel for {eph}")
    # PINT caches the last selected kernel by name. The global Astropy kernel
    # must be reselected even if earlier pulsars used a different ephemeris.
    from astropy.coordinates import solar_system_ephemeris
    solar_system_ephemeris.set(str(kernel))
    load_kernel(eph, path=str(kernel))

    with contextlib.chdir(tim.parent):
        model, toas = get_model_and_toas(
            str(par), str(tim), usepickle=False, allow_T2=True, allow_tcb=True,
            ephem=eph, limits="error", allow_name_mixing=True,
            include_bipm=True, bipm_version=bipm,
        )
    original_name = str(model.PSR.value).upper()
    if PSR.fullmatch(original_name) and ALIASES.get(original_name, original_name) != row["pulsar"]:
        raise ValueError(f"model identifies {original_name}, expected {row['pulsar']}")
    flags = toas.table["flags"]
    good = np.asarray(["cut" not in f for f in flags])
    # Only a window explicitly enabled in the par file excludes TOAs.
    # Plain START/FINISH written by PINT are descriptive fit bounds.
    mjds = np.asarray(toas.get_mjds().value, float)
    for name, predicate in (("START", np.greater_equal), ("FINISH", np.less_equal)):
        param = getattr(model, name, None)
        enabled = re.search(rf"^\s*{name}\s+\S+\s+1\b", text, re.M | re.I)
        if enabled and param is not None and param.value is not None:
            good &= predicate(mjds, float(param.value))
    raw_count, cut_count = len(toas), int((~good).sum())
    toas = toas[good]
    if not len(toas):
        raise ValueError("published quality and timing fit cuts leave zero TOAs")
    residual = Residuals(toas, model)
    y = np.asarray(residual.time_resids.to_value(u.us), float)
    sigma = np.asarray(model.scaled_toa_uncertainty(toas).to_value(u.us), float)
    rms = float(residual.rms_weighted().to_value(u.us))
    tres = model.TRES.quantity
    ref = None if tres is None else float(tres.to_value(u.us))
    if ref is not None and (not np.isfinite(ref) or ref <= 0):
        ref = None
    ratio = None if ref is None else rms / ref
    if not np.isfinite(rms):
        raise ValueError("post-fit residual RMS is nonfinite")
    row.update(raw_toas=raw_count, cut_count=cut_count, used_toas=len(toas),
               rms_postfit_us_before_restoring_dmx=rms, parfile_tres_us=ref,
               wrms_to_tres_ratio=ratio,
               parfile_ntoa=int(model.NTOA.value) if model.NTOA.value is not None else None)
    if ratio is None:
        row["tres_check"] = "not supplied"
    elif not 0.5 <= ratio <= 2.0:
        row["tres_check"] = "excluded"
        raise ValueError(f"pre-DMX-restoration WRMS differs from published TRES by factor {ratio:.4g}")
    elif not 0.8 <= ratio <= 1.25:
        row["tres_check"] = "inspect weighting and model agreement"
    else:
        row["tres_check"] = "within 25 percent by ratio"
    # Gross phase connection problems, e.g. an unsupported EPTA timing model,
    # cannot be admitted as transient events.
    if rms > 0.05 * (1e6 / model.F0.value):
        raise ValueError(f"post-fit timing RMS {rms:g} us above 5% of a rotation")
    dmx = model.components.get("DispersionDMX")
    if dmx is not None:
        y += np.asarray(dmx.DMX_dispersion_delay(toas).to_value(u.us), float)
    for component in model.components.values():
        for name in component.params:
            # DMX is replaced by the chromatic stochastic model downstream.
            # DMJUMP parameters apply to measured wideband DMs, not TOA-only design.
            if name.startswith(("DMX_", "DMJUMP")):
                getattr(model, name).frozen = True
    design, names, _ = model.designmatrix(toas)
    design = np.asarray(design, float)
    scales = np.max(np.abs(design), axis=0)
    keep = scales > 0
    if not np.isfinite(scales).all():
        raise ValueError("nonfinite timing-model Jacobian")
    design = design[:, keep] / scales[keep]
    sky = model.coords_as_ICRS()
    # Use the barycentric arrival epoch for array-coincidence statistics.
    t = np.asarray(model.get_barycentric_toas(toas).to_value(u.day), float)
    freq = np.asarray(model.barycentric_radio_freq(toas).to_value(u.MHz), float)
    if not np.all(np.isfinite([t, y, sigma, freq])) or np.any((sigma <= 0) | (freq <= 0)):
        raise ValueError("nonfinite TOA, frequency or uncertainty")
    flags = toas.table["flags"]
    backend = np.asarray([str(ob) + "/" + f.get("f", f.get("sys", f.get("be", "unknown")))
                          for f, ob in zip(flags, toas.table["obs"])])
    groups, ecorr = ecorr_blocks(model, toas)
    result, order = compress(t, y, sigma, freq, backend, design, groups, ecorr, row["dataset"])
    result.update(ra=sky.ra.rad, dec=sky.dec.rad, design_names=np.array(names)[keep],
                  dm_restored=dmx is not None, dm_retained=True, f0_hz=float(model.F0.value))
    if row["dataset"] == "wb":
        result["pp_dm"] = np.array([float(f["pp_dm"]) for f in flags])[order]
        result["pp_dme"] = np.asarray(model.scaled_dm_uncertainty(toas).value, float)[order]
        # Remove fitted receiver jumps and secular DM while restoring DMX.
        from pint.residuals import WidebandDMResiduals
        dm_res = WidebandDMResiduals(toas, model, subtract_mean=False)
        offset = dmx.dmx_dm(toas).value if dmx is not None else 0.0
        result["dm"] = np.asarray(dm_res.resids.value + offset, float)[order]
        result["dm_err"] = result["pp_dme"]
        if (not np.isfinite(result["dm"]).all() or not np.isfinite(result["dm_err"]).all()
                or (result["dm_err"] <= 0).any()):
            raise ValueError("unusable paired wideband DM or uncertainty")
    dst = out / row["dataset"] / (row["pulsar"] + ".npz")
    dst.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".partial", dir=dst.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            np.savez_compressed(stream, **result)
            stream.flush()
            temporary.replace(dst)
        finally:
            temporary.unlink(missing_ok=True)
    return {"status": "ok", "compressed_rows": len(result["t"]),
            "min_mjd_tdb": float(t.min()), "max_mjd_tdb": float(t.max()),
            "spin_hz": float(model.F0.value),
            "ephemeris": eph, "clock": bipm, "parfile_clock": str(model.CLOCK.value),
            "dmx_restored": dmx is not None, "pint": pint.__version__,
            "components": sorted(model.components), "source_pulsar": str(model.PSR.value),
            "toa_track_mode": str(residual.track_mode)}


def reduce(work, out, dataset, shard, shards):
    from loguru import logger
    from importlib.metadata import version
    from pint.config import runtimefile
    import astropy_iers_data
    logger.remove()
    import sys
    logger.add(sys.stderr, level="WARNING")
    if not 0 <= shard < shards:
        raise ValueError("expected 0 <= shard < shards")
    runtime = Path(runtimefile("."))
    receipt = {
        "packages": {p: version(p) for p in ("pint-pulsar", "astropy", "astropy-iers-data", "pyerfa")},
        "pint_runtime_sha256": {str(p.relative_to(runtime)): digest(p)
                                for p in sorted(runtime.rglob("*")) if p.is_file()},
        "iers_sha256": {p.name: digest(p) for p in
                        sorted((Path(astropy_iers_data.__file__).parent / "data").iterdir())
                        if p.is_file()},
    }
    save(out / f"runtime-{dataset}-{shard}.json", receipt)
    rows = json.loads((work / "inputs.json").read_text())["sources"]
    ng = {r["pulsar"] for r in rows if r["dataset"] == "nb"}
    selected = [r for r in rows if r["dataset"] == dataset and r["pulsar"] in ng]
    log = []
    for row in selected[shard::shards]:
        r = dict(row)
        try:
            if "error" in r:
                raise ValueError(r["error"])
            r.update(extract(work, r, out))
        except Exception as e:
            (out / dataset / (row["pulsar"] + ".npz")).unlink(missing_ok=True)
            r.update(status="inconclusive", error=f"{type(e).__name__}: {e}")
            traceback.print_exc()
        log.append(r)
        save(out / f"coverage-{dataset}-{shard}.json", log)
        print(dataset, r["pulsar"], r["status"], flush=True)
    save(out / f"coverage-{dataset}-{shard}.json", log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["fetch", "reduce"])
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--dataset", choices=["nb", "wb", "epta"])
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    args = parser.parse_args()
    if args.command == "fetch":
        fetch(args.work.resolve(), args.out.resolve())
    else:
        if not args.dataset:
            parser.error("reduce requires --dataset")
        reduce(args.work.resolve(), args.out.resolve(), args.dataset, args.shard, args.shards)
