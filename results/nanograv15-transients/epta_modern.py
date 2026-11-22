"""Try the released EPTA DR2new timing pair only where DR2full could not be validated."""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import traceback

from prepare import ALIASES, PSR, digest, extract, save, timing_includes, unique_best


def published_pairs(work, failures):
    pairs = {name: {".par": [], ".tim": []} for name in failures}
    root = work / "epta"
    folders = set()
    for f in root.rglob("*"):
        if not f.is_file() or f.suffix.lower() not in (".par", ".tim"):
            continue
        parts = f.relative_to(root).parts
        directories = [part.casefold() for part in parts[:-1]]
        # Require the collaboration's named selection, never DR2new+ or a
        # legacy-plus-InPTA set. Aggregate TOA filenames identify a pulsar.
        if "dr2new" not in directories or any("inpta" in part or part.endswith("+")
                                               for part in directories):
            continue
        release = directories.index("dr2new")
        folders.add(str(root.joinpath(*parts[:release+1]).relative_to(work)))
        names = PSR.findall(f.name)
        if len(names) != 1:
            continue
        name = ALIASES.get(names[0].upper(), names[0].upper())
        if name not in pairs:
            continue
        pairs[name][f.suffix.lower()].append(f)
    return sorted(folders), pairs


def select_inputs(work, row, choices):
    row["candidate_paths"] = {k: [str(p.relative_to(work)) for p in sorted(values)]
                              for k, values in choices.items()}
    if not choices[".par"] or not choices[".tim"]:
        raise ValueError("no paired named model and aggregate timing file in DR2new")
    for ext in (".par", ".tim"):
        path = unique_best(choices[ext], "epta", ext)
        key = ext[1:]
        row[key] = str(path.relative_to(work))
        row[key + "_sha256"] = digest(path)
        if key == "tim":
            # An official DR2new aggregate may reference files shared with the
            # full release. Follow only its explicit INCLUDEs inside EPTA.
            release_root = work / "epta"
            included = timing_includes(path, release_root)
            row["include_sha256"] = {str(Path("epta") / name): h for name, h in included.items()}


def run(work, out, shard):
    ledger = out / f"coverage-epta-{shard}.json"
    previous = json.loads(ledger.read_text())
    original = [r.get("dr2full_attempt", r) for r in previous]
    failed = [r["pulsar"] for r in original if r["status"] != "ok"]
    roots, pairs = published_pairs(work, failed) if failed else ([], {})
    updated = []
    for primary in original:
        name = primary["pulsar"]
        full = dict(primary, release_selection="DR2full")
        dst = out / "epta" / f"{name}.npz"
        if full["status"] == "ok":
            if not dst.is_file():
                raise FileNotFoundError(f"validated DR2full product missing for {name}")
            updated.append(full)
            continue
        candidate = {"pulsar":name, "dataset":"epta", "release_selection":"DR2new",
                     "archive_directories":roots, "release_scope":
                     "Publisher-selected modern EPTA backends; checks apply only within this subset's observed span."}
        try:
            select_inputs(work, candidate, pairs[name])
            with TemporaryDirectory(dir=out) as tmp:
                stage = Path(tmp)
                candidate.update(extract(work, candidate, stage))
                product = stage / "epta" / dst.name
                import numpy as np
                with np.load(product, allow_pickle=False) as timing:
                    nights = np.unique(np.floor(timing["t"]))
                candidate["distinct_observing_days"] = len(nights)
                candidate["median_observing_day_separation_days"] = (
                    float(np.median(np.diff(nights))) if len(nights) > 1 else None)
                dst.parent.mkdir(parents=True, exist_ok=True)
                product.replace(dst)
            selected = dict(candidate, dr2full_attempt=full)
        except Exception as exc:
            dst.unlink(missing_ok=True)
            candidate.update(status="inconclusive",
                             error=f"{type(exc).__name__}: {exc}".replace(str(work)+"/", ""))
            traceback.print_exc()
            selected = dict(full, dr2full_attempt=full, dr2new_attempt=candidate)
        updated.append(selected)
        print("EPTA fallback", name, candidate["status"], candidate.get("error",""), flush=True)
    save(ledger, updated)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--shard", required=True, type=int)
    args = parser.parse_args()
    run(args.work.resolve(), args.out.resolve(), args.shard)
