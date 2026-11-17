"""Freeze all assigned DR1 exposures, stellar classifications and Gaia geometric distances."""
import argparse
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from collections import defaultdict
from urllib.parse import urlencode
from urllib.request import urlopen
from time import sleep

DESI = "https://datalab.noirlab.edu/tap/sync"
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"


def freeze(out, name, endpoint, query, columns):
    """Archive a successful CSV response alongside its exact query. Reject TAP error bodies."""
    path = out / (name + ".csv")
    receipt = out / (name + ".json")
    url = endpoint + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                                           MAXREC=100000, QUERY=query))
    for attempt in range(3):
        try:
            if path.exists() and receipt.exists():
                raw = path.read_bytes()
                old = json.loads(receipt.read_text())
                if old["url"] != url or old["sha256"] != hashlib.sha256(raw).hexdigest():
                    raise ValueError("query or content changed; use a new output directory")
            else:
                # TAP accepts either GET or form-encoded POST for synchronous ADQL.
                # Long TARGETID lists can exceed a service's GET request-line limit.
                post = len(url) > 3800
                body = url.partition("?")[2].encode() if post else None
                with urlopen(endpoint if post else url, data=body, timeout=100) as response:
                    raw = response.read()
            reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
            if reader.fieldnames != columns:
                raise RuntimeError(f"{name}: unexpected TAP response {raw[:120]!r}")
            rows = list(reader)
            if len(rows) >= 100000:
                raise RuntimeError("possible TAP truncation")
            path.write_bytes(raw)
            receipt.write_text(json.dumps(dict(url=url, adql=query, method="POST" if len(url)>3800 else "GET", rows=len(rows),
                                                sha256=hashlib.sha256(raw).hexdigest()), indent=2) + "\n")
            return rows
        except (OSError, RuntimeError):
            if attempt == 2:
                raise
            sleep(2 ** attempt)


def has_gaia(row):
    # DR1 mws uses 999999 as a fill value for missing Gaia matches.
    value = row["source_id"]
    return bool(value) and int(value) > 0 and int(value) != 999999


def chunks(ids, size=55):
    ids = sorted(set(int(key) for key in ids))
    return [",".join(str(x) for x in ids[i:i+size]) for i in range(0, len(ids), size)]


def gather(out, prefix, endpoint, ids, columns, query, size=55):
    def part(pair):
        index, group = pair
        return freeze(out, f"{prefix}_{index:03d}", endpoint, query.format(ids=group), columns.split(","))
    with ThreadPoolExecutor(4) as pool:
        return [row for batch in pool.map(part, enumerate(chunks(ids, size))) for row in batch]


def write_rows(path, rows, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def audit(target_file, out, size=55):
    targets = list(csv.DictReader(target_file.open()))
    ids = [r["targetid"] for r in targets]
    out.mkdir(parents=True, exist_ok=True)
    cols = "targetid,source_id,survey,program,rr_spectype"
    mws = gather(out, "mws", DESI, ids, cols, f"SELECT {cols} FROM desi_dr1.mws WHERE targetid IN ({{ids}})", size)
    source_ids = [r["source_id"] for r in mws if has_gaia(r)]
    # A single star may have several targeting identifiers or survey groups.
    aliases = gather(out, "aliases", DESI, source_ids, cols,
                     f"SELECT {cols} FROM desi_dr1.mws WHERE source_id IN ({{ids}})", size)
    all_ids = set(ids) | {r["targetid"] for r in aliases}
    cols = "targetid,tileid,expid,night,survey,program,fiber,petal_loc"
    visits = gather(out, "visits", DESI, all_ids, cols,
                    "SELECT f.targetid,f.tileid,e.expid,e.night,e.survey,e.program,f.fiber,f.petal_loc "
                    "FROM desi_dr1.fiberassign AS f JOIN desi_dr1.exposure AS e ON f.tileid=e.tileid "
                    "WHERE f.targetid IN ({ids})", size)
    cols = "source_id,r_med_geo,r_lo_geo,r_hi_geo"
    ds = gather(out, "distances", GAIA, source_ids, cols,
                f"SELECT {cols} FROM external.gaiaedr3_distance WHERE source_id IN ({{ids}})", size)
    distances = {r["source_id"]:r for r in ds}
    star_sources, source_targets, star_types, star_visits = (defaultdict(set) for _ in range(4))
    for row in mws:
        star_types[row["targetid"]].add(row["rr_spectype"])
        if has_gaia(row):
            star_sources[row["targetid"]].add(row["source_id"])
    for row in aliases:
        source_targets[row["source_id"]].add(row["targetid"])
        star_types[row["targetid"]].add(row["rr_spectype"])
    for row in visits:
        star_visits[row["targetid"]].add(int(row["expid"]))
    summary = []
    for t in targets:
        target = t["targetid"]
        gaia = star_sources[target]
        same = {target}.union(*(source_targets[s] for s in gaia))
        gids = sorted(gaia)
        src = gids[0] if len(gids) == 1 else ""
        expected = set().union(*(star_visits[t] for t in same))
        available = {int(x) for x in t["expids"].split("|")}
        spec = set().union(*(star_types[t] for t in same))
        summary.append(dict(targetid=target, source_id=src, alias_targetids="|".join(sorted(same)),
                            requested_expids=t["expids"], expected_expids="|".join(map(str, sorted(expected))),
                            all_assigned=(expected == available and len(expected) >= 2),
                            star_class_consistent=(spec == {"STAR"}), ambiguous_gaia=len(gaia)>1,
                            **{k:distances.get(src, {}).get(k,"") for k in ("r_med_geo","r_lo_geo","r_hi_geo")}))
    write_rows(out / "audit.csv", summary, list(summary[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--chunk-size", type=int, default=55)
    args = parser.parse_args()
    audit(args.targets, args.output, args.chunk_size)
