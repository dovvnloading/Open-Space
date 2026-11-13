"""Package one or more completed portable queue exports with a verifiable public index."""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


def sha(raw):
    return sha256(raw).hexdigest()


def package(exports, output):
    output.mkdir(parents=True, exist_ok=True)
    archive = output / "discovery-evidence.zip"
    index = output / "discovery-registry.csv"
    report = {
        "schema_version": 1,
        "claim_status": "screening and follow-up measurements",
        "campaigns": [],
    }
    with ZipFile(archive, "w") as bundle, index.open(
        "w", encoding="utf-8", newline=""
    ) as listing:
        table = None
        for name, source in exports:
            registry = json.loads((source / "registry.json").read_bytes())
            if set(registry["counts"]["tasks"]) != {"done"}:
                raise ValueError("Release requires a completed export: " + name)
            inputs = {
                *registry["campaign"]["tess_plan_sha256"],
                *registry["campaign"]["wise_plan_sha256"],
                *registry["campaign"]["ffi_plan_sha256"],
                registry["campaign"]["catalog_sha256"],
            }
            hashes = {t["evidence"]: t["evidence_sha256"] for t in registry["tasks"]}
            hashes.update({f"inputs/{digest}.json": digest for digest in inputs})
            for member in sorted({"registry.csv", "registry.json", *hashes}):
                raw = (source / member).read_bytes()
                if member in hashes and sha(raw) != hashes[member]:
                    raise ValueError(f"Inconsistent evidence in {name}/{member}")
                file = ZipInfo(f"{name}/{member}", date_time=(1980, 1, 1, 0, 0, 0))
                file.compress_type = ZIP_DEFLATED
                file.external_attr = 0o100644 << 16
                bundle.writestr(file, raw)
            with (source / "registry.csv").open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                if table is None:
                    table = csv.DictWriter(
                        listing, ["campaign", *reader.fieldnames], lineterminator="\n"
                    )
                    table.writeheader()
                for row in reader:
                    table.writerow(
                        {
                            "campaign": name,
                            **row,
                            "evidence": name + "/" + row["evidence"],
                        }
                    )
            report["campaigns"].append(
                {
                    "name": name,
                    "project": registry["campaign"],
                    "counts": registry["counts"],
                    "measurement_rows": len(registry["findings"]),
                    "registry_json_sha256": sha(
                        (source / "registry.json").read_bytes()
                    ),
                }
            )
    report.update(
        bundle_sha256=sha(archive.read_bytes()), csv_sha256=sha(index.read_bytes())
    )
    (output / "discovery-run.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--export",
        nargs=2,
        action="append",
        required=True,
        metavar=("NAME", "DIRECTORY"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    names = [n for n, _ in args.export]
    if len(names) != len(set(names)) or any(
        not n.replace("-", "").isalnum() for n in names
    ):
        parser.error("Use distinct alphanumeric or hyphenated archive names")
    package([(name, Path(path)) for name, path in args.export], args.output)


if __name__ == "__main__":
    main()
