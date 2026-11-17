"""Record URLs and hashes for exactly the inputs used by the archived reductions."""

from pathlib import Path
from query import HERE, CACHE
import json, hashlib, gzip, csv


def main():
    optical = json.loads((HERE / "selected.json").read_text())
    ztf = json.loads((HERE / "ztf_results.json").read_text())
    items = []
    extra = []
    for name, rows in optical.items():
        items.extend([f"ls_{name}", f"stationary_{name}"])
        extra.extend([name, "masks_" + name])
        for r in rows:
            items.extend([r["md5sum"], r["md5sum"] + "_mask"])
    for key in [
        "ztf_Chen1_2025",
        "ztf_Chen1_2026",
        "ztf_Chen2_2025",
        "ztf_Chen2_2026",
        "ztf_RR_2025",
        "ztf_RR_2026",
    ]:
        extra.append(key)
    for res in ztf.values():
        for v in res["visits"]:
            items.extend(["ztf_" + v["pid"] + "_" + f for f in ["diff", "mask", "psf"]])
    items.extend(["static_Chen2_412", "Chen1_stationary_primary"])
    receipts = []
    for name in dict.fromkeys(items + extra):
        rec = json.loads((CACHE / (name + ".json")).read_text())
        raw = CACHE / (name + ".raw")
        assert hashlib.sha256(raw.read_bytes()).hexdigest() == rec["sha256"]
        receipts.append(dict(label=name, required_for_reduction=name in items, **rec))
    (HERE / "inputs.json").write_text(json.dumps(receipts, indent=2) + "\n")


if __name__ == "__main__":
    main()
