"""Check complete Fe XIV parent lineage using the deposited FITS headers."""

from argparse import ArgumentParser
import csv
from datetime import datetime
import gzip
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

from astropy.io import fits


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    receipts = json.loads(gzip.decompress((args.packet / "header-receipts.json.gz").read_bytes()))
    receipt_by_line = {row["fe"]: row for row in receipts}
    result = []
    expected_headers = {f"{name}.hdr" for name in receipt_by_line}
    with ZipFile(args.packet / "fe-headers.zip") as archive:
        if set(archive.namelist()) != expected_headers:
            raise ValueError("Header files and receipts differ")
        for orbit in (795, 772):
            products = json.loads(gzip.decompress((args.packet / f"orbit{orbit}.json.gz").read_bytes()))
            donors = {}
            for row in products:
                if row["PROD_ID"] == "Total brightness" and row["QUALITY"] == 0:
                    parents = frozenset(row[f"PARENT{i}"] for i in range(1, 4))
                    if row["N_PARENT"] != 3 or None in parents or len(parents) != 3:
                        raise ValueError(f"Incomplete wideband metadata: {row['name']}")
                    if parents in donors:
                        raise ValueError("Duplicate wideband parent set")
                    donors[parents] = row

            for row in products:
                if row["PROD_ID"] != "Green line":
                    continue
                name = row["name"]
                receipt = receipt_by_line.pop(name)
                raw = archive.read(f"{name}.hdr")
                if hashlib.sha256(raw).hexdigest() != receipt["header_sha256"]:
                    raise ValueError(f"Header digest differs: {name}")
                header = fits.Header.fromstring(raw.decode("ascii"), sep="")
                if (header["FILENAME"] != name or header["ORBIT_ID"] != orbit
                        or header["QUALITY"] != row["QUALITY"]
                        or header["CYCLE_ID"] != row["CYCLE_ID"]
                        or datetime.fromisoformat(header["DATE-AVG"]) != datetime.fromisoformat(row["DATE-AVG"])):
                    raise ValueError(f"Header and archive disagree for {name}")
                parents = [header[f"PARENT{i}"] for i in range(1, header["N_PARENT"] + 1)]
                donor = donors[frozenset(p for p in parents if p.startswith("aspiics_wb_l2_"))]
                lag = (datetime.fromisoformat(donor["DATE-AVG"])
                       - datetime.fromisoformat(header["DATE-AVG"])).total_seconds()
                different_cycle = header["CYCLE_ID"] != donor["CYCLE_ID"]
                if [(donor["name"], lag, not different_cycle)] != [
                    (match["bt"], match["lag_s_bt_minus_fe"], match["same_cycle"])
                    for match in receipt["matches"]
                ]:
                    raise ValueError(f"Independent lineage receipt disagrees for {name}")
                result.append({
                    "orbit": orbit, "line": name, "quality": header["QUALITY"],
                    "n_parent": header["N_PARENT"], "donor": donor["name"],
                    "different_cycle": int(different_cycle), "lag_seconds": lag,
                })
    if receipt_by_line:
        raise ValueError("Some line headers were not checked")
    with args.output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)


if __name__ == "__main__":
    main()
