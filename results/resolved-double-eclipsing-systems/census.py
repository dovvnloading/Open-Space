"""Join measured Gaia eclipse periods to spatially resolved wide pairs."""
import argparse
import csv
import gzip

import numpy as np
import pandas as pd


def read_header(stream):
    header = {}
    while True:
        block = stream.read(2880).decode("ascii")
        for index in range(0, len(block), 80):
            card = block[index:index + 80]
            if card.startswith("END "):
                return header
            if card[8:10] == "= ":
                header[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
        if not block:
            raise ValueError("Incomplete FITS header")


def census(wide, eclipses):
    eb = pd.read_csv(eclipses).set_index("source_id")
    if len(eb) != 2184477 or not eb.index.is_unique:
        raise ValueError("Incomplete Gaia eclipsing binary table")
    joined = []
    with gzip.open(wide, "rb") as stream:
        read_header(stream)
        h = read_header(stream)
        formats = {"E": ">f4", "D": ">f8", "K": ">i8", "4A": "S4"}
        dt = np.dtype([(h[f"TTYPE{i}"], formats[h[f"TFORM{i}"]]) for i in range(1, int(h["TFIELDS"]) + 1)])
        assert dt.itemsize == int(h["NAXIS1"])
        remaining = int(h["NAXIS2"])
        while remaining:
            size = min(50000, remaining)
            block = np.frombuffer(stream.read(size * dt.itemsize), dtype=dt, count=size)
            remaining -= size
            primary = eb.index.get_indexer(block["source_id1"]) >= 0
            secondary = eb.index.get_indexer(block["source_id2"]) >= 0
            for row in block[primary & secondary]:
                item = {k: row[k].item() for k in dt.names}
                item["binary_type"] = item["binary_type"].decode()
                for i in (1, 2):
                    s = eb.loc[item[f"source_id{i}"]]
                    item[f"period{i}"] = 1 / s["frequency"]
                    item[f"rank{i}"] = s["global_ranking"]
                item["sep_arcsec"] = 3600 * item["pairdistance"]
                joined.append(item)
    return joined


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("wide")
    p.add_argument("eclipses")
    p.add_argument("output")
    a = p.parse_args()
    result = census(a.wide, a.eclipses)
    with open(a.output,"w",newline="") as f:
        writer = csv.DictWriter(f, fieldnames=result[0])
        writer.writeheader()
        writer.writerows(result)
    print(f"{len(result)} doubly matched wide pairs")
