"""Reduce the original photometry to replayable input columns."""

from pathlib import Path
import argparse
import re

import pandas as pd

ROOT = Path(__file__).resolve().parent


def prepare(source):
    data = []
    for line in (source / "asascgi.txt").read_text().splitlines():
        if line.startswith("#dataset="):
            field = line.split()[-1]
        if re.match(r"\s*\d", line):
            row = line.split()
            data.append(
                [
                    "ASAS-3",
                    field,
                    float(row[0]) + 2450000,
                    float(row[1]),
                    float(row[6]),
                    row[-2],
                ]
            )
    d = pd.read_parquet(source / "patrollc.parquet")
    assert set(d.asas_sn_id) == {515396329465}
    for r in d.itertuples():
        data.append(["ASAS-SN", r.camera, r.jd, r.mag, r.mag_err, r.quality])
    pd.DataFrame(
        data, columns=["survey", "channel", "jd", "mag", "err", "flag"]
    ).to_csv(
        ROOT / "ground-photometry.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.read_csv(source / "asphoto_csv")[["hjd", "camera"]].to_csv(
        ROOT / "legacy-clock.csv", index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    prepare(parser.parse_args().source)
