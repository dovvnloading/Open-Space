"""Derive a compact, lossless selection table from the JPL JSON response."""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path

COLUMNS = ('pdes', 'epoch', 'e', 'q', 'i', 'om', 'w', 'data_arc', 'condition_code', 'H')


def pack(raw):
    source = json.loads(raw)
    indices = [source['fields'].index(column) for column in COLUMNS]
    text = io.StringIO(newline='')
    writer = csv.writer(text, lineterminator='\n')
    writer.writerow(COLUMNS)
    for row in source['data']:
        writer.writerow([row[index] for index in indices])
    return gzip.compress(text.getvalue().encode(), compresslevel=9, mtime=0)


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('snapshot', type=Path)
    cli.add_argument('output', type=Path)
    args = cli.parse_args()
    args.output.write_bytes(pack(args.snapshot.read_bytes()))
