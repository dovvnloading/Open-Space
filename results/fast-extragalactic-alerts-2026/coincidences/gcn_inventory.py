"""Rebuild the frozen list of GCN circulars; optionally pass a local GCN archive path."""
import csv
from datetime import UTC, datetime
import json
from pathlib import Path
import sys
import tarfile
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
START = datetime(2026, 11, 1, tzinfo=UTC)
END = datetime(2026, 11, 17, tzinfo=UTC)


def read_archive(archive):
    index = {}
    for member in archive:
        if member.isfile() and member.name.endswith('.json'):
            obj = json.load(archive.extractfile(member))
            dt = datetime.fromtimestamp(obj['createdOn'] / 1000, UTC)
            if START <= dt < END:
                index[obj['circularId']] = (dt.isoformat(),obj.get('eventId'),obj['subject'])
    return index


if __name__ == '__main__':
    if len(sys.argv)>1:
        with tarfile.open(sys.argv[1]) as arc:
            index = read_archive(arc)
    else:
        with urlopen('https://gcn.nasa.gov/circulars/archive.json.tar.gz',timeout=90) as raw:
            with tarfile.open(fileobj=raw,mode='r|gz') as arc:
                index = read_archive(arc)
    for n in range(46123, 46310):
        if n not in index:
            with urlopen(f'https://gcn.nasa.gov/circulars/{n}.json',timeout=30) as f:
                obj = json.load(f)
            dt = datetime.fromtimestamp(obj['createdOn']/1000, UTC)
            if not START <= dt < END:
                raise ValueError(f'Unexpected publication time for circular {n}')
            index[n] = (dt.isoformat(),obj.get('eventId'),obj['subject'])
    with (HERE / 'gcn_index.csv').open('w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['circular','issued_utc','event','subject'])
        w.writerows((n, *x) for n,x in sorted(index.items()))
