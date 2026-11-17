"""Recover failed positional ZTF DR24 API queries from IRSA's public HATS collection."""
import csv
import datetime
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

from astropy import units as u
from astropy_healpix import HEALPix
from pyarrow import fs
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "wise"))
from fetch import position, separation, targets

BUCKET = "ipac-irsa-ztf"
CAT = f"{BUCKET}/ztf/enhanced/dr24/lc/hats/ztf_dr24_lc-hats"
COLUMNS = ["objectid", "objra", "objdec", "filterid", "nepochs", "lightcurve"]


def cone(order, ra, dec, radius):
    grid = HEALPix(nside=2**order, order="nested")
    return set(grid.cone_search_lonlat(ra*u.deg, dec*u.deg, radius*u.arcsec).tolist())


def retrieve(t, remote, partitions):
    ra, dec = position(t, 2022)
    included = {}
    for order in sorted({key[0] for key in partitions}):
        for pix in cone(order, ra, dec, 2.5):
            if (order, pix) in partitions:
                leaf = f"{CAT}/dataset/Norder={order}/Dir={pix//10000*10000}/Npix={pix}.parquet"
                included[leaf] = [order, pix]
    narrow = sorted(cone(16, ra, dec, 2.5))
    width = 4**(29-16)
    bounds = [[("filterid", "=", 2), ("_healpix_29", ">=", p*width),
               ("_healpix_29", "<", (p+1)*width)] for p in narrow]
    # Select a small interval in the physically HEALPix-sorted Parquet data before reading arrays.
    points, objects, parts = [], [], []
    for path in included:
        table = pq.read_table(path, filesystem=remote, filters=bounds, columns=COLUMNS)
        parts.append(dict(path="s3://"+path, size=remote.get_file_info(path).size))
        for row in table.to_pylist():
            off = separation(row["objra"], row["objdec"], ra, dec)
            if off > 2.5:
                continue
            objects.append(dict(oid=row["objectid"], offset_arcsec=off, nepochs=row["nepochs"]))
            if off > 1.5:
                continue
            series = row["lightcurve"]
            if isinstance(series, dict):
                series = [dict(zip(series, values)) for values in zip(*series.values())]
            for pt in series:
                points.append(dict(oid=row["objectid"], mjd=pt["hmjd"], mag=pt["mag"], magerr=pt["magerr"],
                                   catflags=pt["catflags"], clrcoeff=pt["clrcoeff"],
                                   objra=row["objra"], objdec=row["objdec"]))
    dest = ROOT / f"hats_{t['label']}.csv.gz"
    out = io.StringIO()
    writer = csv.DictWriter(out, "oid mjd mag magerr catflags clrcoeff objra objdec".split(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(points, key=lambda p: (p["oid"], p["mjd"])))
    data = out.getvalue().encode()
    dest.write_bytes(gzip.compress(data, mtime=0))
    return dict(label=t["label"], collection="s3://"+CAT, query_epoch=2022, center_deg=[ra, dec],
                radius_arcsec=2.5, acceptance_arcsec=1.5, filterid=2, parts=parts, objects=objects,
                time_scale="HMJD (column mjd)", rows=len(points), sha256=hashlib.sha256(data).hexdigest())


def main():
    failed = {r["label"] for r in json.loads((ROOT / "query_receipts.json").read_text()) if "error" in r}
    if not failed:
        return
    remote = fs.S3FileSystem(anonymous=True, region=fs.resolve_s3_region(BUCKET))
    index_path = f"{CAT}/partition_info.csv"
    with remote.open_input_stream(index_path) as stream:
        index = stream.read()
    leaves = {(int(row["Norder"]), int(row["Npix"])) for row in csv.DictReader(io.StringIO(index.decode()))}
    metadata = dict(partition_index="s3://"+index_path, index_sha256=hashlib.sha256(index).hexdigest(),
                    utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), targets=[])
    for target in targets():
        if target["label"] not in failed:
            continue
        row = retrieve(target, remote, leaves)
        print(target["label"], row["rows"], row["objects"], flush=True)
        metadata["targets"].append(row)
    (ROOT / "hats_receipts.json").write_text(json.dumps(metadata, indent=2)+"\n")


if __name__ == "__main__":
    main()
