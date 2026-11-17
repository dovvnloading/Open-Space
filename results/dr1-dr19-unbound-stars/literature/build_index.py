"""Rebuild the published-candidate index from the sources in manifest.tsv."""
import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import zipfile
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
FIELDS = ["catalogue", "source_id", "release", "name", "ra_deg", "dec_deg", "epoch", "pmra", "pmdec"]
ID = re.compile(r"(?<!\d)\d{16,19}(?!\d)")


def vizier_tables(data):
    """The ASU TSV format allows consecutive tables with different column names."""
    lines = data.decode().splitlines()
    for i, line in enumerate(lines):
        if not re.fullmatch(r"[-\t]+", line) or "\t" not in line:
            continue
        header = lines[i-2].split("\t")
        for entry in lines[i+1:]:
            if not entry.strip() or entry.startswith("#"):
                break
            yield {k:v.strip() for k,v in zip(header, entry.split("\t"))}


def first(row, keys):
    return next((row[k] for k in keys if row.get(k)), "")


def parse_vizier(key, release, data):
    for item in vizier_tables(data):
        source = first(item, ["GaiaDR3", "GaiaDR2", "Gaia"])
        source = source.split()[-1] if source else ""
        name = first(item, ["Name", "Target", "HVS"])
        ra = first(item, ["RA_ICRS", "RAICRS", "RAJ2000", "_RA.icrs", "_RA"])
        dec = first(item, ["DE_ICRS", "DEICRS", "DEJ2000", "_DE.icrs", "_DE"])
        epoch = {"DR2": "2015.5", "DR3": "2016.0"}[release]
        if key in {"harvardli23","mar19"}:
            # Coordinate columns marked J2000; exact-ID checks remain primary.
            epoch = ""
        if not (source.isdecimal() or (ra and dec)):
            continue
        yield [key, source, release, name, ra, dec, epoch, item.get("pmRA", ""), item.get("pmDE", "")]


def sexagesimal(value, hours=False):
    parts = re.split(r"[:\s]+", value.strip())
    sign = -1 if parts[0].startswith("-") else 1
    d, m, s = map(float, parts) if len(parts)==3 else (float(parts[0]), 0, 0)
    return f"{sign*(abs(d)+m/60+s/3600)*(15 if hours else 1):.9f}"


def parse_ofsc(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for fname in z.namelist():
            if not fname.endswith(".json"):
                continue
            obj = next(iter(json.loads(z.read(fname)).values()))
            ras, decs = obj["ra"], obj["dec"]
            # Use a coordinate pair sharing a publication source whenever possible.
            pairs = [(x,y) for x in ras for y in decs if set(x["source"].split(",")) & set(y["source"].split(","))]
            x,y = next(iter(pairs), (ras[0],decs[0]))
            ra = sexagesimal(x["value"], hours=x.get("u_value") in ("hours","hms"))
            dec = sexagesimal(y["value"])
            aliases = set()
            for a in obj["alias"]:
                m = re.fullmatch(r"Gaia\s*(DR[123])\s+([\d\s]+)", a["value"])
                if m:
                    aliases.add((m[1], "".join(m[2].split())))
            for release,source in sorted(aliases or {("", "")}):
                yield ["ofsc_pinned", source, release, obj["name"], ra, dec, "", "", ""]


def parse_pdf(key, release, path):
    doc = subprocess.check_output(["pdftotext","-layout",str(path),"-"],stderr=subprocess.DEVNULL).decode()
    if key == "lam25pdf":
        coords = []
        for line in doc.splitlines():
            if "Position (" in line and "ICRS(Epoch J2016.0)" in line:
                coords.extend(re.findall(r"\((-?\d+\.\d+),\s*(-?\d+\.\d+)\)", line))
        if len(coords) != 10:
            raise ValueError("Guo et al. Tables 3–5 changed")
        for idx, (ra, dec) in enumerate(coords, 5):
            yield [key, "", "", f"LAMOST-HVS{idx}", ra, dec, "2016.0", "", ""]
        return
    if key=="bhat":
        ids = re.findall(r"^\s*\d+\s+(\d{16,19})\s+\d\.\d",doc,re.M)
        if len(set(ids))!=53:
            raise ValueError("Bhat published candidate-table extraction changed")
    elif key=="desihvs1pdf":
        ids = re.findall(r"Gaia DR3 source id\s+(\d{16,19})", doc)
    else:
        ids = ID.findall(doc)
    names={"desi312pdf":"DESI-312","desihvs1pdf":"DESI-HVS1"}
    for sid in sorted(set(ids)):
        yield [key,sid,release,names.get(key,""),"","","","",""]


def extract(spec, path):
    key, release = spec["catalogue"], spec["release"]
    if spec["format"]=="vizier":
        return parse_vizier(key,release,path.read_bytes())
    if spec["format"]=="ofsc":
        return parse_ofsc(path.read_bytes())
    if spec["format"]=="pdf":
        return parse_pdf(key,release,path)
    with path.open() as f:
        return [[key,r["GaiaDR3_source_id"],"DR3","","","","","",""] for r in csv.DictReader(f)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--output",type=Path,default=HERE/"index.csv")
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    rows=set()
    with (HERE/"manifest.tsv").open() as manifest:
        for spec in csv.DictReader(manifest,delimiter="\t"):
            path=args.cache/spec["catalogue"]
            if not path.exists():
                path.write_bytes(urlopen(spec["url"],timeout=120).read())
            raw=path.read_bytes()
            # VizieR's response headers include a varying query timestamp.
            scientific_bytes=(b"\n".join(l for l in raw.splitlines() if l and not l.startswith(b"#")) if spec["format"]=="vizier" else raw)
            if hashlib.sha256(scientific_bytes).hexdigest()!=spec["scientific_sha256"]:
                raise ValueError(f"Source content changed: {path}")
            projected=set(tuple(r) for r in extract(spec,path))
            if not projected:
                raise ValueError(f"No source rows from {path}")
            rows.update(projected)
    with args.output.open("w",newline="") as out:
        w=csv.writer(out,lineterminator="\n");w.writerow(FIELDS);w.writerows(sorted(rows))


if __name__=="__main__":
    main()
