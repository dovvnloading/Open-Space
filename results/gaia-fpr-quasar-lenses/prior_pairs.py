"""Supplemental vetoes for previously discussed quasar pairs and HSC lenses."""
import csv
import hashlib
import io
import json
import re
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from scipy.spatial import cKDTree


def vectors(ra, dec):
    a, d = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack((np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)))


def parse_sugohi(text):
    """Retain every valid coordinate in the authors' full-list download."""
    lines = text.lstrip("\ufeff").splitlines(keepends=True)
    first = next((line for line in lines if line.strip() and not line.lstrip().startswith("#")), "")
    # The current release starts directly with name, decimal RA and decimal Dec.
    # Accommodate an explicitly column-labelled release as well.
    delimiter = max([",", "\t", ";", "|"], key=first.count)
    reader = csv.reader(lines, delimiter=delimiter, strict=True)
    name_col, ra_col, dec_col = 0, 1, 2
    positions, exceptions = [], []
    previous, skipped, invalid, agreement = 0, 0, 0, 0
    nonstandard, discrepant = 0, 0
    parsed_to_end = True
    compact = re.compile(r"(?<!\d)(\d{2})(\d{2})(\d{2}(?:\.\d+)?)([+-])"
                         r"(\d{2})(\d{2})(\d{2}(?:\.\d+)?)(?!\d)")
    try:
        for row in reader:
            start, end = previous+1, reader.line_num
            raw = "".join(lines[previous:end]).rstrip("\r\n")
            previous = end
            if not raw.strip():
                skipped += 1
                continue
            names = [re.sub("[^a-z0-9]", "", token.lower()) for token in row]
            ra_head = next((i for i,n in enumerate(names) if n in {"ra","radeg","raj2000"}), None)
            de_head = next((i for i,n in enumerate(names) if n in {"de","dec","decdeg","dej2000"}), None)
            if not positions and not invalid and ra_head is not None and de_head is not None:
                name_col = next((i for i,n in enumerate(names) if n in {"name","id","object","objectname"}),0)
                ra_col, dec_col = ra_head, de_head
                skipped += 1
                continue
            if raw.lstrip().startswith("#"):
                skipped += 1
                continue
            try:
                ra, dec = float(row[ra_col]), float(row[dec_col])
                if not (0 <= ra < 360 and -90 <= dec <= 90):
                    raise ValueError("coordinates outside valid range")
            except (IndexError, ValueError) as exc:
                invalid += 1
                exceptions.append(dict(line=start, end_line=end, reason=str(exc), raw=raw))
                continue
            name = (row[name_col].strip() if name_col<len(row) else "") or f"SuGOHI row {start}"
            positions.append(dict(name=name, ra=ra, dec=dec))
            match = compact.search(name)
            if match:
                hh, mm, ss, sign, dd, dm, ds = match.groups()
                id_ra = 15*(int(hh)+int(mm)/60+float(ss)/3600)
                id_dec = (-1 if sign=="-" else 1)*(int(dd)+int(dm)/60+float(ds)/3600)
                delta_ra = abs((ra-id_ra+180)%360-180)*3600
                delta_dec = abs(dec-id_dec)*3600
                # The name truncates RA to seconds of time and Dec to arcsec.
                if (int(hh)<24 and int(mm)<60 and float(ss)<60 and int(dm)<60 and float(ds)<60
                        and abs(id_dec)<=90 and delta_ra<=16 and delta_dec<=1.5):
                    agreement += 1
                    continue
                reason = "name coordinates differ; included the decimal coordinate columns"
                discrepant += 1
            else:
                reason = "nonstandard designation; included the decimal coordinate columns"
                nonstandard += 1
            exceptions.append(dict(line=start, end_line=end, reason=reason, raw=raw))
    except csv.Error as exc:
        parsed_to_end = False
        exceptions.append(dict(line=previous+1, end_line=len(lines), reason=str(exc),
                               raw="".join(lines[previous:]).rstrip("\r\n")))
    assert len(positions)==agreement+nonstandard+discrepant
    audit = dict(rows=len(positions), data_records=len(positions)+invalid if parsed_to_end else None,
                 parsed_to_end=parsed_to_end, last_parsed_line=previous, physical_lines=len(lines),
                 ignored_headers_comments_blank=skipped, invalid_records=invalid,
                 numeric_with_compact_id_agreement=agreement, numeric_with_discrepant_id=discrepant,
                 numeric_with_nonstandard_id=nonstandard, exceptions=exceptions,
                 complete=(parsed_to_end and len(positions)>=2000 and invalid==0 and agreement>len(positions)/2),
                 completeness_scope="Complete requires valid decimal coordinates for every data record in "
                    "the checksummed author-provided full-list file and agreement with the compact sky name in a majority. "
                    "The 2000-row minimum is a format/truncation guard. No external total is assumed.")
    return pd.DataFrame(positions, columns=["name","ra","dec"]), audit


def extra_catalogues(out):
    out = Path(out)
    rows, receipts = [], []

    def get(url, params=None):
        res = requests.get(url, params=params, timeout=(8, 45))
        res.raise_for_status()
        receipts.append(dict(url=res.url, size=len(res.content),
                             sha256=hashlib.sha256(res.content).hexdigest()))
        return res

    def numeric(frame, ra, dec, name=None):
        x = pd.DataFrame({"ra": pd.to_numeric(frame[ra], errors="coerce"),
                          "dec": pd.to_numeric(frame[dec], errors="coerce")})
        x["name"] = frame[name].astype(str) if name else frame.index.astype(str)
        return x[x.ra.between(0, 360) & x.dec.between(-90, 90)].copy()

    def vizier(table):
        for host in ["vizier.cds.unistra.fr", "vizier.cfa.harvard.edu"]:
            try:
                r = get("https://" + host + "/viz-bin/asu-tsv", {
                    "-source": table, "-out.max": "unlimited", "-out": "**",
                    "-out.add": "_RAJ2000,_DEJ2000", "-oc.form": "d"})
                df = pd.read_csv(io.BytesIO(r.content), sep="\t", comment="#", dtype=str)
                df.columns = df.columns.str.strip()
                for ra, dec in [("_RAJ2000", "_DEJ2000"), ("RAJ2000", "DEJ2000"),
                                ("RAdegA", "DEdegA"), ("RA_ICRS", "DE_ICRS")]:
                    if ra in df and dec in df:
                        return numeric(df, ra, dec, next((x for x in ["Name","ID","recno"] if x in df),None))
                raise ValueError("No coordinate columns in " + str(list(df)[:30]))
            except (ValueError, KeyError, requests.RequestException) as exc:
                receipts.append(dict(table=table, host=host, error=str(exc)))
        raise ValueError("No readable coordinate table")

    for tag, table, expected in [
        ("Chen et al. 2026 MGQPC", "J/A+A/707/A30/mgqpc", 4112),
        ("SuGOHI VI", "J/A+A/642/A148/sugohi", 2000),
    ]:
        try:
            try:
                records = vizier(table)
                if len(records) != expected:
                    raise ValueError("Public release population changed")
            except Exception:
                # The first runner deposited a validated coordinate extract of
                # these same immutable CDS tables. Check its original receipt.
                snapshot = ("https://raw.githubusercontent.com/dovvnloading/Open-Space/"
                            "1c172928c6d266c8e54d4960edf60b286f8f3884/"
                            "results/gaia-fpr-quasar-lenses/evidence/")
                original = get(snapshot+"supplemental_receipts.json").json()
                if not any(r.get("catalog")==tag and r.get("rows")==expected and r.get("complete") for r in original):
                    raise ValueError("Prior input receipt missing")
                extract = get(snapshot+"supplemental_positions.csv.gz")
                all_rows = pd.read_csv(io.BytesIO(extract.content),compression="gzip",dtype={"name":str})
                records = all_rows.loc[all_rows.catalog==tag].copy()
                if len(records)!=expected:
                    raise ValueError("Archived subset differs from source count")
                receipts.append(dict(table=table,source_receipts=original,recovered_from=snapshot))
            receipts.append(dict(catalog=tag, rows=len(records), expected=expected, complete=len(records) == expected))
            records["catalog"] = tag
            rows.append(records)
        except Exception as exc:
            receipts.append(dict(catalog=tag, error=str(exc),complete=False))

    # The authors maintain a version extending the fixed 2020 CDS release.
    tag = "SuGOHI current compilation"
    try:
        page = get("https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/")
        index = BeautifulSoup(page.content, "html.parser")
        link = next(a["href"] for a in index.select("a[href]") if "download the full" in a.get_text(" ",strip=True).lower())
        data = get(urljoin(page.url, link))
        records, audit = parse_sugohi(data.text)
        receipts.append(dict(catalog=tag, **audit))
        # Even an incompletely parsed source can positively veto a known object.
        if not records.empty:
            records["catalog"] = tag
            rows.append(records)
    except Exception as exc:
        receipts.append(dict(catalog=tag,error=str(exc),complete=False))

    (out / "supplemental_receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    result = pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(columns=["ra","dec","name","catalog"])
    result.to_csv(out / "supplemental_positions.csv.gz", index=False)


def matches(components, out):
    refs = pd.read_csv(Path(out) / "supplemental_positions.csv.gz")
    tree = cKDTree(vectors(refs.ra, refs.dec))
    radius = 2 * np.sin(np.deg2rad(6 / 3600) / 2)
    result = []
    for row, nearby in zip(components.itertuples(),tree.query_ball_point(vectors(components.ra,components.dec),radius)):
        for i in nearby:
            ref = refs.iloc[i]
            result.append(dict(system_id=str(row.system_id), catalog=ref.catalog, name=ref["name"], ra=ref.ra, dec=ref.dec))
    result = pd.DataFrame(result,columns=["system_id","catalog","name","ra","dec"]).drop_duplicates()
    result.to_csv(Path(out) / "supplemental_matches.csv",index=False)
    return set(result.system_id)
