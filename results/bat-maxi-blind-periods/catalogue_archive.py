"""Decode VizieR binary or ASCII table serializations and archived CDS tables."""

import gzip
import io
import re
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import ascii, votable


CDS_HOSTS = ("https://cdsarc.cds.unistra.fr/ftp", "https://cdsarc.u-strasbg.fr/ftp",
             "https://cdsarc.cds.unistra.fr/ftp/cats", "https://vizier.cfa.harvard.edu/ftp/cats")


def display(value):
    if value is None or np.ma.is_masked(value):
        return ""
    if isinstance(value,bytes):
        return value.decode("utf-8",errors="replace").strip()
    return str(value).strip()


def binary_tables(raw):
    result = []
    for table in votable.parse(io.BytesIO(raw),verify="ignore").iter_tables():
        if table.array is None or not len(table.array):
            continue
        fields = [dict(name=f.name or f.ID,ucd=str(f.ucd or ""),unit=str(f.unit or ""),
                       desc=str(f.description or "")) for f in table.fields]
        rows = [{f["name"]:display(v) for f,v in zip(fields,row)} for row in table.array]
        result.append(dict(name=table.name or table.ID or "",fields=fields,rows=rows))
    return result


def readme(archive,vid):
    errors = []
    for base in CDS_HOSTS:
        try:
            blob,m=archive.fetch(f"{base}/{vid}/ReadMe")
            txt=blob.decode("utf-8",errors="replace")
            if not re.search("Byte-by-byte Description of file",txt,re.I):
                raise ValueError("unexpected CDS ReadMe content")
            return txt,m,base
        except ValueError as exc:
            errors.append(str(exc))
    raise ValueError(" | ".join(errors))


def readme_title(txt):
    return " ".join(txt.split("Description:")[0].split("File Summary:")[0].split()[:250])


def cds_tables(archive,vid,check_title):
    txt,metadata,base=readme(archive,vid)
    title=readme_title(txt)
    if not check_title(title):
        raise ValueError("unexpected official catalogue title "+title[:200])
    # The File Summary lists the actual data files; following byte-description
    # headings can wildcard several files at once.
    file_summary = re.split("File Summary:", txt, flags=re.I)
    if len(file_summary)<2:
        raise ValueError("no CDS File Summary")
    section=re.split("Byte-by-byte|See also:|Description:",file_summary[1],flags=re.I)[0]
    names=list(dict.fromkeys(re.findall(r"^\s*([A-Za-z0-9_.+-]+\.dat(?:\.gz)?)\s+", section, re.M|re.I)))
    if not names:
        raise ValueError("no CDS table files in File Summary")
    results,sources=[],[]
    errors=[]
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        rp=root/"ReadMe"
        rp.write_text(txt)
        for filename in names:
            try:
                table_meta=None
                url=f"{base}/{vid}/{filename}"
                for alt in ([url] if filename.endswith(".gz") else [url,url+".gz"]):
                    try:
                        data,table_meta=archive.fetch(alt)
                        if data.startswith(b"\x1f\x8b"):
                            data=gzip.decompress(data)
                        if not data.strip() or data.lstrip().startswith(b"<"):
                            raise ValueError("unexpected data file type")
                        break
                    except (ValueError,OSError):
                        table_meta=None
                if table_meta is None:
                    raise ValueError("missing table file")
                fname=filename.removesuffix(".gz")
                path=root/fname
                path.write_bytes(data)
                t=ascii.read(str(path),format="cds",readme=str(rp),guess=False)
                if not len(t):
                    raise ValueError("empty data table")
                fields=[dict(name=k,unit=str(t[k].unit or ""),ucd="",desc=str(t[k].description or "")) for k in t.colnames]
                rows=[{k:display(r[k]) for k in t.colnames} for r in t]
                results.append(dict(name=vid+"/"+fname.removesuffix(".dat"),fields=fields,rows=rows))
                sources.append(dict(file=filename,query=table_meta))
            except (ValueError,TypeError,OSError) as exc:
                errors.append(f"{filename}: {exc}")
    if not results or not any(any(re.search("orbital period|orb[.]?per|p_?orb",f["name"]+" "+f["desc"],re.I)
                                  for f in t["fields"]) for t in results):
        raise ValueError("no orbital-period-bearing CDS source table; "+" | ".join(errors))
    return results, metadata, title, sources, errors
