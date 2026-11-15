"""Compare both members and both periods with the 781-system Kolář compilation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from check_catalogs import vectors

SOURCE = "https://is.muni.cz/th/m79qa/phd_thesis_kolar_jakub.pdf"
ROW = re.compile(r"(?P<ra>\d+\.\d+)\s+(?P<dec>[+−–-]?\d+\.\d+)\s+"
                 r"(?P<loc>Gal|SMC|LMC)\s+[01]\s+(?P<pa>\d+\.?\d*)\s+(?P<pb>\d+\.?\d*)")


def read_appendix(text):
    # The earlier article reprint shows a short extract. Use the complete final table.
    start = text.rindex("TIC 418699570")
    entries = []
    for line in text[start:].splitlines():
        m = ROW.search(line)
        if m:
            dec = m["dec"].replace("−","-").replace("–","-")
            entries.append((line[:m.start()].strip(),float(m["ra"]),float(dec),float(m["pa"]),float(m["pb"])))
    if len(entries) != 781:
        raise ValueError(f"Expected 781 entries, found {len(entries)}")
    return entries


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("thesis_pdf",type=Path)
    p.add_argument("pairs",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    layout=subprocess.check_output(["pdftotext","-layout",str(a.thesis_pdf),"-"]).decode()
    entries=read_appendix(layout)
    sky=vectors([r[1] for r in entries],[r[2] for r in entries])
    clocks=np.sort(np.array([[r[3],r[4]] for r in entries]),axis=1)
    out={}
    for r in csv.DictReader(a.pairs.open()):
        stars=vectors([float(r["ra1"]),float(r["ra2"])],[float(r["dec1"]),float(r["dec2"])])
        chord=np.linalg.norm(sky[:,None,:]-stars[None,:,:],axis=2).min(axis=1)
        i = int(np.argmin(chord))
        periods=np.sort([float(r["period1"]),float(r["period2"])])
        error=np.max(np.abs(clocks/periods-1),axis=1)
        j=int(np.argmin(error))
        out[r["source_id1"]]={
            "nearest_position":entries[i][0], "nearest_position_arcsec":round(np.degrees(2*np.arcsin(chord[i]/2))*3600,2),
            "closest_period_pair":entries[j][0], "closest_period_max_fractional_error":round(error[j],7),
        }
    receipt={"source":SOURCE,"pdf_sha256":hashlib.sha256(a.thesis_pdf.read_bytes()).hexdigest(),
             "catalog_size":len(entries),"pairs":out}
    a.output.write_text(json.dumps(receipt,indent=2)+"\n")
