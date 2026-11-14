"""Refresh OGLE/EROS inputs in a new directory. The evidence package works offline."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
OGLE2 = "https://ogledb.astrouw.edu.pl/~ogle/photdb/"
VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
REPO_REV = "d3cec82599d4bf332debb9c7f1e1522e2893acc5"

def ogle2_url():
    fields = {"db_target":"smc", "dbtyp":"dia2", "use_field":"on", "val_field":"SMC_SC9",
              "use_starid":"on", "val_starid":"141913", "disp_ra":"on", "disp_decl":"on",
              "disp_imean":"on", "sort":"field", "sorting":"ASC", "pagelen":"50"}
    with urlopen(OGLE2+"query.php?qtype=phot&first=1",urlencode(fields).encode(),timeout=90) as resp:
        listing = resp.read().decode()
    detail = re.search(r"href='(getobj.php\?[^']+)'",listing).group(1)
    with urlopen(OGLE2+detail,timeout=90) as resp:
        entry = resp.read().decode()
    return OGLE2 + re.search(r"href=['\"](data/[^'\"]+\.dat)['\"]",entry).group(1)

def urls():
    base = "https://www.astrouw.edu.pl/ogle/"
    out = {
        "iii-04602-I.dat":base+"ogle3/OIII-CVS/smc/ecl/phot/I/OGLE-SMC-ECL-4602.dat",
        "iv-04602-I.dat":base+"ogle4/OCVS/smc/ecl/phot/I/OGLE-SMC-ECL-04602.dat",
        "epoch_measurements.csv":"https://raw.githubusercontent.com/dovvnloading/Open-Space/"
                +REPO_REV+"/results/smc-06090/population/epoch_measurements.csv",
    }
    for tab in ("xcat","xtime"):
        fields = {"-source": "II/390/"+tab, "-c": "15.7324167 -72.8732778",
                  "-c.rs": "8", "-out.max": "30000"}
        if tab == "xcat":
            fields["-out.all"] = "1"
        else:
            fields["-out"] = "EROS,HJD,Rmag,e_Rmag,Bmag,e_Bmag"
        out[f"04602_{tab}.tsv"] = VIZIER + urlencode(fields)
    return out

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output",type=Path,help="New directory for refreshed files")
    out = parser.parse_args().output
    out.mkdir(parents=True,exist_ok=False)
    sources = urls()
    sources["ii-04602-I.dat"] = ogle2_url()
    receipts = {}
    for name,url in sources.items():
        data=urlopen(url,timeout=120).read()
        (out/name).write_bytes(data)
        receipts[name] = {"url":url,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}
    (out/"receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")

if __name__ == "__main__":
    main()
