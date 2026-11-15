"""Compare this dated catalogue snapshot with the public archive services."""
import csv
import io
import json
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

IDS = ["4318465062103815552", "4318465066397284608", "4318465066396437504"]
RA, DEC = 294.82786250815144, 14.930979608612361
COLS = "source_id ra dec pmra pmdec parallax phot_g_mean_mag ruwe duplicated_source"
EXTRA = ("ra_error dec_error pmra_error pmdec_error parallax_error ra_dec_corr "
         "ra_pmra_corr ra_pmdec_corr dec_pmra_corr dec_pmdec_corr pmra_pmdec_corr "
         "ipd_frac_multi_peak ipd_frac_odd_win visibility_periods_used")


def get(url, params):
    with urlopen(url + "?" + urlencode(params), timeout=120) as handle:
        return handle.read()


def query(adql):
    raw = get("https://gaia.aip.de/tap/sync",
              dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=adql))
    root = ET.fromstring(raw)
    status = root.findall(".//{*}INFO[@name='QUERY_STATUS']")
    if not status or any(info.get("value") != "OK" for info in status):
        raise ValueError("Gaia query incomplete: " + str([i.text for i in status]))
    names = [f.get("name") for f in root.findall(".//{*}FIELD")]
    rows = [dict(zip(names, [c.text or "" for c in row]))
            for row in root.findall(".//{*}TR")]
    return sorted(rows, key=lambda r: r["source_id"])


def snapshot():
    footprint = (f"1=CONTAINS(POINT('ICRS',ra,dec),"
                 f"CIRCLE('ICRS',{RA},{DEC},{80/3600}))")
    field = query(f"SELECT {','.join(COLS.split())} FROM gaiadr3.gaia_source WHERE {footprint}")
    extra = query(f"SELECT source_id,{','.join(EXTRA.split())} FROM gaiadr3.gaia_source "
                  f"WHERE source_id IN ({','.join(IDS)})")
    dr2ids = ["4318465062103815552", "4318465066397284608", "4318465066396437632"]
    dr2 = query("SELECT source_id,ra,dec,pmra,pmdec FROM gaiadr2.gaia_source "
                f"WHERE source_id IN ({','.join(dr2ids)})")
    pscols = "[objID,raMean,decMean,epochMean,nDetections]"
    ps = get("https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv",
             dict(ra=294.82785023, dec=14.93091909, radius=.012,
                  pagesize=10000, columns=pscols))
    ps = list(csv.DictReader(io.StringIO(ps.decode())))
    psids = ["125922948284312089", "125912948273943217", "125912948273282660"]
    ps = [next(r for r in ps if r["objID"] == i) for i in psids]
    tsv = get("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",
              {"-source":"I/352/gedr3dis", "-c":"294.82785 14.930919", "-c.rs":30,
               "-out.max":1000, "-out.all":""}).decode()
    lines = [line for line in tsv.splitlines() if line and not line.startswith("#")]
    distances = {r["Source"].strip(): {k:r[k].strip() for k in ["rgeo","b_rgeo","B_rgeo","rpgeo"]}
                 for r in csv.DictReader(lines, delimiter="\t") if r["Source"].strip() in IDS}
    assert len(extra) == len(dr2) == len(distances) == len(ps) == 3
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=COLS.split(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(field)
    return {"field_csv":out.getvalue(), "extra":extra, "ps1":ps, "distances":distances, "dr2":dr2}


if __name__ == "__main__":
    import sys
    from pathlib import Path
    data = snapshot()
    if "--emit" in sys.argv:
        print(json.dumps(data, separators=(",", ":")))
    else:
        here = Path(__file__).resolve().parent
        checked = json.loads((here / "followup.json").read_text())
        assert data.pop("field_csv").encode() == (here / "field.csv").read_bytes()
        assert data == checked
        print("Archive rows agree with the deposited snapshot.")
