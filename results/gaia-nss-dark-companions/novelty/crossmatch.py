"""Exact-ID prior-art audit. Missing identifications do not establish novelty."""
import argparse
import csv
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

from lxml import html

HERE = Path(__file__).resolve().parent
ARXIV = "https://arxiv.org/html/"
VIZ = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source="
ZEN = "https://zenodo.org/api/records/19181131/files/"
# Each tuple gives a table element ID and the zero-based Gaia ID column.
PAPERS = {
    "ands": ("Andrews et al. 2022", "2207.00680", [("S2.T1",0)]),
    "arengaia": ("Gaia Collaboration DR3 multiplicity", "2206.05595", [("S7.T9",0),("S7.T10",0)]),
    "ns": ("El-Badry et al. 2024 NS sample and prior candidates", "2405.00089v2", [("S2.T1",8),("A1.T3",0)]),
    "stewart": ("Stewart et al. 2026 follow-up", "2603.20371v1", [("S2.T1",0)]),
    "el": ("El-Badry et al. 2026 all selected targets", "2608.06453v2", [("A2.tab1",0),("A3.tab1",0)]),
    "zhang": ("Li et al. 2026 field sample, printed excerpt only", "2610.34715v1", [("S2.T1",0),("S3.T2",0)]),
    "liu": ("Liu et al. 2026 cluster samples, printed IDs", "2609.39925", [("S2.T1",0),("S3.T2",0),("S5.T3",0),("A1.T4",0)]),
    "lam26": ("LAMOST wide-orbit candidates 2026, DR3 intersection", "2602.05421", [("S3.T1",1)]),
}
URLS = {k:ARXIV+v[1] for k,v in PAPERS.items()}
URLS.update(andrew="https://zenodo.org/records/6779374/files/msco_candidates.csv?download=1",
            gaiaell="https://gea.esac.esa.int/tap-server/tap/sync?"+urllib.parse.urlencode(
                dict(request="doQuery",lang="ADQL",format="csv",query="SELECT source_id FROM gaiadr3.vari_compact_companion")),
            sfull=VIZ+"J/MNRAS/518/2991&-out.max=177",
            jfull=VIZ+"J/MNRAS/521/5927&-out.max=234")
for k in ("MS","RGB"):
    URLS[k] = ZEN+k+"+BH_candidates_from_GaiaDR3_summary_diagnostics.csv/content"


def tables(data, specs):
    doc = html.fromstring(data)
    for name, col in specs:
        tab, = doc.xpath('//*[@id=$name]', name=name)
        ids = set()
        for row in tab.xpath('.//tr'):
            cells = row.xpath('td|th')
            if len(cells)<=col:
                continue
            m = re.match(r"^\s*(\d{11,19})(?!\d)", cells[col].text_content())
            if m:
                ids.add(m.group(1))
        if not ids:
            raise ValueError(f"No Gaia IDs in {name}")
        yield name, ids


def asu(data, table, col, n):
    section = data.decode().split("#Table\t"+table+":",1)[1].split("#Table\t",1)[0]
    lines = [s for s in section.splitlines() if s.strip() and not s.startswith('#')]
    idx = lines[0].split('\t').index(col)
    ids = {s.split('\t')[idx].strip() for s in lines[3:]}
    if len(ids)!=n or not all(s.isdecimal() for s in ids):
        raise ValueError(f"VizieR table {table} failed completeness check")
    return ids


def main(cache):
    cache.mkdir(parents=True,exist_ok=True)
    with (HERE.parent/'gaia/mass_screen.csv').open() as f:
        target = {r['source_id'] for r in csv.DictReader(f)}
    matches, manifest = [], []
    for key,url in URLS.items():
        path=cache/key
        if not path.exists():
            with urllib.request.urlopen(url,timeout=90) as response:
                path.write_bytes(response.read())
        data=path.read_bytes()
        if key in PAPERS:
            title, _, specs = PAPERS[key]
            groups = list(tables(data,specs))
        elif key in ('sfull','jfull'):
            shahaf = key=='sfull'
            title = 'Shahaf et al. 2023' if shahaf else 'Jayasinghe et al. 2023 all screened SB1s'
            table = 'J_MNRAS_518_2991_table2' if shahaf else 'J_MNRAS_521_5927_table'
            groups=[(table,asu(data,table,'GaiaDR3' if shahaf else 'Source',177 if shahaf else 234))]
        else:
            title = {'andrew':'Andrew et al. 2022 all tiers','gaiaell':'Gaia DR3 ellipsoidal compact candidates',
                     'MS':'Müller-Horn et al. 2026 MS parent','RGB':'Müller-Horn et al. 2026 RGB parent'}[key]
            ids = {r['source_id'] for r in csv.DictReader(io.StringIO(data.decode()))}
            assert all(s.isdecimal() for s in ids)
            groups=[('parent' if key in ('MS','RGB') else 'catalog',ids)]
        for tab,ids in groups:
            manifest.append(dict(study=title,table=tab,unique_ids=len(ids),overlap=len(ids & target),
                input_key=key,url=url,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
            for sid in sorted(ids & target,key=int):
                matches.append(dict(source_id=sid,study=title,table=tab))
    with (HERE/'prior_matches.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['source_id','study','table'],lineterminator='\n')
        w.writeheader();w.writerows(sorted(matches,key=lambda r:(int(r['source_id']),r['study'],r['table'])))
    (HERE/'inputs.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path)
    main(p.parse_args().cache)
