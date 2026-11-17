"""Search the public DIVA+ reduction/characterization tables by catalogued DIVA+ identity."""
import argparse,csv,hashlib,io,json
from urllib.parse import urlencode
from urllib.request import urlopen
from datetime import datetime,timezone
from pathlib import Path
def main(out):
    urls={
        'DIVA_HD38206_identity':('stars','1;2;4;6;42','1::eq::1150'),
        'DIVA_HD38206_reductions':('reductions','1;2;4;5;11;12;13;14;16;20;21;23;24;25;27;46','46::eq::1150'),
        'DIVA_HD38206_images':('images','1;2;3;4;5;9;13;14;17;19;21;23','2::in::7211|7713|8026'),
        'DIVA_HD38206_characters':('characterisations','1;2;3;6;7;8;9;10;11;12;13;37;38','2::eq::HIP 26966')}
    rec=[]
    for name,(tab,a,c) in urls.items():
        url='https://data.lam.fr/server/search/divap/'+tab+'?'+urlencode({'a':a,'c':c,'f':'csv'})
        with urlopen(url,timeout=60) as f:b=f.read()
        rows=list(csv.DictReader(io.StringIO(b.decode())))
        if name.endswith('identity'):assert len(rows)==1 and rows[0]['star_name']=='HIP 26966'
        if name.endswith('reductions'):assert all(r['star_name']=='HIP 26966' for r in rows) and len(rows)==3
        (out/(name+'.csv')).write_bytes(b)
        rec.append(dict(name=name,url=url,rows=len(rows),sha256=hashlib.sha256(b).hexdigest(),utc=datetime.now(timezone.utc).isoformat()))
    (out/'DIVA_HD38206_receipts.json').write_text(json.dumps(rec,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);main(p.parse_args().out)
