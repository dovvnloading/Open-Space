"""Export a frozen manifest of pySIS bundle URLs, stream members and content hashes."""
import argparse
import csv
import json
from pathlib import Path
from freeze import HERE


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args()
    rows=list(map(json.loads,(args.cache/'index.jsonl').read_text().splitlines()))
    fields=('event','status','utc','size','sha256','n_I','n_V','members','url','error')
    for year in (2024,2025):
        with (HERE/f'pySIS-{year}.tsv').open('w') as h:
            writer=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            for r in rows:
                event=r['event']
                if str(year) not in event:continue
                receipt=args.cache/(event+'.receipt.json');err=args.cache/(event+'.error.json')
                out=dict(event=event,url=r['url'],status='pending')
                if receipt.exists():
                    s=json.loads(receipt.read_text());curves=s['curves'];out.update({
                        k:s[k] for k in ('status','utc','size','sha256')})
                    out.update(n_I=sum('_I.' in n for n in curves),n_V=sum('_V.' in n for n in curves),members='|'.join(curves))
                elif err.exists():
                    s=json.loads(err.read_text());out.update(status=s['status'] or 'failed',error=s['error'],utc=s['utc'])
                writer.writerow(out)


if __name__=='__main__':
    main()
