"""Extract the actual additive flux calibration from the hash-verified alert pages.

A positive reported reference flux permits a source-flux sensitivity bound in
that survey's own photometric units. A zero or missing reference remains unknown.
"""
import csv
import json
import re
from pathlib import Path
from lxml import html
from fetch import OUT,SHA,write_csv,get

CACHE=Path('.survey-cache')
PAGES=CACHE/'pages'
PATTERN=re.compile(r'([IHK])\s*=\s*([0-9.]+)\s*-\s*2\.5\s*log\s*10\s*\(\s*Δ\s*Flux\s*([+-])\s*(-?[\d.]+)\s*\)',re.I)

def main():
    end={r['name']:r for r in csv.DictReader((OUT/'index_receipts.csv').open())}
    raw_p=(CACHE/'prime2024idx').read_bytes();raw_m=(CACHE/'moa2024').read_bytes()
    assert SHA(raw_p)==end['prime2024idx']['input_sha256'] and SHA(raw_m)==end['moa2024']['input_sha256']
    notes={k:v.get('remark') for k,v in json.loads(raw_p).items()}
    notes.update({r[0]:r[-1] for r in json.loads(raw_m)['entries']})
    out=[]
    for r in csv.DictReader((OUT/'photometry.csv').open()):
        name=r['event'];page=PAGES/(name+'.html')
        data=page.read_bytes() if page.exists() else get(r['event_url'])
        if SHA(data)!=r['page_sha256']:raise ValueError('source page differs from frozen version '+name)
        doc=html.fromstring(data)
        blocks=doc.xpath('//*[@id="calib"]')
        match=PATTERN.search(blocks[0].text_content()) if blocks else None
        if blocks and not match:raise ValueError((name,blocks[0].text_content()))
        entry=dict(event=name,page_sha256=r['page_sha256'],archive_index_remark=notes.get(name) or '',band='',
                   zero_point='',additive_reference_counts='',reference_status='no transform on page')
        if match:
            band,z,sgn,c=match.groups()
            ref=float(c)*(1 if sgn=='+' else -1)
            entry.update(band=band,zero_point=float(z),additive_reference_counts=ref,
                         reference_status='positive reference' if ref>0 else 'zero or negative reference')
        out.append(entry)
    write_csv(OUT/'page_calibrations.csv',out)


if __name__=='__main__':main()
