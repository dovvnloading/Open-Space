"""Metadata-only guard; no image-reading entry point is provided."""
import datetime as dt
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXPECTED = ('202611010000', '202611010008', '202611010016')

def verify(provenance, parents, phase='metadata'):
    if phase != 'metadata':
        raise ValueError('science phases are not authorized or implemented')
    if len(provenance) != 3 or len(parents) != 9:
        raise ValueError('unexpected population')
    byname = {p['name']: p for p in parents}
    if len(byname) != 9:
        raise ValueError('duplicate parent')
    allnames = []
    windows = []
    for item, timestamp in zip(provenance, EXPECTED):
        if 'PUNCH_L2_CTM_'+timestamp not in item['url'] or len(item['parents']) != 3:
            raise ValueError('unexpected product')
        rows=[]
        for i,name in enumerate(item['parents'],1):
            if name not in byname or f'_CR{i}_'+timestamp not in name or not name.endswith('_v0l.fits'):
                raise ValueError('unexpected parent')
            h=byname[name]['headers']
            if h['LEVEL'] != '1' or h['OBSCODE'] != str(i) or h['OBSTYPE'] != 'Unpolarized WFI image':
                raise ValueError('invalid instrument')
            begin=dt.datetime.fromisoformat(h['DATE-BEG']); end=dt.datetime.fromisoformat(h['DATE-END'])
            if not begin < end: raise ValueError('invalid interval')
            rows.append((begin,end)); allnames.append(name)
        windows.append((min(x[0] for x in rows),max(x[1] for x in rows)))
    if len(set(allnames)) != 9 or any(windows[i][1] >= windows[i+1][0] for i in range(2)):
        raise ValueError('overlapping exposures or parents')
    return {'phase':'metadata', 'epochs':3, 'parents':9, 'windows':[[a.isoformat(),b.isoformat()] for a,b in windows]}

if __name__ == '__main__':
    print(json.dumps(verify(json.loads((ROOT/'provenance.json').read_text()),json.loads((ROOT/'parent_headers.json').read_text())),indent=2))
