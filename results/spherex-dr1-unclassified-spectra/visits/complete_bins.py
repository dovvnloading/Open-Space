"""Try remaining wide-field footprints when an entire native wavelength bin failed quality cuts."""
import argparse
import csv
import json
from pathlib import Path
from astropy.io import fits
from measure import choose_images, fetch_one, analyse, cut_url

def run(a):
    with a.targets.open() as f:
        targets=list(csv.DictReader(f))
    with a.sources.open() as f:
        sources={r['source_id']:r for r in csv.DictReader(f)}
    for target in targets:
        name=target['source_id'];suffix=target.get('label','')
        stem=name+('_'+suffix if suffix else '')
        path=a.out/f'{stem}.csv'
        if not path.exists():
            continue
        with path.open() as f:
            rows=list(csv.DictReader(f))
        for r in rows:
            for k in ('mjd','wave_um','flux','sigma'):
                if r.get(k):
                    r[k]=float(r[k])
        existing={r['url'] for r in rows}
        info=json.loads((a.out/f'{stem}.json').read_text())
        ra,dec=info['ra'],info['dec'];center=info['center_um'];scale=info['scale_um'];d=int(target['detector'])
        cal=json.loads((a.calibration/f'calibration_D{d}.json').read_text())
        options=choose_images(a.index_cache/f'footprints-D{d}.csv',cal,ra,dec,center,scale,6,all_options=True)
        ok={(r['release'],r['bin']) for r in rows if r['status']=='ok'}
        extra=[]
        with fits.open(a.science_cache/f'sci_D{d}.fits') as h:
            maps=[h[k].data for k in ('CWAVE','CBAND')]
            for record in options:
                group=record['release'],record['bin']
                url=cut_url(record['uri'],ra,dec)
                if group in ok or url in existing:continue
                row=fetch_one(record,ra,dec,center,scale,a.index_cache/'cutouts',maps)
                row['selection_reason']='missing_bin_alternate'
                extra.append(row);existing.add(url)
                if row['status']=='ok':
                    ok.add(group)
        if extra:
            for row in rows:
                row['selection_reason']=row.get('selection_reason') or 'initial'
            rows+=extra
            keys=list(dict.fromkeys(k for row in rows for k in row))
            with path.open('w') as f:
                w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
            info.update(attempted=len(rows),status_counts={s:sum(r['status']==s for r in rows) for s in sorted({r['status'] for r in rows})},groups=analyse(rows))
            (a.out/f'{stem}.json').write_text(json.dumps(info,indent=2,allow_nan=False)+'\n')
        print(stem,'additional images',len(extra),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('targets','sources','index_cache','calibration','out','science_cache'):p.add_argument(name,type=Path)
    run(p.parse_args())
