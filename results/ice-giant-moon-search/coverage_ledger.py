"""Join each MAST observation group to actual detector work and its published evidence."""
from collections import defaultdict
import csv
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent


def read(path):
    with (ROOT/path).open() as f:return list(csv.DictReader(f))


def ledger():
    obs=read('manifests/observations.csv')
    parent=defaultdict(set)
    for r in read('manifests/products.csv'):
        if r['productSubGroupDescription'] in ('FLT','CAL'):
            for p in r['parent_obsids'].split('|'):parent[p].add(r['productFilename'])
    evidence=defaultdict(lambda:defaultdict(set))
    def add(file,level,href):evidence[file][level].add(href)
    for r in read('hst_multiband/exposure_hashes.csv'):
        file=r['uri'].split('/')[-1]
        add(file,'header_or_registration','hst_multiband/exposure_hashes.csv')
        if r['status']=='used' and r['selection_status']=='screened':add(file,'orbit_screen','hst_multiband/epoch_index.csv')
    for path in sorted(ROOT.glob('**/*inputs*.json')):
        if 'epochs' in path.parts:continue
        record=json.loads(path.read_text())
        if isinstance(record,dict):record=record.get('tile_inputs',[])
        for row in record:
            if not isinstance(row,dict):continue
            n=row.get('filename','')
            if n.endswith('_cal.fits') and 'exposure_id' in row:
                add(n,'science_pixels',str(path.relative_to(ROOT)))
    for r in read('miri/footprints.csv'):add(r['file'],'footprint','miri/footprints.csv')
    for r in read('outer6379/inputs.csv'):add(r['filename'],'orbit_screen','outer6379/README.md')
    for fold in ('jwst2739','jwst2739_uranus','jwst2768','jwst8975','jwst6379_9482','miri'):
        for p in (ROOT/fold).glob('*inputs*.json'):
            rec=json.loads(p.read_text()); rec=rec.get('tile_inputs',[]) if isinstance(rec,dict) else rec
            for r in rec:
                if isinstance(r,dict) and r.get('filename','').endswith('_cal.fits'):
                    # A chip used in a mosaic may contribute off-center pixels only.
                    add(r['filename'],'orbit_screen',str(p.relative_to(ROOT)))
    for r in json.loads((ROOT/'jwst9482/inputs.json').read_text()):add(r['filename'],'orbit_screen','jwst6379_9482/README.md')
    for p in (ROOT/'nircam_supplement').glob('*cutouts.csv'):
        for r in read(p):
            add(r['source'],'footprint',str(p.relative_to(ROOT)))
            if r['filename']:add(r['source'],'science_pixels',str(p.relative_to(ROOT)))
    for r in read('nircam_supplement/registrations.csv'):
        add(r['filename'],'orbit_screen','nircam_supplement/visits.csv')
    # Dedicated 2023 and extra MIRI reductions may also record detections and frame membership.
    for directory in ['jwst2739_uranus','miri/extra']:
        for p in (ROOT/directory).rglob('*.csv'):
            if 'ranges' in p.name:continue
            for row in read(p):
                for k in ('filename','source','file','exposure_id'):
                    name=row.get(k,'')
                    if name.endswith('_cal'):name+='.fits'
                    if not name.endswith('_cal.fits'):continue
                    if 'cutout' in p.name or 'input' in p.name:
                        add(name,'footprint',str(p.relative_to(ROOT)))
                if 'registration' in p.name or 'centers' in p.name:
                    name=row.get('filename',row.get('exposure_id',''))
                    if name.endswith('_cal'):name+='.fits'
                    if name.endswith('_cal.fits'):add(name,'science_pixels',str(p.relative_to(ROOT)))
    for p in (ROOT/'miri/extra').glob('*_images.csv'):
        for row in read(p):
            if int(row['usable_pixels'])>0:add(row['file'],'orbit_screen',str(p.relative_to(ROOT)))
    result=[]
    for r in obs:
        files=parent[r['obsid']]
        levels=defaultdict(set);links=set()
        for name in files:
            for kind,refs in evidence[name].items():
                levels[kind].add(name);links.update(refs)
        work=next((stage for stage in ('orbit_screen','science_pixels','footprint','header_or_registration') if levels[stage]),'inventory_only')
        result.append(dict(planet=r['planet'],obsid=r['obsid'],observation=r['obs_id'],instrument=r['instrument_name'],filter=r['filters'],
            n_detector_exposures=len(files),highest_recorded_stage=work,detectors_recorded_in_orbit_screen=len(levels['orbit_screen']),
            evidence='|'.join(sorted(links))))
    with (ROOT/'manifests/processing_ledger.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)


if __name__=='__main__':ledger()
