"""Verify the released measurements and all file checksums offline, with Python 3."""
import csv
import hashlib
import json
import math
import zipfile
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def rows(path):
    with (ROOT/path).open() as f:
        return list(csv.DictReader(f))


def ensure(test, failure):
    if not test:
        raise ValueError(failure)


def main():
    inventory=(ROOT/'SHA256SUMS').read_text().splitlines()
    for line in inventory:
        expected,path=line.split('  ',1)
        ensure(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path)
    sources=rows('catalog/slde_screen.csv')
    pos={r['id_str'] for r in sources if r['model_state']=='positive_radius'}
    ensure(len(sources)==500 and len(pos)==312,'catalogue reconciliation')
    core={r['id_str']:float(r['source_subtracted_snr']) for r in rows('imaging/apertures.csv')
          if r['band']=='VIS' and r['aperture']=='core' and r['id_str'] in pos}
    ensure(set(core)==pos and 19.7<min(core.values())<19.8,'modelled VIS core detections')
    unmatched=set(r['id_str'] for r in sources)-pos
    targeted=json.loads((ROOT/'arc_search/unmodelled.json').read_text())
    fallback=json.loads((ROOT/'arc_search/fallback.json').read_text())
    tested={r['id_str'] for r in targeted+fallback if r['status']=='examined'}
    ensure(tested==unmatched and len(tested)==188,'parents without valid mass models')
    ensure(sum(r.get('pass_vis',0) for r in targeted+fallback)==1,'targeted candidates')
    strips=rows('arc_search/strips.csv')
    raw=rows('arc_search/vis_triggers.csv')
    ensure(len(strips)==32 and sum(int(r['eligible_binned_centers']) for r in strips)==27013228,'area denominator')
    ensure(len(raw)==1258 and sum(int(r['pass_vis']) for r in raw)==6,'independent shortlist')
    with zipfile.ZipFile(ROOT/'arc_search/frozen_stages.zip') as z:
        ensure(z.testzip() is None and len(z.namelist())==64,'frozen per-strip receipts')
        for r in strips:
            name=r['tile']+'_'+r['row']+'.json'
            one=json.loads(z.read('runs/'+name))
            ensure(one['sha256']==r['sha256'] and len(one['triggers'])==int(r['n_triggers']),name)
            ensure(sum(s['pass_vis'] for s in json.loads(z.read('profiles/'+name))['profiles'])==int(r['n_profiles_passing']),name)
    ray=rows('arc_search/lens_injections.csv')
    ensure(len(ray)==1491 and sum(int(r['preexisting']) for r in ray)==0,'ray-traced injections')
    good=Counter((r['source_offset_over_theta'],r['intrinsic_peak_sigma']) for r in ray if int(r['recovered']))
    ensure(good=={('0.02','3'):127,('0.02','8'):213,('0.18','3'):28,('0.18','8'):150},'recovery')
    models=rows('arc_search/model/cross_validation.csv')
    ordinary=[r for r in models if r['center_bound_arcsec']=='0.35']
    ensure(len(ordinary)==7 and sum(float(r['mean_square_knots'])<float(r['mean_square_lens']) for r in ordinary)==5,'held-out image fits')
    ensure(len(rows('arc_search/model/cross_validation_folds.csv'))==45,'fold count')
    native=rows('arc_search/validation/native_inputs.csv')
    ensure(len(native)==84 and len({r['sha256'] for r in native})==84,'native cutout provenance')
    bank=rows('arc_search/validation/conditional_limits.csv')
    ensure(len(bank)==700 and all(float(r['conditional_ML_floor'])>0 for r in bank),'conditional limits')
    b='102158888_2659753895653871052'
    p={r['z_l_model']:float(r['ML_p50']) for r in rows('inference/conditional_ML_posterior.csv')
       if r['id_str']==b and r['band']=='h' and r['zs']=='2'}
    ensure(math.isclose(p['Euclid_PHZ'],105.2,abs_tol=.1) and
           math.isclose(p['Legacy_DR9'],54.4,abs_tol=.1),'low-z conditional H posteriors')
    print(f'{len(inventory)} file hashes verified; catalogue, cutout, model, recovery and geometry checks passed.')

if __name__=='__main__': main()
