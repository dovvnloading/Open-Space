"""Analyze one disjoint inventory shard, checkpointing complete visit records."""
import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

# Small least-squares fits are repeated per product and per injection.
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import archive_io
import field
import photometry
import vetting


def _error(exc):
    return type(exc).__name__ + ': ' + str(exc)[:500]


def _family(path):
    name = path.name.lower()
    if ('pipe' in name or 'psf' in name or
            name.removesuffix('.gz').endswith(('_sa.fits','_im.fits'))):
        return 'PIPE'
    if 'sci_cor_lc' in name or 'drp' in name or 'lightcurve' in name:
        return 'DRP'
    return 'unclassified'


def _rank(path):
    name = path.name.lower()
    return (not ('default' in name), not ('optimal' in name), name)


def _close(x, y):
    return (x['source_id'] == y['source_id']
            and abs(x['epoch_bjd'] - y['epoch_bjd']) < .4 * max(x['duration_days'],y['duration_days'])
            and .4 <= x['duration_days'] / y['duration_days'] <= 2.5)


def _coordinate(meta, visit, header, key, target):
    v = meta.get(header)
    if v is None and target:
        v = visit.get(key)
    return v


def _write_refs(literature, folder):
    # Titles and abstracts recur across events. Store the originals once by digest.
    refs = {k: literature.pop(k) for k in ('associated_references', 'arxiv_records') if k in literature}
    if not refs:
        return
    # Identifiers and hashes allow independent reconciliation without republishing abstracts.
    for paper in refs.get('arxiv_records', []):
        abstract = paper.pop('abstract', None)
        if abstract:
            paper['abstract_sha256'] = hashlib.sha256(abstract.encode()).hexdigest()
    data = json.dumps(refs, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode()
    key = hashlib.sha256(data).hexdigest()
    path = folder / (key + '.json.gz')
    folder.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with path.open('wb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', mtime=0) as f:
            f.write(data)
    literature.update(reference_data=folder.name + '/' + path.name, uncompressed_sha256=key,
                      reference_counts={k: len(v) for k, v in refs.items()})


def _classify(group, reductions, paths, reference, catalogs, cache, out, alternative):
    group = sorted(group, key=lambda e: (e['photometry_pass'], e['snr']), reverse=True)
    event = dict(group[0])
    matches = [g for g in group if g['photometry_pass']]
    families = sorted(set(g['reduction'] for g in matches))
    independent = {'DRP', 'PIPE'}.issubset(families)
    if independent:
        drp = next(g for g in matches if g['reduction'] == 'DRP')
        pipe = next(g for g in matches if g['reduction'] == 'PIPE')
        if .5 <= drp['depth']/pipe['depth'] <= 2.:
            cross = 'DRP and PIPE agree'
        else:
            cross = 'DRP and PIPE depth discrepancy'
    elif {'DRP', 'PIPE'}.issubset(reductions):
        cross = 'inconclusive across available DRP and PIPE products'
    else:
        cross = 'single reduction family'
    event.update(matching_products=[dict(product=g['product'], reduction=g['reduction'],
                                         epoch_bjd=g['epoch_bjd'],duration_days=g['duration_days'],
                                         depth=g['depth'],snr=g['snr'],checks=g['checks']) for g in group],
                 cross_reduction=cross)
    try:
        event['localization'] = field.locate(event, paths, reference,
                                             (event['ra_deg'], event['dec_deg']))
    except Exception as exc:
        event['localization'] = dict(status='incomplete', error=_error(exc))
    try:
        match = catalogs.match(event['star'], event['ra_deg'], event['dec_deg'])
        event['alternate_aperture_depths'] = {}
        for aperture in alternative if matches and event['source_kind']=='target' else []:
            try:
                check = photometry.forced_depth(aperture, event, match['ephemerides'], reference)
            except Exception as exc:
                check = dict(status='incomplete', error=_error(exc))
            event['alternate_aperture_depths'][aperture.name] = check
        v = vetting.vet_event(event, match, cache)
        _write_refs(v.get('checks', {}).get('literature', {}), out)
    except Exception as exc:
        v = dict(disposition='requires review', errors=[_error(exc)], checks={})
    event['vetting'] = v
    checks = v.get('checks', {})
    gaia = checks.get('field', {})
    tess = checks.get('tess', {})
    literature = checks.get('literature', {})
    # Unknown coverage or missing phase information remains visible as unresolved.
    if v.get('known_ephemeris_overlaps'):
        disposition = 'catalog ephemeris overlap'
    elif not matches and all(any(answer is False for answer in g['checks'].values()) for g in group):
        disposition = 'photometric rejection'
    elif (not matches or event['localization'].get('status') != 'consistent'
          or 'inconclusive' in cross or 'discrepancy' in cross
          or any(c['status']!='agreement' for c in event.get('alternate_aperture_depths',{}).values())
          or gaia.get('status') != 'screened'
          or tess.get('status') not in ('compared','no TESS sector footprint')
          or 'reference_data' not in literature or literature.get('truncated')
          or v.get('unresolved') or v.get('catalog_errors') or v.get('errors')):
        disposition = 'unresolved'
    else:
        # Reading actual publications, adjudicating Gaia neighbors and TESS comparison
        # remain human review steps before publishing a claim of novelty.
        disposition = 'requires review'
    event['disposition'] = disposition
    return event


def _run_visit(visit, folder, catalogs, cache, out):
    rec = dict(visit_id=visit['visit_id'], star=visit.get('star'), ra_deg=visit.get('ra'),
               dec_deg=visit.get('dec'), receipts=[], field={}, products=[], events=[],
               masks={}, errors=[], status='unsearched')
    raw_events = []
    families = {}
    paths = []
    try:
        paths = archive_io.fetch_visit(visit, folder, include_field=True)
        source_meta = {}
        for path in paths:
            if 'subarray' in path.name.lower() or 'starcatalog' in path.name.lower():
                continue
            try:
                source_meta[path] = photometry.metadata(path)
            except Exception as exc:
                rec['errors'].append(dict(product=path.name, error=_error(exc)))
        reference = next(iter(sorted((p for p in source_meta if _family(p) == 'DRP'), key=_rank)), None)

        def process(path, scope, meta, extracted=False):
            family = 'field aperture' if extracted else _family(path)
            star = meta.get('star') or (visit.get('star') if scope=='target' else None)
            ra = _coordinate(meta, visit, 'ra_deg', 'ra', scope=='target')
            dec = _coordinate(meta, visit, 'dec_deg', 'dec', scope=='target')
            sid = 'target' if scope == 'target' else str(star) + ' ' + str(ra) + ' ' + str(dec)
            # Distinct photometric extractions inherit the same source mask when matched.
            match = catalogs.match(star or '', ra, dec)
            rec['masks'][sid] = match
            result = photometry.analyze(path, match['ephemerides'], reference_path=reference)
            result.update(name=path.name, reduction=family, source_kind=scope, source_id=sid,
                          star=star, ra_deg=ra, dec_deg=dec)
            if extracted:
                result['derived_sha256'] = archive_io.digest(path)
            for ev in result.pop('events', []):
                raw_events.append(dict(ev, star=star, ra_deg=ra, dec_deg=dec,
                                       source_id=sid,source_kind=scope,product=path.name,reduction=family))
            rec['products'].append(result)
            if result['status']=='searched':
                families.setdefault(sid,set()).add(family)
        alternatives = []
        primary_by_source = {}
        for path,meta in sorted(source_meta.items(), key=lambda item:_rank(item[0])):
            try:
                scope = 'target'
                # Non-target released light curves must carry a different source
                # identification and position; central pointing alone is ambiguous.
                if all(v is not None for v in (visit.get('ra'),visit.get('dec'),
                                               meta.get('ra_deg'),meta.get('dec_deg'))):
                    dra = ((meta['ra_deg']-visit['ra']+180)%360-180)*math.cos(math.radians(visit['dec']))
                    sep = 3600*math.hypot(dra,meta['dec_deg']-visit['dec'])
                    if sep>8:
                        scope='field'
                if _family(path)=='DRP':
                    sid = (scope,meta.get('star'),meta.get('ra_deg'),meta.get('dec_deg')) if scope=='field' else ('target',)
                    if sid in primary_by_source:
                        rec.setdefault('auxiliary_aperture_names',[]).append(path.name)
                        if scope=='target' and any(n in path.name.lower() for n in ('rinf','rsup')):
                            alternatives.append(path)
                        continue
                    primary_by_source[sid] = path
                process(path, scope, meta)
            except Exception as exc:
                rec['errors'].append(dict(product=path.name,error=_error(exc)))
        try:
            for path in field.extract(paths, reference, folder, rec['field']):
                try:
                    process(path, 'field', photometry.metadata(path), extracted=True)
                except Exception as exc:
                    rec['errors'].append(dict(product=path.name,error=_error(exc)))
                finally:
                    path.unlink(missing_ok=True)
        except Exception as exc:
            rec['field']['status']='extraction incomplete'
            rec['field']['error']=_error(exc)
        grouped = []
        for ev in sorted(raw_events, key=lambda e: (e['photometry_pass'],e['snr']), reverse=True):
            group = next((members for members in grouped if _close(ev, members[0])
                          and all(m['product'] != ev['product'] for m in members)), None)
            if group is None:
                grouped.append([ev])
            else:
                group.append(ev)
        if any(e['source_kind']=='target' and e['photometry_pass'] for e in raw_events) and not alternatives:
            try:
                alternatives = archive_io.fetch_validation_apertures(visit, folder)
            except Exception as exc:
                rec['errors'].append(dict(stage='alternate-aperture retrieval',error=_error(exc)))
        if any(e['photometry_pass'] for e in raw_events):
            try:
                extra = archive_io.fetch_localization_images(visit, folder)
                paths = list(dict.fromkeys(paths + extra))
            except Exception as exc:
                rec['errors'].append(dict(stage='event image retrieval',error=_error(exc)))
        for group in grouped:
            rec['events'].append(_classify(group, families.get(group[0]['source_id'],set()),
                                          paths, reference, catalogs, cache, out, alternatives))
        searched = [p for p in rec['products'] if p['status']=='searched']
        rec['pipe_searched'] = any(p['reduction']=='PIPE' for p in searched)
        rec['drp_searched'] = any(p['reduction']=='DRP' for p in searched)
        if searched:
            rec['status']='searched'
            if (not rec['pipe_searched'] or not rec['drp_searched']
                    or len(searched)!=len(rec['products']) or rec['errors']
                    or not rec['field'].get('status','').startswith('catalogue apertures extracted')
                    or any(m['errors'] or m['unresolved'] for m in rec['masks'].values())):
                rec['status']='partial'
    except Exception as exc:
        rec['errors'].append(dict(stage='visit',error=_error(exc)))
        if any(p['status']=='searched' for p in rec['products']):
            rec['status']='partial'
    finally:
        manifest = folder/'manifest.jsonl'
        if manifest.exists():
            rec['receipts'] = [json.loads(line) for line in manifest.read_text().splitlines() if line]
    return rec


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inventory',type=Path)
    parser.add_argument('out',type=Path)
    parser.add_argument('cache',type=Path)
    parser.add_argument('index',type=int)
    parser.add_argument('shards',type=int)
    args=parser.parse_args()
    if not 0 <= args.index < args.shards:
        parser.error('index must be zero-based and smaller than the shard count')
    args.cache.mkdir(parents=True,exist_ok=True)
    args.out.mkdir(parents=True,exist_ok=True)
    catalog=None
    with (args.out / f'visits-{args.index:03d}.jsonl').open('w') as log, args.inventory.open() as inv:
        for i,line in enumerate(inv):
            if i % args.shards != args.index:
                continue
            if catalog is None:
                catalog=vetting.prepare_catalogs(args.cache)
            with tempfile.TemporaryDirectory(dir=args.cache) as tmp:
                rec=_run_visit(json.loads(line),Path(tmp),catalog,args.cache,
                               args.out / f'evidence-{args.index:03d}')
            log.write(json.dumps(rec,allow_nan=False,separators=(',',':'))+'\n')
            log.flush()
            shutil.rmtree(args.cache/'tess', ignore_errors=True)
            print(rec['visit_id'], rec['status'], len(rec['products']), len(rec['events']),flush=True)


if __name__ == '__main__':
    main()
