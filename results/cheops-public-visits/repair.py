"""Complete only frozen field extractions that failed on detector/window coordinates."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import sys
import tempfile
import time

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                    separators=(',',':')).encode()).hexdigest()


def open_lines(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def frozen_shards(outputs):
    files=sorted(list(outputs.glob('visits-*.jsonl.gz'))+list(outputs.glob('visits-*.jsonl')))
    if not files:
        raise ValueError('No original visit shards')
    return files


def affected(row):
    return (row.get('field',{}).get('status')=='extraction incomplete' and
            'drp centroid outside subarray' in row['field'].get('error','').lower())


def event_key(e):
    return digest([e['source_id'],e['product'],e['epoch_bjd'],e['depth'],e['duration_days']])


def disposition(e):
    vet=e.get('vetting',{})
    if vet.get('known_ephemeris_overlaps'):
        return 'catalog ephemeris overlap'
    if not e['photometry_pass'] and any(v is False for v in e['checks'].values()):
        return 'photometric rejection'
    outcomes=vet.get('checks',{})
    lit=outcomes.get('literature',{})
    if (not e['photometry_pass'] or e['localization'].get('status')!='consistent'
            or 'inconclusive' in e['cross_reduction'] or 'discrepancy' in e['cross_reduction']
            or any(r['status']!='agreement' for r in e.get('alternate_aperture_depths',{}).values())
            or outcomes.get('field',{}).get('status')!='screened'
            or outcomes.get('tess',{}).get('status') not in ('compared','no TESS sector footprint')
            or lit.get('status')!='requires publication and event epoch review' or lit.get('search_errors')
            or not lit.get('reference_data') or lit.get('truncated')
            or vet.get('unresolved') or vet.get('catalog_errors') or vet.get('errors')):
        return 'unresolved'
    return 'requires review'


def target_retry_keys(e):
    """Only unfinished services on already detected, newly localized target events."""
    vet=e.get('vetting',{})
    if (e['source_kind']!='target' or not e['photometry_pass']
            or e.get('localization',{}).get('status')!='consistent'
            or vet.get('known_ephemeris_overlaps')
            or 'inconclusive' in e['cross_reduction'] or 'discrepancy' in e['cross_reduction']
            or any(c['status']!='agreement' for c in e.get('alternate_aperture_depths',{}).values())):
        return []
    checks=vet.get('checks',{})
    keys=[]
    if checks.get('field',{}).get('status')!='screened':
        keys.append('field')
    if checks.get('tess',{}).get('status') not in ('compared','no TESS sector footprint'):
        keys.append('tess')
    lit=checks.get('literature',{})
    if not lit.get('reference_data') or lit.get('truncated'):
        keys.append('literature')
    return keys


def summarize_overlap(match, vet, overlaps):
    return dict(vet, catalog_errors=match['errors'], unresolved=match['unresolved'],
                known_objects=match['known_objects'], nearby_variables=match['nearby_variables'],
                sources=match['sources'], known_ephemeris_overlaps=overlaps,
                disposition='catalog ephemeris overlap' if overlaps else 'requires review')


def process_visit(visit, original, temp, catalogs, checkpoint, queue_event):
    import archive_io
    import field
    import photometry
    import vetting
    import worker

    checkpoint('begin', input_sha256=digest(original), prior_field=original['field'])
    log={}
    receipt_count=0
    for previous in original['events']:
        try:
            match=catalogs.match(previous['star'],previous['ra_deg'],previous['dec_deg'])
            old_set=set(previous.get('vetting',{}).get('known_ephemeris_overlaps',[]))
            added=set(vetting.known_overlaps(previous,match))-old_set
            if added:
                checkpoint('additional_overlap',key=event_key(previous),overlaps=sorted(old_set|added),
                           sources=match['sources'])
        except Exception as exc:
            checkpoint('error',stage='revised known-phase screen',key=event_key(previous),
                       error=worker._error(exc))
    try:
        paths, reference=archive_io.fetch_repair_inputs(visit,original,temp)
        manifest=temp/'manifest.jsonl'
        receipts=list(open_lines(manifest)) if manifest.exists() else []
        receipt_count=len(receipts)
        checkpoint('inputs',receipts=receipts)
        if reference is None or not paths:
            raise ValueError('Frozen image, catalogue and reference could not be recovered')
        localization_queue=[]
        seen={p['source_id'] for p in original['products'] if p['source_kind']=='field'}
        try:
            for path in field.extract(paths,reference,temp,log):
                try:
                    meta=photometry.metadata(path)
                    name,ra,dec=meta['star'],meta['ra_deg'],meta['dec_deg']
                    sid=f'{name} {ra} {dec}'
                    if sid in seen:
                        continue
                    seen.add(sid)
                    match=catalogs.match(name,ra,dec)
                    result=photometry.analyze(path,match['ephemerides'],reference_path=reference)
                    result.update(name=path.name,reduction='field aperture',source_kind='field',
                                  source_id=sid,star=name,ra_deg=ra,dec_deg=dec,
                                  derived_sha256=archive_io.digest(path))
                    new=[]
                    for trigger in result.pop('events',[]):
                        e=dict(trigger,star=name,ra_deg=ra,dec_deg=dec,source_id=sid,source_kind='field',
                               product=path.name,reduction='field aperture',cross_reduction='single reduction family',
                               alternate_aperture_depths={},localization={'status':'pending image check'})
                        e['matching_products']=[{k:e[k] for k in ('product','reduction','epoch_bjd','duration_days',
                                                                   'depth','snr','checks')}]
                        e['vetting']=summarize_overlap(match,{'checks':{},'errors':['external checks pending']},
                                                       vetting.known_overlaps(e,match))
                        e['disposition']=disposition(e)
                        new.append(e)
                    # Preserve all extracted curves and crossings before optional external work.
                    checkpoint('source', mask=match,product=result,events=new)
                    localization_queue.extend(new)
                except Exception as exc:
                    checkpoint('error',stage='field source',product=path.name,error=worker._error(exc))
                finally:
                    path.unlink(missing_ok=True)
        except Exception as exc:
            log.update(status='repair incomplete',error=worker._error(exc))
        # Extract/search every source before reopening the cube once for localization.
        # A previously rejected signal needs no further pixel interpretation when a
        # measured veto or known ephemeris already identifies the exclusion.
        original_survivors=[e for e in original['events'] if e.get('photometry_pass')]
        if localization_queue or original_survivors:
            from contextlib import ExitStack
            with ExitStack() as stack:
                try:
                    locate=stack.enter_context(field.locator(paths,reference))
                except Exception as exc:
                    reason=worker._error(exc)
                    locate=lambda event,coords: dict(status='incomplete',error=reason)
                for old in original_survivors:
                    try:
                        coords=(old['ra_deg'],old['dec_deg'])
                        match=catalogs.match(old['star'],*coords)
                        old_v=old.get('vetting',{})
                        overlap=vetting.known_overlaps(old,match)
                        v=summarize_overlap(match,old_v,list(dict.fromkeys(
                            old_v.get('known_ephemeris_overlaps',[])+overlap)))
                        try:
                            loc=locate(old,coords)
                        except Exception as exc:
                            loc=dict(status='incomplete',error=worker._error(exc))
                        updated=dict(old,localization=loc,vetting=v)
                        checkpoint('target_evidence',key=event_key(old),localization=loc,
                                   vetting={k:v[k] for k in ('catalog_errors','unresolved','known_objects',
                                        'nearby_variables','sources','known_ephemeris_overlaps','disposition')},
                                   disposition=disposition(updated))
                        if target_retry_keys(updated):
                            queue_event(visit['visit_id'],updated)
                    except Exception as exc:
                        checkpoint('error',stage='existing photometric event image or binary check',
                                   event_key=event_key(old),error=worker._error(exc))
                for e in localization_queue:
                    if e['disposition'] in ('catalog ephemeris overlap','photometric rejection'):
                        failed=[k for k,v in e['checks'].items() if v is False]
                        e['localization']=dict(status='excluded before pixel checking',
                            reason=e['disposition'],failed_photometric_checks=failed,
                            known_overlaps=e['vetting'].get('known_ephemeris_overlaps',[]))
                    else:
                        try:
                            e['localization']=locate(e,(e['ra_deg'],e['dec_deg']))
                        except Exception as exc:
                            e['localization']=dict(status='incomplete',error=worker._error(exc))
                    e['disposition']=disposition(e)
                    checkpoint('image',key=event_key(e),localization=e['localization'],
                               disposition=e['disposition'])
                    queue_event(visit['visit_id'],e)
    except Exception as exc:
        log.update(status='repair incomplete',error=worker._error(exc))
    finally:
        manifest=temp/'manifest.jsonl'
        receipts=list(open_lines(manifest)) if manifest.exists() else []
        checkpoint('end',field=log,receipts=receipts[receipt_count:])


def external_checks(queue_path,catalogs,cache,out,write,deadline):
    import vetting
    import worker

    # Prioritize coherent events. Process every queued event that fits the same run.
    class ServiceBudget(BaseException):
        pass
    def timeout_alarm(signum,frame):
        raise ServiceBudget('service time budget reached')
    signal.signal(signal.SIGALRM,timeout_alarm)
    for priority in (0,1,2):
        for q in open_lines(queue_path):
            e=q['event']
            existing_target=e['source_kind']=='target'
            actual_priority=0 if existing_target else (1 if e['photometry_pass'] else 2)
            if actual_priority!=priority:
                continue
            v=e['vetting']
            actions=[('field',vetting._gaia,45),('literature',vetting._literature,70),
                     ('tess',vetting._tess,120)]
            eligible=set(target_retry_keys(e)) if existing_target else {a[0] for a in actions}
            v.setdefault('checks',{})
            for key,action,seconds in actions:
                if key not in eligible:
                    continue
                remain=int(deadline-time.monotonic())
                if remain<5:
                    return
                try:
                    signal.alarm(min(seconds,remain))
                    v['checks'][key]=action(e,cache)
                    if key=='literature':
                        worker._write_refs(v['checks'][key],out)
                except (Exception,ServiceBudget) as exc:
                    v['checks'][key]=dict(status='incomplete',error=worker._error(exc))
                finally:
                    signal.alarm(0)
                missing=[k for k,_,_ in actions if k not in v['checks']]
                v['errors']=[f'{k} check pending' for k in missing]
                if existing_target:
                    write(q['visit_id'],'target_external',key=event_key(e),service=key,
                          result=v['checks'][key],errors=v['errors'],disposition=disposition(e))
                else:
                    write(q['visit_id'],'external',key=event_key(e),vetting=v,disposition=disposition(e))
            shutil.rmtree(cache/'tess',ignore_errors=True)



def run(args):
    import vetting
    source=args.study/'outputs'
    args.out.mkdir(parents=True,exist_ok=True)
    args.cache.mkdir(parents=True,exist_ok=True)
    if not 0<=args.index<args.shards:
        raise ValueError('Invalid shard assignment')
    idx={r['visit_id']:r for r in open_lines(source/'inventory.jsonl')}
    catalogs=vetting.prepare_catalogs(args.cache)
    deadline=time.monotonic()+int(os.getenv('CHEOPS_REPAIR_SECONDS','19800'))
    queue=args.cache/f'field-event-queue-{args.index:03d}.jsonl'
    with (args.out/f'repair-{args.index:03d}.jsonl').open('w') as log:
        def write(visit,op,**vals):
            log.write(json.dumps(dict(visit_id=visit,op=op,**vals),allow_nan=False,
                                 separators=(',',':'))+'\n')
            log.flush()
        with queue.open('w') as todo:
            def enqueue(vid,event):
                todo.write(json.dumps(dict(visit_id=vid,event=event),allow_nan=False,
                                      separators=(',',':'))+'\n')
            for shard in frozen_shards(source):
                shard_id=int(shard.name.partition('-')[2].partition('.')[0])
                if shard_id % args.shards!=args.index:
                    continue
                for row in open_lines(shard):
                    if not affected(row):
                        continue
                    if time.monotonic()>=deadline:
                        break
                    vid=row['visit_id']
                    with tempfile.TemporaryDirectory(dir=args.cache) as td:
                        process_visit(idx[vid],row,Path(td),catalogs,
                                      lambda op,**kw:write(vid,op,**kw),enqueue)
                    print(vid,'field checkpoint complete',flush=True)
        external_checks(queue,catalogs,args.cache,args.out/f'repair-evidence-{args.index:03d}',
                        write,deadline)


def apply_patch(row, ops):
    if not affected(row) or ops[0]['op']!='begin' or ops[0]['input_sha256']!=digest(row):
        raise ValueError('Repair does not match the frozen visit '+row['visit_id'])
    original_products=digest(row['products'])
    nproducts=len(row['products'])
    original_measurements=[event_key(e) for e in row['events']]
    events={event_key(e):e for e in row['events']}
    target_allowed={}
    source_ids={p['source_id'] for p in row['products'] if p['source_kind']=='field'}
    repair=dict(prior_record_sha256=ops[0]['input_sha256'],prior_field=row['field'],
                prior_status=row['status'],prior_events={},added_products=0,errors=[])
    row['field_repair']=repair
    row['field']={'status':'repair incomplete; worker checkpoint ended before completion'}
    for op in ops[1:]:
        kind=op['op']
        if kind=='source':
            p=op['product']
            if p['source_kind']!='field' or p['source_id'] in source_ids:
                raise ValueError('Repair attempted to replace an existing source')
            source_ids.add(p['source_id'])
            row['products'].append(p)
            row['masks'][p['source_id']]=op['mask']
            repair['added_products']+=1
            for e in op['events']:
                key=event_key(e)
                if key in events or e['source_id']!=p['source_id']:
                    raise ValueError('Inconsistent repair event ID')
                events[key]=e
                row['events'].append(e)
        elif kind=='additional_overlap':
            e=events[op['key']]
            if op['key'] not in original_measurements:
                raise ValueError('Updated known-phase screen must apply to an original event')
            repair.setdefault('revised_phase_checks',{})[op['key']]=dict(
                prior_disposition=e['disposition'],prior_overlaps=e['vetting'].get('known_ephemeris_overlaps',[]),
                sources=op['sources'])
            e['vetting']['known_ephemeris_overlaps']=op['overlaps']
            e['vetting']['disposition']='catalog ephemeris overlap'
            e['disposition']='catalog ephemeris overlap'
        elif kind=='target_evidence':
            e=events[op['key']]
            if not e.get('photometry_pass') or op['key'] not in original_measurements:
                raise ValueError('Unrequested existing-event change')
            repair['prior_events'][op['key']]=dict(localization=e.get('localization'),
                disposition=e.get('disposition'),vetting={k:e.get('vetting',{}).get(k) for k in op['vetting']})
            e.update(localization=op['localization'],disposition=op['disposition'])
            e['vetting'].update(op['vetting'])
            target_allowed[op['key']]=set(target_retry_keys(e))
        elif kind=='target_external':
            key, service=op['key'],op['service']
            if service not in target_allowed.get(key,set()):
                raise ValueError('Target service update lacks an eligible corrected pixel check')
            target_allowed[key].remove(service)
            e=events[key]
            v=e['vetting']
            history=repair['prior_events'][key]
            history.setdefault('external_errors',v.get('errors',[]))
            history.setdefault('external_checks',{})[service]=v.setdefault('checks',{}).get(service)
            v['checks'][service]=op['result']
            v['errors']=op['errors']
            e['disposition']=disposition(e)
        elif kind in ('image','external'):
            if op['key'] in original_measurements:
                raise ValueError('External lookup attempted to replace a prior event')
            e=events[op['key']]
            key='vetting' if kind=='external' else 'localization'
            e.update({key:op[key],'disposition':op['disposition']})
        elif kind in ('inputs','end'):
            repair['added_receipts']=repair.get('added_receipts',0)+len(op['receipts'])
            row['receipts'].extend(op['receipts'])
            if kind=='end':
                row['field']=op['field']
        elif kind=='error':
            repair['errors'].append({k:v for k,v in op.items() if k not in ('op','visit_id')})
        else:
            raise ValueError('Unknown repair operation '+kind)
    if digest(row['products'][:nproducts])!=original_products:
        raise ValueError('Original products or injection trials changed')
    if [event_key(e) for e in row['events'][:len(original_measurements)]]!=original_measurements:
        raise ValueError('Original event measurements changed')
    # Restore truthful coverage status. Original errors (including missing PIPE) persist.
    if any(p['status']=='searched' for p in row['products']):
        row['status']='partial'
        if (row.get('drp_searched') and row.get('pipe_searched') and not row['errors']
                and not repair['errors'] and not any(m['errors'] for m in row['masks'].values())
                and all(p['status']=='searched' for p in row['products'])
                and row['field'].get('status','').startswith('catalogue apertures extracted')):
            row['status']='searched'
    return row


def merge(study, patches):
    output=study/'outputs'
    # SQLite keeps the per-visit deltas indexed on disk for large field catalogues.
    import sqlite3
    with tempfile.TemporaryDirectory() as td, sqlite3.connect(str(Path(td)/'repair.sqlite')) as db:
        db.execute('CREATE TABLE delta (vid TEXT, seq INTEGER PRIMARY KEY, value TEXT)')
        for path in sorted(list(patches.rglob('repair-*.jsonl'))+list(patches.rglob('repair-*.jsonl.gz'))):
            db.executemany('INSERT INTO delta (vid,value) VALUES (?,?)',
                           ((r['visit_id'],json.dumps(r,separators=(',',':'))) for r in open_lines(path)))
        db.execute('CREATE INDEX visit_delta ON delta(vid)')
        count=0
        for shard in frozen_shards(output):
            tmp=shard.with_name(shard.name+'.tmp')
            changed=False
            with tmp.open('wb') as fp:
                stream=gzip.GzipFile(filename='',fileobj=fp,mode='wb',mtime=0) if shard.suffix=='.gz' else fp
                try:
                    for row in open_lines(shard):
                        vid=row['visit_id']
                        ops=[json.loads(r[0]) for r in db.execute(
                             'SELECT value FROM delta WHERE vid=? ORDER BY seq',(vid,))]
                        if ops:
                            row=apply_patch(row,ops)
                            count+=1
                            changed=True
                            db.execute('DELETE FROM delta WHERE vid=?',(vid,))
                        stream.write((json.dumps(row,separators=(',',':'),allow_nan=False)+'\n').encode())
                finally:
                    if stream is not fp:
                        stream.close()
            if changed:
                tmp.replace(shard)
            else:
                tmp.unlink()
        if db.execute('SELECT COUNT(*) FROM delta').fetchone()[0]:
            raise ValueError('Repair mentions a visit absent from frozen input')
    for path in patches.rglob('repair-evidence-*'):
        if path.is_dir():
            shutil.copytree(path,output/path.name,dirs_exist_ok=True)
    print(f'Merged {count} field repair checkpoints; original product measurements preserved')


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='merge':
        if len(sys.argv)!=4:
            raise SystemExit('usage: repair.py merge STUDY_DIR PATCH_DIR')
        merge(Path(sys.argv[2]),Path(sys.argv[3]))
    else:
        parser=argparse.ArgumentParser(description=__doc__)
        for name in ('study','out','cache'):
            parser.add_argument(name,type=Path)
        parser.add_argument('index',type=int)
        parser.add_argument('shards',type=int)
        run(parser.parse_args())
