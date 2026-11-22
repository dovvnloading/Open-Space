"""Summarize the frozen CHEOPS visit records and optional reviewed event dispositions."""
import csv
import gzip
import hashlib
import os
import json
import io
from collections import Counter, defaultdict
from pathlib import Path

import publish

OUT = Path(os.environ.get('CHEOPS_OUTPUTS',publish.OUT)).resolve()
publish.OUT = OUT
STUDY = publish.STUDY
REVIEW_CHECKS = ('known_planet', 'known_eclipse', 'roll_angle', 'cosmic_ray',
                 'background_star', 'tess', 'literature')


def lines(path):
    with (gzip.open(path,'rt',encoding='utf8') if path.suffix=='.gz' else path.open(encoding='utf8')) as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def cell(value):
    return str(value).replace('|', '/').replace('\n', ' ')


def event_id(visit, event):
    identity = [visit, event['source_id'], event['product'], event['epoch_bjd'],
                event['depth'], event['duration_days']]
    return hashlib.sha256(json.dumps(identity, separators=(',', ':')).encode()).hexdigest()[:24]


def decisions():
    path = STUDY / 'review_decisions.json'
    if not path.exists():
        return {}
    review = json.loads(path.read_text())
    if not isinstance(review, dict):
        raise ValueError('Expected review decisions keyed by event ID')
    for key, item in review.items():
        if item['decision'] not in ('candidate','rejected') or not str(item['reason']).strip():
            raise ValueError(f'Invalid review decision {key}')
        if not item.get('references') or not all(isinstance(v,str) and v.strip()
                                                 for v in item['references']):
            raise ValueError(f'Missing review evidence {key}')
        if item['decision']=='candidate' and not all(
                isinstance(item.get('checks',{}).get(k),str) and item['checks'][k].strip()
                for k in REVIEW_CHECKS):
            raise ValueError(f'Missing required candidate check {key}')
    return review


def event_summary(visit, e, review):
    key = event_id(visit, e)
    decision = review.get(key)
    stage = e['disposition']
    if decision:
        if decision['decision']=='candidate' and stage!='requires review':
            raise ValueError(f'Candidate {key} failed the automated evidence gates')
        stage = 'reviewed candidate' if decision['decision']=='candidate' else 'reviewed rejection'
    vet = e.get('vetting', {})
    checks = vet.get('checks', {})
    return dict(event_id=key, visit_id=visit, star=e.get('star',''), source_kind=e['source_kind'],
                source_id=e['source_id'], epoch_bjd_tdb=e['epoch_bjd'], depth_ppm=1e6*e['depth'],
                duration_hours=24*e['duration_days'], snr=e['snr'], disposition=stage,
                product=e['product'], corroboration=e['cross_reduction'],
                automated_photometry_pass=e['photometry_pass'],
                photometry_checks=json.dumps(e['checks'],sort_keys=True),
                image_localization=e.get('localization',{}).get('status','missing'),
                gaia=checks.get('field',{}).get('status','missing'),
                tess=checks.get('tess',{}).get('status','missing'),
                literature=checks.get('literature',{}).get('status','missing'),
                known_overlaps=' / '.join(vet.get('known_ephemeris_overlaps',[])),
                review_reason=decision['reason'] if decision else '',
                review_checks=json.dumps(decision.get('checks',{}),sort_keys=True) if decision else '',
                review_references=' '.join(decision['references']) if decision else '')


def preferred(product):
    n = product['name'].lower()
    return (not ('default' in n), not ('optimal' in n), n)


def render():
    status = json.loads((OUT/'archive_status.json').read_text())
    if not status['verified_lightcurve']:
        return publish.publish()
    inventory = {r['visit_id'] for r in lines(OUT/'inventory.jsonl')}
    if len(inventory)!=status['visits']:
        raise ValueError('Inventory count differs from access receipt')
    review = decisions()
    done, searched, reviewed_ids = set(), set(), set()
    products, visit_status, limitations, event_counts = Counter(), Counter(), Counter(), Counter()
    per_family, source_visits = defaultdict(set), set()
    injections = defaultdict(Counter)
    events, candidates, pending = [], [], []
    bounds, sectors = [], Counter()
    field_limits = Counter()
    photometry_exclusions = Counter()
    repair_counts = Counter()
    distinct_stars = defaultdict(set)
    seen_digests = set()
    with (OUT/'input_checksums.jsonl.gz').open('wb') as raw, gzip.GzipFile(
            filename='', fileobj=raw, mode='wb', mtime=0) as comp, io.TextIOWrapper(comp, encoding='utf8') as inputs:
        def receipt(visit, name, sha, context):
            if sha and sha not in seen_digests:
                seen_digests.add(sha)
                inputs.write(json.dumps(dict(visit_id=visit, name=name, sha256=sha, source=context),
                                        sort_keys=True, separators=(',',':'))+'\n')
        for p in lines(OUT/'archive_probes.jsonl'):
            receipt(p.get('file_key',''),p.get('file',p.get('operation')),p.get('sha256'),p.get('url','archive'))
        for shard in sorted(list(OUT.rglob('visits-*.jsonl'))+list(OUT.rglob('visits-*.jsonl.gz'))):
            for row in lines(shard):
                visit = row['visit_id']
                if visit not in inventory or visit in done:
                    raise ValueError(f'Unexpected or duplicate visit record {visit}')
                done.add(visit)
                visit_status[row['status']] += 1
                for r in row.get('receipts',[]):
                    receipt(visit,r.get('file',r.get('operation')),r.get('sha256'),r.get('url','archive'))
                if row.get('errors'):
                    limitations['Visit or product interpretation error'] += 1
                field_status = row.get('field',{}).get('status','missing field products')
                if not field_status.startswith('catalogue apertures extracted'):
                    limitations['Field photometry unavailable or incomplete'] += 1
                if not row.get('pipe_searched'):
                    limitations['No searched PIPE light curve'] += 1
                if not row.get('drp_searched'):
                    limitations['No searched DRP light curve'] += 1
                for sid,m in row.get('masks',{}).items():
                    for k,src in m.get('sources',{}).items():
                        receipt(visit,k,src.get('sha256'),src.get('url',k))
                    if m.get('errors') or m.get('unresolved'):
                        limitations['Unresolved catalog query or known ephemeris (source-visits)'] += 1
                diag = row.get('field', {})
                for name in ('near_target','outside','too_few','extracted','unresolved_core','quality_excluded','undetected','blended'):
                    if isinstance(diag.get(name), int):
                        field_limits[name] += diag[name]
                repair_counts[diag.get('status', 'no field status')] += 1
                selected = {}
                for p in row['products']:
                    if p.get('derived_sha256'):
                        receipt(visit,p['name'],p['derived_sha256'],'field aperture extraction')
                    products[p['status']] += 1
                    if p['status']!='searched':
                        photometry_exclusions[(p.get('source_kind','unknown'),p['status'])] += 1
                        continue
                    distinct_stars[p['source_kind']].add((p.get('star') or str(p['source_id'])).strip())
                    searched.add(visit)
                    family = p['reduction']
                    scope = p['source_kind']
                    per_family[(scope,family)].add((visit,p['source_id']))
                    source_visits.add((visit,p['source_id']))
                    key = (scope,family,p['source_id'])
                    if key not in selected or preferred(p)<preferred(selected[key]):
                        selected[key]=p
                    if p.get('span_bjd'):
                        bounds.append(p['span_bjd'])
                for (scope,family,sid),product in selected.items():
                    for b in product['injections']:
                        c = injections[(scope,family,b['depth_ppm'],b['duration_hours'])]
                        for k in ('trials','sampled','confused','recovered'):
                            c[k] += b[k]
                for ev in row['events']:
                    checks=ev.get('vetting',{}).get('checks',{})
                    external=[checks.get('field',{}).get('source',{})]
                    external+=checks.get('literature',{}).get('sources',[])
                    external+=checks.get('tess',{}).get('products',[])
                    for data in external:
                        receipt(visit,data.get('file',data.get('uri','external query')),
                                data.get('sha256'),data.get('url',data.get('uri','vetting')))
                    item = event_summary(visit,ev,review)
                    events.append(item)
                    stage = item['disposition']
                    event_counts[stage] += 1
                    if item['event_id'] in review:
                        reviewed_ids.add(item['event_id'])
                    if stage=='reviewed candidate':
                        candidates.append(item)
                    elif stage in ('requires review','unresolved'):
                        pending.append(item)
                        sectors['Image localization: '+item['image_localization']] += 1
                        sectors['Comparison between reductions: '+item['corroboration']] += 1
                    elif stage not in ('photometric rejection','catalog ephemeris overlap','reviewed rejection'):
                        raise ValueError(f'Unrecognized event disposition {stage}')
    if set(review) != reviewed_ids:
        raise ValueError('Review file refers to events absent from this evidence snapshot')
    columns = list(events[0]) if events else ['event_id','visit_id','star','source_kind','source_id',
                'epoch_bjd_tdb','depth_ppm','duration_hours','snr','disposition','product',
                'corroboration','automated_photometry_pass','photometry_checks','image_localization',
                'gaia','tess','literature','known_overlaps','review_reason','review_checks','review_references']
    with (OUT/'events.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t')
        w.writeheader()
        w.writerows(sorted(events,key=lambda x:(x['visit_id'],x['epoch_bjd_tdb'],str(x['star']))))
    unprocessed=sorted(inventory-done)
    (OUT/'missing_visits.txt').write_text(''.join(x+'\n' for x in unprocessed))
    complete = bool(inventory) and done==inventory and bool(status.get('inventory_complete'))
    bounded = complete and bool(searched) and not candidates and any(c['trials'] for c in injections.values())
    if candidates:
        title=f"{len(candidates):,} CHEOPS dimming event{'s' if len(candidates)!=1 else ''} passed candidate review"
        label='Candidate'
        catalog_text=f"**{len(candidates):,} CHEOPS dimming event{'s' if len(candidates)!=1 else ''} passed the archive transit search and published-event checks.**"
    elif bounded:
        title=f"No new transit candidate survived in {len(searched):,} CHEOPS visits"
        label='Null result'
        catalog_text=f"**No new transit passed every check in {len(searched):,} searchable CHEOPS visits.**"
    elif searched:
        title=f"The CHEOPS search covered {len(searched):,} visits"
        label='Measured'
        catalog_text=f"**The transit search measured sensitivity in {len(searched):,} CHEOPS visits.**"
    else:
        title='No CHEOPS visit reached the transit search'
        label='Null result'
        catalog_text='**The CHEOPS data checks yielded zero searchable visits in this run.** No transit sensitivity was measured.'
    if pending:
        catalog_text += f" {len(pending):,} signals remain unresolved."
    elif searched and not complete and not candidates:
        catalog_text += ' Incomplete coverage limits the result.'
    published = [f'# {title}', '',
        f'{len(candidates):,} new candidates passed automated checks and publication review. '
        f'{len(pending):,} additional signals require evidence or review. '
        f'The inventory contained {len(inventory):,} public visit identifiers. '
        f'The workers returned {len(done):,} visit records and searched {len(searched):,} visits '
        f'covering {len(source_visits):,} distinct visit/source pairs.']
    if bounded:
        published.append('This is a null for validated candidates in the searched light curves. '
                         'The injection table measures recovery at the photometric stage. '
                         'Open signals and incomplete checks prevent a limit on the number of planets or eclipsing systems.')
    elif not candidates:
        published.append('The available observations and unresolved checks do not establish an archive-wide astrophysical null.')
    if candidates or pending:
        published += ['', '## Events', '',
            '[The event table](outputs/events.tsv) gives the star, visit, BJD TDB epoch, fractional-depth equivalent in ppm, duration, reduction checks, image result, TESS check and disposition for every distinct signal. Events from multiple apertures at the same source and epoch are grouped.']
        showcase = sorted(candidates,key=lambda e:(e['star'],e['epoch_bjd_tdb']))
        if not showcase and len(pending)<=25:
            showcase=pending
        if showcase:
            published += ['', '| Star | Visit | Epoch (BJD TDB) | Depth (ppm) | Duration (h) | Status |',
                           '| --- | --- | ---: | ---: | ---: | --- |']
            for item in showcase:
                published.append(f"| {cell(item['star'])} | {cell(item['visit_id'])} | {item['epoch_bjd_tdb']:.6f} | {item['depth_ppm']:.0f} | {item['duration_hours']:.2f} | {cell(item['disposition'])} |")
        if candidates:
            for e in candidates:
                published += ['',f"**{cell(e['star'])}, BJD {e['epoch_bjd_tdb']:.6f}.** {e['review_reason']}",
                    ' '.join(f'{k.replace("_"," ")} — {v}.' for k,v in json.loads(e['review_checks']).items()),
                    e['review_references']]
    published += ['', '## Coverage and exclusions', '',
        '| Extraction | Distinct visit/source pairs |','| --- | ---: |']
    for (scope,fam),pairs in sorted(per_family.items()):
        published.append(f'| {scope}, {fam} | {len(pairs):,} |')
    if distinct_stars:
        published += ['', 'The searched light curves have the following distinct source identifiers, '+', '.join(
            f'{len(names):,} {kind}' for kind,names in sorted(distinct_stars.items()))+'. '
            'Target aliases can refer to the same star across visits.']
    if repair_counts:
        published += ['', '| Field extraction outcome | Visits |','| --- | ---: |']
        published += [f'| {cell(k)} | {n:,} |' for k,n in repair_counts.most_common()]
    if field_limits:
        published += ['', 'Across visit catalogues, '+', '.join(f'{k.replace("_"," ")} {v:,}' for k,v in field_limits.items())+'.']
    if photometry_exclusions:
        published += ['', '| Product excluded from search | Count |', '| --- | ---: |']
        published += [f'| {k[0]}, {k[1]} | {v:,} |' for k,v in photometry_exclusions.most_common()]

    if not per_family:
        published.append('| Searched | 0 |')
    if bounds:
        published+=['',f'The first and last accepted product spans run from BJD TDB {min(b[0] for b in bounds):.6f} to {max(b[1] for b in bounds):.6f}. These endpoints include gaps between and within visits.']
    published+=['',f'{len(unprocessed):,} inventory visits have no worker record. '
                    f'{len(done-searched):,} recorded visits have no searched light curve. '
                    f'{visit_status["partial"]:,} visits have partial product, field or ephemeris coverage. '
                    f'DACE listing completeness was {status.get("dace_complete")}. ESA listing completeness was {status.get("esa_complete")}.',
                '', '| Search limitation | Count |', '| --- | ---: |']
    for reason,count in limitations.most_common():
        published.append(f'| {cell(reason)} | {count:,} |')
    if not limitations:
        published.append('| Recorded visit-level limitations | 0 |')
    if event_counts:
        published+=['','| Signal disposition | Count |','| --- | ---: |']
        published += [f'| {k} | {n:,} |' for k,n in sorted(event_counts.items())]
    if sectors:
        published+=['','Unresolved signals can have multiple evidence limitations.','','| Unresolved signal check | Count |','| --- | ---: |']
        published += [f'| {cell(k)} | {n:,} |' for k,n in sectors.most_common()]
    published+=['','## Injection recovery','',
        'Each selected visit/source product received trapezoidal dips before detrending, at epochs drawn uniformly across its quality-filtered visit span, including gaps and known-event masks. Two trials were drawn for each depth and duration in each product. Every trial counts in the denominator. Recovery requires a consistent epoch and duration and all automatic photometric cuts. The table reports recovered/tried and the observed recovery fraction.',
        '',
        'Each source contributes one product per visit and reduction family, preferring DEFAULT, then OPTIMAL, then filename order. DRP and PIPE can cover the same photons. Curves from distinct field stars in the same visit can share noise. Injections into the extracted field-star curves measure detectability after aperture extraction. A field star needed at least 48 valid cadences and a measurable core in the sampled frames. Field apertures used a radius of 10 pixels and a local sky annulus of 14 to 18 pixels. The smear and contamination regressors come from the visit’s DRP target telemetry. The recovery fractions measure the automated photometric stage. Catalog completeness, image localization and publication review are separate filters.']
    durations=sorted({k[3] for k in injections})
    if durations:
        published += ['', '| Source and extraction | Depth (ppm) | '+' | '.join(f'{d:g} h' for d in durations)+' |',
                      '| --- | ---: | '+' | '.join('---:' for _ in durations)+' |']
        for scope,fam,depth in sorted({k[:3] for k in injections}):
            vals=[]
            for dur in durations:
                c=injections.get((scope,fam,depth,dur),{})
                total=c.get('trials',0)
                vals.append(f"{c['recovered']:,}/{total:,} ({c['recovered']/total:.1%})" if total else '—')
            published.append(f"| {scope}, {fam} | {depth:,} | "+' | '.join(vals)+' |')
    else:
        published += ['', 'Zero injection trials produced a usable recovery denominator.']
    published+=['','## Search and reproduction','',
        'The visit search masked published transit and eclipse windows, then fit a time trend, four Fourier harmonics of spacecraft roll, background, smear and contamination with iterative outlier rejection. A sliding box compared its central depth against both adjacent flanks. Durations ranged from 0.5 to 8 hours, with a nominal threshold of 8 times the larger of the white-noise floor and the empirical scatter of the trial depths at that duration. The search required at least four cadences in each interval, including both flanks.',
        '',
        'The signal checks tested persistence after deleting the deepest cadence, baseline agreement, recurrence at the same roll angle and auxiliary changes. Independent extraction and pixel-difference centroid checks, nearby variables, Gaia neighbors, NASA planet and TOI ephemerides, available TESS data and the literature constrain identification. The recorded depth and duration describe the best box. The sampling cuts limit sensitivity at the edges of a visit and across data gaps. Field apertures require calibrated images, catalog positions, enough valid cadences and separation from the target. Unresolved blends, catalog omissions, faint uncataloged sources and signals outside the searched durations remain possible.',
        '',
        f"The frozen visit list has SHA-256 `{status['inventory_sha256']}`. "
        '[Input digests](outputs/input_checksums.jsonl.gz) identify the retrieved archive responses, source masks, light curves and derived apertures. The visit records carry request URLs, product checksums, masks, diagnostics and per-event evidence. '
        '[Output checksums](outputs/SHA256SUMS) cover every deposited output. '
        '[Missing visit IDs](outputs/missing_visits.txt) and [all distinct events](outputs/events.tsv) are machine-readable.',
        '', '## Reproduce the result', '',
        'The [first run](https://github.com/dovvnloading/Open-Space/actions/runs/45054151730) deposited the frozen inventory. The [target survey](https://github.com/dovvnloading/Open-Space/actions/runs/45054699737) contains the original target measurements and catalog snapshot. The final field pass translated detector centroids into subarray coordinates and retained the original target measurements and injection trials. The public archive and external catalogs can change between downloads.',
        '', 'To regenerate the report from the deposited measurements in a clean checkout, run', '',
        '```sh','study=results/'+STUDY.name,'(cd "$study/outputs" && sha256sum -c SHA256SUMS)',
        'python "$study/science_report.py"','```','',
        'To repeat the field pass against the same target measurements, use the following commands from the repository root. The base commit fixes the original visit ledgers. The final catalog seed fixes the published masks used for field photometry.', '',
        '```sh','study=results/'+STUDY.name,'scratch=$(mktemp -d)',
        'git archive 6f0b77d2c85cb5e422d12e81e1c698f6f9f5ed92 "$study" | tar -x -C "$scratch"',
        'base="$scratch/$study"','python -m pip install -r "$study/requirements.txt"',
        'export CHEOPS_CATALOG_SEED="$PWD/$study/outputs/field-ephemerides/catalog-seed.json.gz"',
        'python "$study/repair.py" "$base" "$scratch/patches" "$scratch/cache" 0 1',
        'python "$study/repair.py" merge "$base" "$scratch/patches"',
        'CHEOPS_OUTPUTS="$base/outputs" python "$study/science_report.py"','```','',
        'The workflow at commit `69fb1aa8c973b9ef1d749a2dbef5cdede3682a9d` provides the commands for the original DRP survey. The input receipts identify the actual bytes for comparison.']
    replay = 'CHEOPS_OUTPUTS' in os.environ
    destination = OUT/'README.md' if replay else STUDY/'README.md'
    text = '\n'.join(published)+'\n'
    destination.write_text(text.replace('](outputs/', '](') if replay else text)
    if not replay:
        publish.catalog(catalog_text, label)
    publish.checksums()


if __name__=='__main__':
    render()
