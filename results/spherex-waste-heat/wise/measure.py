"""Measure WISE counterparts and per-visit statistics from the archived CSVs."""
import argparse
import csv
import hashlib
from pathlib import Path
import gzip
import io
import itertools
import json
import math
from statistics import median
from fetch import ROOT, TABLES, position, separation, targets, request


def read(target, kind, archive, allow_legacy_status=False):
    key=target['label']+'_'+kind
    raw=gzip.decompress((archive/(key+'.csv.gz')).read_bytes())
    rec=json.loads((archive/(key+'.json')).read_text(encoding='utf-8'))
    query,url=request(kind,target)
    reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    rows=list(reader)
    status_ok=(rec.get('status')==200 if 'status' in rec else allow_legacy_status)
    if ('error' in rec or not status_ok or rec.get('url')!=url or rec.get('query')!=query
        or len(raw)!=rec['bytes'] or hashlib.sha256(raw).hexdigest()!=rec['sha256']
        or reader.fieldnames!=TABLES[kind][1].split(',') or len(rows)!=rec.get('rows')):
        raise ValueError(f'Archived response mismatch: {key}')
    return rows


def write(output, name, rows, fields=None):
    with (output/name).open('x', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields or list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def angular_distance(row,target,epoch):
    return separation(float(row['ra']),float(row['dec']),*position(target,epoch))


def arc(row, target, kind):
    if kind == 'catwise' and row['meanobsmjd']:
        epoch=2000+(float(row['meanobsmjd'])-51544.5)/365.25
    else:
        epoch=TABLES[kind][2]
    return angular_distance(row,target,epoch)


def summarize_catalogs(target,archive,allow_legacy_status=False):
    fields=dict(label=target['label'],gaia_source_id=target['source_id'])
    for kind in ('catwise','unwise','allwise'):
        rows=read(target,kind,archive,allow_legacy_status)
        ranked=sorted(rows,key=lambda r:arc(r,target,kind))
        pref=kind+'_'
        fields[pref+'neighbors_12arcsec']=len(rows)
        fields[pref+'next_sep_as']=round(arc(ranked[1],target,kind),3) if len(rows)>1 else ''
        row=ranked[0] if ranked and arc(ranked[0],target,kind)<=2 else {}
        fields[pref+'sep_as']=round(arc(row,target,kind),3) if row else ''
        if kind=='unwise':
            fields[pref+'id']=row.get('unwise_objid','')
            for i,cor in ((1,0.004),(2,0.032)):
                flux=float(row.get('flux_'+str(i)) or 0)
                fields[pref+f'w{i}']=round(22.5-2.5*math.log10(flux)-cor,4) if flux>0 else ''
                fields[pref+f'w{i}_staterr']=round(1.085736*float(row['dflux_'+str(i)])/flux,4) if flux>0 else ''
                fields[pref+f'w{i}_fracflux']=row.get('fracflux_'+str(i),'')
                fields[pref+f'w{i}_flags']=row.get('flags_unwise_'+str(i),'')
        else:
            cols=['source_name','cc_flags','ab_flags'] if kind=='catwise' else ['designation','cc_flags','ph_qual','ext_flg']
            for c in cols:
                fields[pref+c]=row.get(c,'')
            bands=2 if kind=='catwise' else 4
            for i in range(1,bands+1):
                for field in ['mpro','sigmpro']+(['snr'] if kind=='allwise' else []):
                    c=f'w{i}'+field
                    fields[pref+c]=row.get(c,'')
    return fields


def single_exposures(target, band,archive,allow_legacy_status=False):
    selected={}
    for r in read(target,'neo',archive,allow_legacy_status):
        mjd=float(r['mjd'])
        epoch=2000+(mjd-51544.5)/365.25
        sep=angular_distance(r,target,epoch)
        if sep>2 or float(r['qual_frame'])<1 or float(r['qi_fact'])<=0 or float(r['saa_sep'])<=0:
            continue
        if r['cc_flags'][band-1:band]!='0' or r['moon_masked'][band-1:band]!='0':
            continue
        key=f'w{band}'
        if r['ph_qual'][band-1:band] not in ('A','B') or not r[key+'sigmpro']:
            continue
        if float(r[key+'sat'] or 0)>0 or float(r[key+'rchi2'] or 'inf')>=3:
            continue
        frame=(r['scan_id'],r['frame_num'])
        if frame not in selected or sep<selected[frame][0]:
            selected[frame]=(sep,mjd,float(r[key+'mpro']),float(r[key+'sigmpro']))
    return sorted(selected.values(),key=lambda p:p[1])


def summarize_visits(target,band,archive,allow_legacy_status=False):
    visits=[]
    for _,mjd,mag,err in single_exposures(target,band,archive,allow_legacy_status):
        if not visits or mjd-visits[-1][-1][0]>20:
            visits.append([])
        visits[-1].append((mjd,mag,err))
    for visit in visits:
        mags=[v[1] for v in visit]
        center=median(mags)
        scatter=1.4826*median(abs(m-center) for m in mags)
        # Median standard error with a floor set by individual quoted errors.
        stat=1.2533*max(scatter,median(v[2] for v in visit))/math.sqrt(len(visit))
        yield dict(label=target['label'],band='W'+str(band),mjd_start=round(visit[0][0],5),mjd_end=round(visit[-1][0],5),n=len(visit),median_mag=round(center,4),scaled_mad_mag=round(scatter,4),stat_err_mag=round(stat,4))


def main(archive,output,allow_legacy_status=False):
    archive=Path(archive);output=Path(output)
    if not archive.is_dir():raise ValueError('archive directory missing')
    if output.exists():raise FileExistsError('derived output directory must be new')
    selection=targets()
    unknown=[]
    for t in selection:
        for kind in TABLES:
            read(t,kind,archive,allow_legacy_status)
            key=t['label']+'_'+kind
            rec=json.loads((archive/(key+'.json')).read_text(encoding='utf-8'))
            if 'status' not in rec:unknown.append(key)
    output.mkdir(parents=True)
    cats=[]; visits=[]; tests=[]
    for t in selection:
        cats.append(summarize_catalogs(t,archive,allow_legacy_status))
        for b in (1,2):
            v=list(summarize_visits(t,b,archive,allow_legacy_status))
            visits.extend(v)
            usable=[r for r in v if r['n']>=3]
            diffs=[]
            for a,c in itertools.combinations(usable,2):
                diff=c['median_mag']-a['median_mag']
                sigma=math.hypot(a['stat_err_mag'],c['stat_err_mag'],(0.02,0.06)[b-1])
                diffs.append((abs(diff)/sigma,diff,sigma))
            score,d,s=max(diffs) if diffs else ('','','')
            tests.append(dict(label=t['label'],band=f'W{b}',raw_neighbors=len(read(t,'neo',archive,allow_legacy_status)),n_kept=sum(r['n'] for r in v),visits_ge3=len(usable),signed_later_minus_earlier_mag=round(d,4) if diffs else '',difference_error_incl_floor=round(s,4) if diffs else '',max_pair_score=round(score,2) if diffs else ''))
    write(output,'catalog_photometry.csv',cats)
    write(output,'visits.csv',visits,['label','band','mjd_start','mjd_end','n','median_mag','scaled_mad_mag','stat_err_mag'])
    write(output,'variability.csv',tests)
    recs={p.stem:json.loads(p.read_text()) for p in sorted(archive.glob('*.json'))}
    (output/'replay-provenance.json').write_text(json.dumps(dict(allow_legacy_status=allow_legacy_status,unknown_status_receipts=unknown,note='Missing HTTP status remains unknown; a zero-row legacy response does not establish successful empty query or completeness.'),indent=2)+'\n',encoding='utf-8')
    (output/'queries.json').write_text(json.dumps(recs,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--archive',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--allow-legacy-status',action='store_true',help='Permit receipts with absent HTTP status for offline replay only')
    args=ap.parse_args()
    main(args.archive,args.output,args.allow_legacy_status)
