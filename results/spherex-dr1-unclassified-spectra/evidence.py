"""Reconcile each original DR1 trigger with its independent native-exposure test.

Run from any directory after producing the three native target selections. Missing
follow-up records fail loudly; native contrast is left blank for incomplete bands.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def check_input_hashes(field,run):
    with (run/'analysis.json').open() as f:
        hashes=json.load(f)['input_sha256']
    paths=[field/'spectra.csv',ROOT/'templates/relative_fnu.csv',ROOT/'analyze.py']
    if field!=ROOT:
        paths.append(field/'catalogs/model_redshifts.csv')
    for p in paths:
        if digest(p)!=hashes[p.name]:
            raise ValueError(f'Input hash differs for {p}')
    for p in (field/'catalogs').glob('*_query.json'):
        receipt=json.loads(p.read_text())
        if digest(p.with_name(p.name.replace('_query.json','.csv')))!=receipt['sha256']:
            raise ValueError(f'Catalog query digest differs for {p}')

def visits(field_id):
    locations=['nep'] if field_id=='nep' else ['wide','wide_remaining']
    out={}
    for d in locations:
        folder=ROOT/'visits'/d
        with (folder/'targets.csv').open() as f:
            targets=list(csv.DictReader(f))
        for t in targets:
            name=t['source_id']
            key=(name,int(t['channel']))
            if key in out:raise ValueError(f'Duplicate exposure follow-up selection {key}')
            stem=name+('_'+t['label'] if t['label'] else '')
            j=json.loads((folder/(stem+'.json')).read_text())
            with (folder/(stem+'.csv')).open() as f:
                obs=list(csv.DictReader(f))
            if len(obs)!=j['attempted'] or j['source_id']!=name:raise ValueError(f'Lost receipts {key}')
            data=dict(receipt_directory='visits/'+d,receipt=stem+'.csv',selection=t['selection'])
            for release in ('dr1','qr3'):
                rows=[r for r in obs if r['release']==release]
                data[release+'_attempted']=len(rows)
                data[release+'_usable']=sum(r['status']=='ok' for r in rows)
                data[release+'_failed_read']=sum(r['status'].startswith(('read_error','archive_bytes_changed')) for r in rows)
                group=j['groups'][release+'_all']
                data[release+'_core_dates']=group['bands'].get('core',{}).get('dates',0)
                data[release+'_full_bracket']=len(group['bands'])==3
                data[release+'_contrast_fraction']=group.get('fractional_contrast')
                data[release+'_diagnostic_score']=group.get('signed_score')
            out[key]=data
    return out

def measure(field_id,field,run):
    check_input_hashes(field,run)
    ranks=pd.read_csv(run/'ranks.csv').set_index('source_id')
    ext=pd.read_csv(field/'spectra.csv')
    parent=pd.read_csv(field/'sources.csv')
    eligibility=pd.read_csv(run/'eligibility.csv')
    filt=pd.read_csv(run/'filter_scores.csv')
    inj=pd.read_csv(run/'injections.csv')
    summary=json.loads((run/'analysis.json').read_text())
    ids=set(parent.loc[parent.disposition=='extracted','source_id'])
    if set(ext.source_id)!=ids or ext.groupby('source_id').channel.nunique().min()!=102 or len(ext)!=102*len(ids):
        raise ValueError('A source-channel record is missing or duplicated')
    if set(ranks.index)!=set(eligibility.loc[eligibility.status=='eligible','source_id']):
        raise ValueError('Ranked and eligible source lists disagree')
    if len(inj)!=summary['n_injections']:raise ValueError('Missing injection trials')
    grouped=inj.groupby('amplitude_frac').recovered.agg(['sum','count'])
    pairs=visits(field_id)
    triggers=set((key,int(row.peak_channel)) for key,row in ranks.iterrows() if row.score>=5)
    if not triggers<=pairs.keys():raise ValueError('Unaccounted discovery-spectrum trigger')
    records=[]
    for (name,ch),data in pairs.items():
        rank=ranks.loc[name]
        original=(name,ch) in triggers
        if original:
            score=rank.peak_signed
        else:
            adjusted=pd.read_csv(run/'peer_channel_check.csv').set_index('source_id').loc[name]
            if int(adjusted.channel_adjusted)!=ch:raise ValueError(f'Unexpected follow-up center {name} {ch}')
            width=int(adjusted.fwhm_ch_adjusted)
            matches=filt[(filt.source_id==name)&(filt.channel==ch)&(filt.fwhm_ch==width)]
            score=float(matches.score.item())
        q=data['qr3_diagnostic_score']
        records.append(dict(field=field_id,source_id=name,channel=ch,original_trigger=original,
            coadd_local_score=score,qr3_same_sign_ge3=(None if q is None else bool(np.sign(score)*q>=3)),**data))
    rec=pd.DataFrame(records)
    raw=rec[rec.original_trigger]
    return records,dict(parent=len(parent),extracted=len(ids),ranked=len(ranks),raw_triggers=len(triggers),
        native_feature_locations=len(rec),native_sources=rec.source_id.nunique(),qr3_complete_raw=int(raw.qr3_full_bracket.sum()),
        dr1_complete_raw=int(raw.dr1_full_bracket.sum()),qr3_same_sign_ge3_raw=int(raw.qr3_same_sign_ge3.eq(True).sum()),
        native_checks_with_qr3_triplet=int(rec.qr3_full_bracket.sum()),
        target_exposure_tests=int(rec.qr3_attempted.sum()+rec.dr1_attempted.sum()),
        failed_reads=int(rec.qr3_failed_read.sum()+rec.dr1_failed_read.sum()),
        eligible_trials={str(a):{'recovered':int(r['sum']),'trials':int(r['count'])} for a,r in grouped.iterrows()})

def main():
    records=[];counts={}
    for name,field,rankdir in [('wide',ROOT/'wide','full'),('nep',ROOT,'nep_two_band')]:
        result,counts[name]=measure(name,field,field/rankdir)
        records.extend(result)
    with (ROOT/'followup_by_feature.csv').open('w', encoding='utf-8', newline='') as f:
        pd.DataFrame(records).to_csv(f, index=False, float_format='%.9g', lineterminator='\n')
    with (ROOT/'search_counts.json').open('w', encoding='utf-8', newline='\n') as f:
        f.write(json.dumps(counts, indent=2) + '\n')

if __name__=='__main__':
    main()
