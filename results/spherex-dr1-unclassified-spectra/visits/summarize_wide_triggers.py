"""Audit each release and sign for all 52 original six-detector peak triggers."""
import csv
import json
from pathlib import Path
from collections import Counter

HERE=Path(__file__).resolve().parent
OUT=HERE/'wide_remaining'

def rows(file):
    with file.open() as f:return list(csv.DictReader(f))

def check(partial=False):
    peaks=[x for x in rows(HERE.parent/'wide/full/ranks.csv') if float(x['score'])>=5]
    local={x['source_id']:float(x['fraction']) for x in rows(HERE.parent/'wide/coadd_checks.csv')
           if x['mask']=='strict' and x['radius_arcsec']=='19' and x['sky_inner']=='36' and x['fraction']}
    select={}
    for sub in ('wide','wide_remaining'):
        folder=HERE/sub
        for t in rows(folder/'targets.csv'):
            k=(t['source_id'],t['channel'])
            if k in select:raise ValueError('overlap between follow-up lists')
            select[k]=folder/(t['source_id']+('_'+t['label'] if t['label'] else ''))
    out=[];not_done=0
    for r in peaks:
        key=(r['source_id'],r['peak_channel'])
        base=select[key]
        summary=Path(str(base)+'.json')
        if not summary.exists():
            if partial:not_done+=1;continue
            raise FileNotFoundError(summary)
        s=json.loads(summary.read_text())
        report=Path(str(base)+'.csv')
        exposures=rows(report) if s['attempted'] else []
        if len(exposures)!=s['attempted']:raise ValueError('missing cutout receipt')
        sign=1 if float(r['peak_signed'])>0 else -1
        for rel in ('dr1','qr3'):
            g=s['groups'][rel+'_all']
            a=[x for x in exposures if x['release']==rel]
            counts=Counter(x['status'].split(':')[0] for x in a)
            n=[g['bands'].get(k,{}).get('dates',0) for k in ('blue','core','red')]
            z=g.get('signed_score');f=g.get('fractional_contrast')
            out.append(dict(source_id=key[0],template_rank=r['template_rank'],channel=key[1],
                    center_um=r['peak_wave_um'],coadd_signed_score=r['peak_signed'],coadd_local_fraction=local.get(key[0],''),release=rel,
                    records_directory=str(base.parent.name),attempted=len(a),usable=counts['ok'],
                    aperture_flags=counts['aperture_flags'],few_blanks=counts['few_blanks'],
                    off_bin=counts['off_bin'],edge=counts['edge'],read_errors=counts['read_error'],
                    blue_dates=n[0],core_dates=n[1],red_dates=n[2],three_bins=int(z is not None),
                    all_bins_2plus=int(min(n)>=2),exposure_fraction=f if f is not None else '',
                    exposure_signed_score=z if z is not None else '',
                    same_sign_ge3=int(z*sign>=3) if z is not None else '',
                    local_sign_ge3=int(z*(1 if local[key[0]]>0 else -1)>=3) if z is not None and key[0] in local else '',
                    early_score=s['groups'][rel+'_early'].get('signed_score',''),
                    late_score=s['groups'][rel+'_late'].get('signed_score','')))
    stats={}
    for rel in ('dr1','qr3'):
        rr=[x for x in out if x['release']==rel]
        stats[rel]={k:sum(int(x[k] or 0) for x in rr) for k in ['three_bins','all_bins_2plus','attempted','usable',
                    'aperture_flags','few_blanks','off_bin','edge','read_errors','same_sign_ge3','local_sign_ge3']}
        stats[rel]['incomplete_source_windows']=len(rr)-stats[rel]['three_bins']
    pairs={r['source_id']:{x['release']:x for x in out if x['source_id']==r['source_id']} for r in out}
    stats['both_complete']=sum(all(x['three_bins'] for x in g.values()) for g in pairs.values())
    stats['both_same_sign_ge3']=[n for n,g in pairs.items() if all(x['same_sign_ge3']==1 for x in g.values())]
    stats['both_local_sign_ge3']=[n for n,g in pairs.items() if all(x['local_sign_ge3']==1 for x in g.values())]
    stats['same_sign_ge3_by_release']={rel:[x['source_id'] for x in out if x['release']==rel and x['same_sign_ge3']==1] for rel in ('dr1','qr3')}
    stats['not_done']=not_done
    if not partial:
        with (OUT/'all_raw_trigger_checks.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=out[0],lineterminator='\n');w.writeheader();w.writerows(out)
        (OUT/'all_raw_trigger_summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    return stats

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');args=p.parse_args()
    print(json.dumps(check(args.partial),indent=2))
