"""Aggregate actual native cutout coverage, flags and signed line-to-flank contrasts."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
records=[]
for field in ['nep','wide']:
    directory=ROOT/field
    with (directory/'targets.csv').open() as f:
        targets=list(csv.DictReader(f))
    for target in targets:
        name=target['source_id']
        suffix=target['label']
        stem=name+('_'+suffix if suffix else '')
        summary=json.loads((directory/f'{stem}.json').read_text())
        with (directory/f'{stem}.csv').open() as f:
            exposures=list(csv.DictReader(f))
        for release in ['dr1','qr3']:
            group=summary['groups'][release+'_all']
            sample=[r for r in exposures if r['release']==release]
            record=dict(field=field,source_id=name,channel=target['channel'],center_um=target['wave_um'],selection=target['selection'],release=release,
                        attempted=len(sample),usable=sum(r['status']=='ok' for r in sample),
                        flagged_apertures=sum(r['status']=='aperture_flags' for r in sample),
                        few_blanks=sum(r['status']=='few_blanks' for r in sample),
                        failures=sum(r['status'].startswith('read_error') for r in sample))
            for label,mask in [('overflow',2),('bloom',1<<26),('persist',1<<17),
                               ('ghost',(1<<22)|(1<<23)|(1<<24)),('halo',(1<<28)|(1<<29))]:
                record[label+'_apertures']=sum(bool(int(r.get('flag_or') or 0)&mask) for r in sample)
            for bin_name in ['blue','core','red']:
                b=group['bands'].get(bin_name,{})
                record[bin_name+'_dates']=b.get('dates',0)
                record[bin_name+'_flux']=b.get('median','')
            record['contrast_fraction']=group.get('fractional_contrast','')
            record['diagnostic_score']=group.get('signed_score','')
            for epoch in ['early','late']:
                record[epoch+'_score']=summary['groups'][release+'_'+epoch].get('signed_score','')
            records.append(record)
with (ROOT/'summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
