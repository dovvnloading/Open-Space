"""Offline cross-checks for the frozen inventory, model outputs and positive control."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent


def table(p):
    if p.suffix=='.gz':
        with gzip.open(p,'rt') as f:return list(csv.DictReader(f))
    with p.open() as f:return list(csv.DictReader(f))


def run(partial=False,models=False):
    result=table(HERE/'tables/grid-fits.csv.gz')
    by_id={r['event']:r for r in result}
    assert len(by_id)==len(result)
    c=collections.Counter(r['survey'] for r in result)
    assert c=={'OGLE':2429,'KMT':6789,'MOA':699,'PRIME':970,'Gaia':89},c
    known=table(HERE/'other/known_matches.csv')
    assert len({r['known_id'] for r in known if r['status']=='known by published ID'})==5
    fs=[r for r in result if r['fit_status']=='ok']
    for r in fs:
        assert re.fullmatch('[0-9a-f]{64}',r['input_sha256']),r['event']
        assert 2450000<float(r['first_hjd'])<=float(r['last_hjd'])<2470000,r['event']
        assert int(r['n'])>=6 and float(r['ps_te'])>0 and float(r['fs_te'])>0 and float(r['fs_rho'])>0,r['event']
        assert all(float(r[k])>=0 for k in ['ps_gain','fs_gain','bin_gain']),r['event']
    assert all(r['survey']=='Gaia' or partial for r in result if r['fit_status']!='ok')
    bounds=table(HERE/'tables/kmt-reference-flux-grid.csv.gz')
    assert len(bounds)==len({r['event'] for r in bounds})
    for r in bounds:
        baseline=by_id[r['event']]
        assert r['status']=='ok' and baseline['fit_status']=='ok'
        assert int(r['n'])==int(baseline['n']),r['event']
        # The capped grid has the same geometric templates and a smaller
        # admissible source flux interval. Its best gain cannot increase.
        for k in ('ps_gain','fs_gain','bin_gain'):
            if float(r[k])>float(baseline[k])+max(.05,abs(float(baseline[k]))*1e-8):
                raise AssertionError((r['event'],k,r[k],baseline[k]))
    if not partial:assert len(bounds)==6789
    cad=table(HERE/'tables/cadence.csv.gz')
    assert len(cad)==len({r['event'] for r in cad})
    for r in cad:
        if r.get('error'):continue
        d=[float(r[k]) for k in ('days_0.03','days_0.1','days_0.25','days_0.5','days_1.5')]
        assert 0<=min(d) and d==sorted(d) and max(d)<=float(r['season_d'])+.00002,r['event']
    data=gzip.decompress((HERE/'controls/earth-figure1.tsv.gz').read_bytes())
    ctrl=json.loads((HERE/'controls/earth-recovery.json').read_text())
    assert hashlib.sha256(data).hexdigest()==ctrl['response_sha256']
    assert .02<ctrl['finite_no_blend']['te_d']<.04 and ctrl['finite_no_blend']['rho']>2.
    if models:
        import earth_control
        replay=earth_control.run(data)
        assert abs(replay['finite_no_blend']['te_d']-ctrl['finite_no_blend']['te_d'])<.00015
        assert abs(replay['finite_no_blend']['chi2']-ctrl['finite_no_blend']['chi2'])<.1
    aliases=table(HERE/'tables/survey-id-groups.csv.gz')
    seen=[n for row in aliases for n in row['alert_ids'].split('|')]
    assert len(seen)==len(set(seen)) and set(seen)==set(by_id)
    return dict(alert_ids=len(by_id),successful_fits=len(fs),kmt_bounded=len(bounds),
                cadence_rows=len(cad),literature_objects_matched=5,control_rows=ctrl['input_rows'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--partial',action='store_true')
    p.add_argument('--models',action='store_true');a=p.parse_args();print(json.dumps(run(a.partial,a.models)))
