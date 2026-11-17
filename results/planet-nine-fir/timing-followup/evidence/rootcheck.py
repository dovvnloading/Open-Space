"""Root's independent receipt, coverage and arithmetic audit; no pixel search."""
import collections
import csv
import gzip
import hashlib
import json
import math
import subprocess
from pathlib import Path

BASE = Path('D:/Open-Space/campaigns')
REPO = Path('D:/Open-Space/worktrees/pr143-mapping-20261117')
ROOT = REPO / 'results/planet-nine-fir'
CODE = 'f6b97d8d2154f5ba9912cea2cb6abefb51c6e473'
CONFIG = {
    'dec1': ('296a210', '296a2102743f2c53229e11e9a819c8b6493444b3ac21d4c5d0440d7b70ae9e6a', '274576c27c69e9bce3a9600c50ade37bbcf688c904d33f9a60d80c1d1b681ec0'),
    'dec2': ('b000d53', 'b000d53af6f035b8012000fbfe05f32cbc64f36c2203b01956bef8d81156402f', 'cb74c32202f613792b675e246292e22641d710e6ca6e62b4b8a5c5315c7d192c'),
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads(path.read_bytes())

def table(path):
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    indexed = {(r['experiment'], int(r['draw'])): r for r in rows}
    assert len(indexed) == len(rows)
    return indexed

def wilson(k, n):
    z = 1.959963984540054
    p = k/n
    d = 1+z*z/n
    center = (p+z*z/(2*n))/d
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [center-half, center+half]

old = table(ROOT / 'archives/planet_draws.csv.gz')
historical = table(BASE / 'pr143-des-historical-581c9c4-ea2d60a-20261117/run/planet_draws.csv.gz')
report = {'kind': 'independent arithmetic and provenance on the same two native exposures', 'cases': {}}
for phase, (short, manifest_sha, auth_sha) in CONFIG.items():
    campaign = BASE / f'pr143-des-{phase}-f6b97d8-{short}-20261117'
    out = campaign / 'run'
    manifest_path = BASE / f'pr143-des-{phase}-packet-f6b97d8-20261117/manifest.json'
    assert sha(manifest_path) == manifest_sha
    manifest = read(manifest_path)
    assert manifest['phase'] == phase and manifest['git_head'] == CODE
    auth_path = campaign / 'root-authorization.json'
    assert sha(auth_path) == auth_sha
    auth = read(auth_path)
    assert auth['phase'] == phase and auth['manifest_sha256'] == manifest_sha
    assert Path(auth['out']).resolve() == out.resolve()
    assert sha(Path(auth['mapping'])) == auth['mapping_sha256'] == '001fc291a796cccfb547a65378d7fe1abc2baf7c875ca21d715efffbd6763404'
    receipt = read(campaign / 'execution-receipt.json')
    assert receipt['code_head'] == CODE
    args = receipt['arguments']
    assert args[args.index('--scenario')+1] == phase
    assert args[args.index('--manifest-sha256')+1] == manifest_sha
    assert args[args.index('--auth-sha256')+1] == auth_sha
    for option, expected in [('--manifest',manifest_path),('--auth',auth_path),('--out',out),('--mapping',Path(auth['mapping'])),('--cache',Path(auth['cache']))]:
        assert Path(args[args.index(option)+1]).resolve() == expected.resolve()
    assert sha(campaign / 'console.log') == receipt['console_sha256']
    for name, expected in receipt['outputs'].items():
        assert sha(out / name) == expected, name
    for name, expected in manifest['files'].items():
        blob = subprocess.check_output(['git','show',CODE+':results/planet-nine-fir/'+name], cwd=REPO)
        assert hashlib.sha256(blob).hexdigest() == expected, name
    run_manifest = read(out / 'run_manifest.json')
    assert run_manifest['scenario'] == phase and run_manifest['retrospective']
    for name, expected in run_manifest['input_sha256'].items():
        assert sha(Path(auth['cache']) / (name+'.response')) == expected
    for name, expected in run_manifest['code_sha256'].items():
        assert manifest['files'][name.replace('\\','/')] == expected
    assert read(out / 'run_status.json')['status'] == 'complete'
    rows = table(out / 'planet_draws.csv.gz')
    assert set(rows) == set(old) == set(historical) and len(rows) == 800
    outcomes = read(out / 'draw_outcomes.json')
    indexed = {(r['experiment'],int(r['draw'])):r for r in outcomes}
    assert set(indexed) == set(rows) and len(outcomes) == 800
    assert all(r['state_status']=='bound' and r['image_status']!='not_tested' and len(r['models'])==21 for r in outcomes)
    recovery = read(out / 'model_recovery.json')
    assert len(recovery) == 42
    summaries = []
    for exp in sorted({k[0] for k in rows}):
        o = [r for r in outcomes if r['experiment']==exp]
        models = [r for r in recovery if r['experiment']==exp]
        assert len(o) == 400 and len(models) == 21
        for j, model in enumerate(models):
            m = [r['models'][j] for r in o]
            for item in m:
                assert item['reason']!='not_tested'
                if item['postplant_eligible']:
                    assert isinstance(item['recovered'],bool)
                    assert item['reason']==('recovered' if item['recovered'] else 'not_recovered')
            k = sum(x['recovered'] is True for x in m)
            assert model['n_drawn']==400
            assert model['n_eligible_before_plant']==sum(r['preplant_eligible'] for r in o)
            assert model['n_eligible_after_plant']==sum(x['postplant_eligible'] for x in m)
            assert model['n_recovered']==k and model['fraction']==k/400
            assert max(abs(x-y) for x,y in zip(model['ci95_wilson'],wilson(k,400)))<1e-6
        flags = {}
        for col in [c for c in next(iter(rows.values())) if c.startswith('hit_')]:
            keys = [key for key in rows if key[0]==exp]
            flags[col] = {'recovered':sum(int(rows[key][col]) for key in keys),
                          'changed_vs_keyed_historical':[key[1] for key in keys if rows[key][col]!=historical[key][col]],
                          'changed_vs_original':[key[1] for key in keys if rows[key][col]!=old[key][col]]}
        summaries.append({'experiment':exp,'coverage':dict(collections.Counter(r['image_status'] for r in o)),'flags':flags})
    stats = read(out / 'difference_result.json')
    peaks = read(out / 'dense_peaks.json')
    assert stats['n_above']==0 and stats['dense_peaks']==len(peaks)==1
    assert stats['dense_max_snr']==peaks[0]['score']
    assert stats['exposure_ids']==[605285,605571]
    report['cases'][phase] = {'outputs_verified':len(receipt['outputs']),'source_blobs_verified':len(manifest['files']),
                             'draw_ids':800,'model_outcomes':16800,'conditional_model_checks':42,
                             'summaries':summaries,'dense_peak':peaks[0],
                             'grid_max_score':stats['max_snr'],'grid_threshold':stats['threshold'],
                             'grid_above_threshold':0,'exposure_ids':stats['exposure_ids'],
                             'manifest_sha256':manifest_sha,'authorization_sha256':auth_sha}
report['limits'] = 'No new independent observation or significance calibration; conditional injected-model fractions are not planet-existence or population probabilities. Both scenarios reuse the same two exposures. Earlier three historical dark-stress flips remain unexplained individually.'
destination = BASE / 'pr143-des-december-rootcheck-20261117.json'
with destination.open('x',encoding='utf-8',newline='\n') as f:
    json.dump(report,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','cases':{k:{x:v[x] for x in ['outputs_verified','source_blobs_verified','draw_ids','model_outcomes','conditional_model_checks','grid_above_threshold']} for k,v in report['cases'].items()},'report':str(destination)}))
