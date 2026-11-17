"""Reproduce and compare the frozen development-only measurements offline."""
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent

def compare(a, b, where):
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys(): raise AssertionError(f'{where}: keys differ')
        for k in a: compare(a[k], b[k], f'{where}.{k}')
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): raise AssertionError(f'{where}: row count differs')
        for i, (x, y) in enumerate(zip(a, b)): compare(x, y, f'{where}[{i}]')
    elif isinstance(a, bool) or isinstance(b, bool):
        if a != b or type(a) != type(b): raise AssertionError(f'{where}: classification differs')
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not (math.isnan(a) and math.isnan(b)) and not math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-8):
            raise AssertionError(f'{where}: {a} != {b}')
    else:
        if a != b: raise AssertionError(f'{where}: {a!r} != {b!r}')

def csv_rows(path):
    with path.open(newline='', encoding='utf-8') as f:
        result=list(csv.DictReader(f))
    for row in result:
        for key, value in row.items():
            if value is None: raise ValueError(f'Malformed row in {path}')
            # IDs, flags and categorical values are exact; floating measurements use tolerance.
            if key.lower() in ('sdss_id', 'id', 'mjd', 'kind', 'pass', 'qualifies', 'valid', 'conditional_threshold', 'baseline_success', 'dev_success', 'degree', 'baseline_dof', 'dev_dof', 'oiii_pixels') or value in ('', 'True', 'False'):
                continue
            try:
                row[key] = float(value)
            except ValueError:
                pass
    return result

def run(script, *args):
    subprocess.run([sys.executable, str(HERE/script), *map(str,args)], check=True)

with tempfile.TemporaryDirectory(prefix='dr20-replay-') as tmp:
    out=Path(tmp); manifest=out/'derived-manifest.json'; qa=out/'development-qa.csv'
    run('rebuild_development.py',out)
    run('qa_dev.py',out,manifest,qa)
    run('wing_proxy.py',out,manifest,qa,out/'wing-proxy.csv')
    run('joint_dev.py',out,manifest,qa,out/'joint-fits.csv',out/'real-injections.csv')
    run('joint_profile.py',out,manifest,qa,out/'profile-fits.csv')
    run('residual_blocks.py',out,manifest,out/'residual-blocks.json')
    run('rm17_spectral_controls.py',out,manifest,out/'rm17-spectral-controls.json')
    run('ztf_check.py'); run('ztf_dr24_check.py')
    for name in ('development-qa.csv','wing-proxy.csv','joint-fits.csv','real-injections.csv','profile-fits.csv'):
        compare(csv_rows(HERE/name),csv_rows(out/name),name)
        print('verified',name)
    for name in ('residual-blocks.json','rm17-spectral-controls.json'):
        compare(json.loads((HERE/name).read_text()),json.loads((out/name).read_text()),name)
        print('verified',name)
print('Offline replay passed; temporary products removed.')
