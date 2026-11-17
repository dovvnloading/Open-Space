"""Offline consistency check for deposited source and control stamps."""
import argparse,csv,json,subprocess,sys,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true',required=True);args=p.parse_args();root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as tmp:
 subprocess.run([sys.executable,str(root/'replay_stamps.py'),'--stamps',str(root/'pilot-v6/stamps.npz'),'--centers',str(root/'pilot-v5/selected-centers.csv'),'--output',tmp],check=True)
 got=json.loads((Path(tmp)/'summary.json').read_text());expected=json.loads((root/'pilot-v6/summary.json').read_text());assert got==expected,(got,expected)
 with (Path(tmp)/'objects.csv').open(newline='') as f: a=list(csv.DictReader(f))
 with (root/'pilot-v6/objects.csv').open(newline='') as f: b=list(csv.DictReader(f))
 assert a==b,'per-source replay drift'
 out=subprocess.check_output([sys.executable,str(root/'replay_controls.py'),'--stamps',str(root/'pilot-v6/stamps.npz')],text=True);print(out.strip());assert 'count 677' in out and '0.044814546038339294' in out and '0.03725048279592218' in out
 print('PASS: frozen source and control stamps reproduce the deposited measurements')
