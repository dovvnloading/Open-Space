"""Offline consistency check for deposited source and control stamps."""
import argparse,csv,json,math,subprocess,sys,tempfile
from pathlib import Path

REL_TOL=1e-10
ABS_TOL=1e-12
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true',required=True);p.parse_args();root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as tmp:
 subprocess.run([sys.executable,str(root/'replay_stamps.py'),'--stamps',str(root/'pilot-v6/stamps.npz'),'--centers',str(root/'pilot-v5/selected-centers.csv'),'--output',tmp],check=True)
 got=json.loads((Path(tmp)/'summary.json').read_text());expected=json.loads((root/'pilot-v6/summary.json').read_text());assert got==expected,(got,expected)
 with (Path(tmp)/'objects.csv').open(newline='') as f:a=list(csv.DictReader(f))
 with (root/'pilot-v6/objects.csv').open(newline='') as f:b=list(csv.DictReader(f))
 assert len(a)==len(b),'source count drift'
 exact={'x','y','n_controls','passes'}
 for index,(actual,wanted) in enumerate(zip(a,b)):
  assert actual.keys()==wanted.keys(),f'columns differ in row {index}'
  for key,value in actual.items():
   reference=wanted[key]
   if key in exact or not value or not reference:
    assert value==reference,f'row {index} {key}: {value!r} != {reference!r}'
   else:
    assert math.isclose(float(value),float(reference),rel_tol=REL_TOL,abs_tol=ABS_TOL),f'row {index} {key}: {value} != {reference}'
 out=subprocess.check_output([sys.executable,str(root/'replay_controls.py'),'--stamps',str(root/'pilot-v6/stamps.npz')],text=True);print(out.strip())
 values=out.split();assert values[0]=='count' and int(values[1])==677 and values[2]=='difference_sigma_mad' and values[4]=='formal_median'
 assert math.isclose(float(values[3]),0.044814546038339294,rel_tol=REL_TOL,abs_tol=ABS_TOL)
 assert math.isclose(float(values[5]),0.03725048279592218,rel_tol=REL_TOL,abs_tol=ABS_TOL)
 print('PASS: frozen source and control stamps reproduce the deposited measurements')
