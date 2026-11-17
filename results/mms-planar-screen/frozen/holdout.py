"""Separate authorization required before reading heldout spacecraft FGM."""
import argparse,json
from pathlib import Path
import cdflib,numpy as np
from manifest import verify,sha,RAW
from development import authorize,EXPECTED
from plane import compare

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--manifest-sha256',required=True);ap.add_argument('--authorization',required=True);ap.add_argument('--authorization-sha256',required=True);ap.add_argument('--prediction',required=True);ap.add_argument('--prediction-sha256',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 m=verify(a.manifest,a.manifest_sha256)
 if json.loads((Path(__file__).parent/'config.json').read_text())!=EXPECTED:raise ValueError('configuration differs from implemented gates')
 if sha(a.prediction)!=a.prediction_sha256:raise ValueError('prediction digest')
 pred=json.loads(Path(a.prediction).read_text());authorize(a.authorization,a.authorization_sha256,a.manifest_sha256,'holdout',a.prediction_sha256)
 if pred['manifest_sha256']!=a.manifest_sha256 or pred['status']!='qualified':raise ValueError('no qualifying prediction')
 if a.output.resolve().drive.upper()!='D:' or a.output.exists():raise ValueError('new D output required')
 result={}
 for i in (2,3,4):
  f=cdflib.CDF(str(RAW/next(n for n in m['inputs'] if n.startswith(f'mms{i}_fgm_'))));t=(f.varget('Epoch')-pred['origin_tt2000'])/1e9;b=f.varget(f'mms{i}_fgm_b_gse_brst_l2')[:,:3];good=(f.varget(f'mms{i}_fgm_flag_brst_l2')==0)&np.isfinite(b).all(axis=1)&(np.abs(b)<=20000).all(axis=1)
  result[str(i)]=compare(t,b,good,pred['selected'],pred['sensitivity']['lag_low_s'][i-2],pred['sensitivity']['lag_high_s'][i-2])
 statuses=[v['status'] for v in result.values()]
 disposition='inconclusive' if any(s.startswith('inconclusive') for s in statuses) else ('inconsistent_with_declared_plane' if 'no_match' in statuses else ('insensitive_zero_lag_alternative' if any(v['zero_lag_alternative'] for v in result.values()) else 'conditional_plane_match'))
 with a.output.open('x',encoding='utf-8') as h:json.dump({'prediction_sha256':a.prediction_sha256,'results':result,'disposition':disposition,'scope':'conditional engineering screen; not physical significance'},h,indent=2,allow_nan=False);h.write('\n')
if __name__=='__main__':main()
