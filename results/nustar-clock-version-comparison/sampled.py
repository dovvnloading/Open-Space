from astropy.io import fits
from astropy.time import Time
import numpy as np, pathlib, json, hashlib, os
p=pathlib.Path(os.environ["NUSTAR_CLOCK_INPUTS"])
out=pathlib.Path(os.environ["NUSTAR_CLOCK_OUTPUT"])
out.mkdir(parents=True, exist_ok=True)
expected={234:'fe872bf396be44c6e4ee20551e4fde89c9cff0809e16057b5a35bcc2b0015e8e',235:'24a440096719c2f8c992d6c03e3c3ac769bc7f7638826b1a422d64a1a9d75d19',236:'f30744159e4dcff43a525d15318bd4fb88e561c1fda9fcc508e0db752bb30f96'}
for version,digest in expected.items():
 path=p/f'v{version}.fits.gz'
 if hashlib.sha256(path.read_bytes()).hexdigest()!=digest: raise ValueError(f'input digest mismatch {path.name}')
start=Time('2026-10-10T02:00:01',scale='utc'); stop=Time('2026-10-10T07:50:00',scale='utc')
duration=(stop-start).to_value('sec'); seconds=np.r_[np.arange(0,duration,60),duration]
curves={}; evidence={}
for version in (234,235,236):
 path=p/f'v{version}.fits.gz'
 with fits.open(path,memmap=False) as f:
  h=f['NU_FINE_CLOCK'].header; rows=f['NU_FINE_CLOCK'].data
  reference=Time(h['MJDREFI'],h['MJDREFF'],format='mjd',scale=h['TIMESYS'].lower())
  ts=(start.tt-reference).to_value('sec')+seconds
  knots=np.asarray(rows['TIME'],dtype=float); offset=np.asarray(rows['CLOCK_OFF_CORR'],dtype=float); freq=np.asarray(rows['CLOCK_FREQ_CORR'],dtype=float)
  i=np.searchsorted(knots,ts,side='right')-1
  covered=(i>=0)&(i<len(knots)-1)
  values=np.full(len(ts),np.nan); ii=i[covered]; gap=knots[ii+1]-knots[ii]; x=(ts[covered]-knots[ii])/gap
  values[covered]=(2*x**3-3*x**2+1)*offset[ii]+(x**3-2*x**2+x)*gap*freq[ii]+(-2*x**3+3*x**2)*offset[ii+1]+(x**3-x**2)*gap*freq[ii+1]
  curves[version]=values
  evidence[str(version)]={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'header_tstart':h['TSTART'],'header_tstop':h['TSTOP'],'min_knot':float(knots.min()),'max_knot':float(knots.max()),'covered_samples':int(covered.sum()),'uncovered_samples':int((~covered).sum()),'max_bracket_seconds':float(gap.max()) if len(gap) else None,'monotone_knots':bool(np.all(np.diff(knots)>0))}
result={'duration_s':duration,'samples':len(seconds),'versions':evidence,'pairs':{}}
for a,b in ((234,235),(235,236),(234,236)):
 d=(curves[b]-curves[a])*1e6; d=d[np.isfinite(d)]
 result['pairs'][f'{b}-{a}']={'common_samples':len(d),'min_us':float(d.min()) if len(d) else None,'max_us':float(d.max()) if len(d) else None,'max_abs_us':float(np.max(abs(d))) if len(d) else None,'max_abs_change_from_first_us':float(np.max(abs(d-d[0]))) if len(d) else None}
(out/'sampled.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
