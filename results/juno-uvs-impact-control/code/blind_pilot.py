"""Frozen first-500k-row PJ74 engineering pilot. No discovery significance."""
import pathlib,json,hashlib
import numpy as np
from astropy.io import fits
from astropy.time import Time
P=pathlib.Path(__file__).resolve().parents[1]/'newer-pj74'
h=fits.Header.fromstring((P/'header-2.txt').read_text(),sep='');types={'J':'>i4','E':'>f4','D':'>f8'}
dtype=np.dtype([(h[f'TTYPE{i}'],types[h[f'TFORM{i}'].strip()]) for i in range(1,h['TFIELDS']+1)])
assert dtype.itemsize==h['NAXIS1']==84
x=np.fromfile(P/'pilot.bin',dtype=dtype);assert len(x)==500000
et=x['EPHEMERIS_TIME'];assert np.all(np.diff(et)>=0)
t0=float(et[0]); bins=np.floor((et-t0)/.020).astype(int)
w=x['WAVELENGTH'];lat=x['JUPITER_LAT_400km'];lon=x['JUPITER_LON_400km'];emi=x['EMISSION_ANGLE']
on=np.isfinite(lat)&(abs(lat)<=90)&np.isfinite(lon)&(lon>=0)&(lon<360)&(emi>=0)&(emi<90)
long=np.isfinite(w)&(w>=170)&(w<200)
quality=np.isfinite(x['WEIGHTED_COUNT'])&(x['WEIGHTED_COUNT']>0)&np.isfinite(x['DIGITAL_DEADTIME_FACTOR'])&(x['DIGITAL_DEADTIME_FACTOR']>=1)&np.isfinite(x['ANALOG_DEADTIME_FACTOR'])&(x['ANALOG_DEADTIME_FACTOR']>=1)
cells=np.floor((np.clip(lat,-90,90)+90)/2).astype(int)*180+np.floor(np.clip(lon,0,359.999)/2).astype(int)
def scan(b):
 m=on&long; labels,counts=np.unique(b[m],return_counts=True); out=[]
 for label,count in zip(labels[counts>=20],counts[counts>=20]):
  q=m&(b==label); cc,nc=np.unique(cells[q],return_counts=True);j=nc.argmax();frac=float(nc[j]/count)
  if frac>=.6:out.append({'bin':int(label),'seconds':float(label*.02),'long_counts':int(count),'dominant_cell_fraction':frac,'cell':int(cc[j]),'positive_weight_deadtime_fraction':float(quality[q].mean()),'unique_detector_x':int(len(np.unique(x['DETECTOR_X'][q]))),'unique_detector_y':int(len(np.unique(x['DETECTOR_Y'][q]))),'pulse_median':float(np.median(x['PULSE_HEIGHT'][q]))})
 return out,{'bins_over_count_cut':int((counts>=20).sum()),'max_long_count':int(counts.max()) if len(counts) else 0}
triggers,counts=scan(bins);null=[]
for seed in range(20):
 a,b=scan(np.random.default_rng(seed).permutation(bins));null.append({'seed':seed,'triggers':len(a),**b})
off_counts=np.unique(bins[~on&long],return_counts=True)[1]
summary={'status':'exploratory frozen pilot, not impact classification','source':json.load(open(P/'pilot-receipt.json')),'utc_start':Time(2451545+t0/86400,format='jd',scale='tdb').utc.isot,'seconds_span':float(et[-1]-et[0]),'rows':len(x),'on_disc_rows':int(on.sum()),'on_disc_longwave_rows':int((on&long).sum()),'positive_weight_deadtime_rows':int(quality.sum()),'off_disc_longwave_count_bins_ge20':int((off_counts>=20).sum()),'off_disc_max_longwave_count':int(off_counts.max()) if len(off_counts) else 0,'geometry_diagnostic':{'on_disc_longwave_dayside_fraction':float(np.mean(x['INCIDENCE_ANGLE'][on&long]<90)),'incidence_angle_p01_p50_p99':np.percentile(x['INCIDENCE_ANGLE'][on&long],[1,50,99]).tolist()},'counts':counts,'triggers':triggers,'shuffle_diagnostics':null,'limits':['Photon mask and housekeeping interpretation required before physical inference','Shuffles do not preserve scan geometry/exposure','No completeness or global significance measured','Fixed angular cells can miss cell-boundary events','No spectrum model or adjacent-rotation confirmation']}
(P/'pilot-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
