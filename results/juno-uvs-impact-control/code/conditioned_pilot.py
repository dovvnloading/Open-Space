"""Metadata-selected interval, exploratory conditional-label diagnostic only."""
import pathlib,json,numpy as np
p=pathlib.Path(__file__).resolve().parents[1]/'clean-screen';x=np.load(p/'clean-photons.npz')['photons'];a=json.load(open(p/'clean-acquisition.json'));t=x['EPHEMERIS_TIME']-a['et_bounds'][0]
lat=x['JUPITER_LAT_400km'];lon=x['JUPITER_LON_400km'];inc=x['INCIDENCE_ANGLE'];em=x['EMISSION_ANGLE'];w=x['WAVELENGTH']
on=(abs(lat)<=90)&(lon>=0)&(lon<360)&(em>=0)&(em<90);night=on&(inc>90)&(inc<=180);long=(w>=170)&(w<200)
bin20=np.floor(t/.02).astype(int);cells=np.floor((np.clip(lat,-90,90)+90)/2).astype(int)*180+np.floor(np.clip(lon,0,359.999)/2).astype(int)
def scan(l):
 counts=np.bincount(bin20[night&l],minlength=1000);triggers=[]
 for j in np.flatnonzero(counts>=20):
  q=night&l&(bin20==j);v,c=np.unique(cells[q],return_counts=True);frac=c.max()/counts[j]
  if frac>=.6:triggers.append({'seconds':float(j*.02),'long_counts':int(counts[j]),'fraction':float(frac),'cell':int(v[c.argmax()])})
 return {'triggers':triggers,'count_cut_bins':int((counts>=20).sum()),'max_long_counts':int(counts.max())}
# Conditional labels preserve 100-ms time strata, 16-pixel slit segments and night status.
# This preserves empirical local count/background structure, not a physical calibrated null.
groups=np.floor(t/.1).astype(int)*10000+(x['DETECTOR_Y']//16)*10+night.astype(int);indices=[np.flatnonzero(groups==g) for g in np.unique(groups)]
null=[]
for seed in range(20):
 rng=np.random.default_rng(seed);z=long.copy()
 for ids in indices:z[ids]=rng.permutation(z[ids])
 s=scan(z);null.append({'seed':seed,'triggers':len(s['triggers']),'count_cut_bins':s['count_cut_bins'],'max_long_counts':s['max_long_counts']})
r={'selection':'metadata-first fixed 20-second interval','source':a['source'],'rows':len(x),'seconds_sampled':[float(t.min()),float(t.max())] if len(x) else [],'on_disc_photons':int(on.sum()),'nightside_photons':int(night.sum()),'nightside_longwave_photons':int((night&long).sum()),'nightside_incidence_percentiles':np.percentile(inc[night],[1,50,99]).tolist() if night.any() else [],'measured':scan(long),'conditional_label_diagnostics':null,'limits':['No calibrated false-alarm probability or completeness','Binary wavelength-label randomization breaks exact detector-X/wavelength relationship; diagnostic only','STIM active; no independent per-photon STIM exclusion yet','Overflow-clear sampled HK is necessary but not proof of loss-free exposure','Any physical candidate needs neighboring rotations and response-model validation']};(p/'conditioned-summary.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
