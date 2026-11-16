"""Targeted known-event recovery diagnostic, not a significance estimator."""
import json,pathlib
import numpy as np
P=pathlib.Path(__file__).resolve().parents[1]
x=np.load(P/'control-photons.npz')['photons'];meta=json.load(open(P/'acquisition.json')); t=x['et']-meta['target_et']
valid=np.isfinite(x['wave'])&(x['wave']>0)&np.isfinite(x['lat'])&np.isfinite(x['lon'])
near=valid&(abs(x['lat']-53)<3)&(abs((x['lon']-200+180)%360-180)<3)
# Counts, not calibrated irradiance. Geometry in this RDR is at 400 km.
edges=np.arange(-2,2.0001,.001); hist,_=np.histogram(t[near&(x['wave']>=170)&(x['wave']<200)],edges)
k=int(np.argmax(hist));peak=(edges[k]+edges[k+1])/2
on=near&(abs(t-peak)<.04);off=near&(abs(t-peak)>.2)
def stats(mask):
 return {'n':int(mask.sum()),'short_130_160':int((mask&(x['wave']>=130)&(x['wave']<160)).sum()),'long_170_200':int((mask&(x['wave']>=170)&(x['wave']<200)).sum())}
result={'status':'targeted published-control diagnostic only','photon_count':len(x),'targeted_region_count':int(near.sum()),'peak_1ms_longwave_counts':int(hist[k]),'peak_time_minus_published_second_s':float(peak),'on_80ms':stats(on),'off_excluding_400ms':stats(off),'coordinate_medians_at_on':{key:float(np.median(x[key][on])) if on.any() else None for key in ['lat','lon','angle','spin','digital','analog']},'unresolved':['Masks and housekeeping not yet applied','No null-trials or injection recovery','On/off exposure not normalized','400 km RDR geometry differs from fitted source altitude','Not a new event or independent confirmation']}
json.dump(result,open(P/'control-summary.json','w'),indent=2);print(json.dumps(result,indent=2))
np.savetxt(P/'control-lightcurve.csv',np.c_[(edges[1:]+edges[:-1])/2,hist],delimiter=',',header='time_relative_to_published_second_s,longwave_counts_1ms',comments='')
