"""Independent same-data count/bin/power audit; no physical interpretation."""
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

BASE=Path(__file__).resolve().parent
RAW=Path('D:/Open-Space/campaigns/gbm-tte-metadata-20261117/raw')
MANIFEST=Path('D:/Open-Space/campaigns/gbm-exploratory-review-e8cc5a76/manifest.json')
PRIMARY=BASE/'characterization.json'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(MANIFEST)=='6240381db20cbb3e46fe7fc2b7530267aace0f343ad7f1f7d8c45e8482fab33f'
assert sha(PRIMARY)=='affc36bc8bc27bad9e44152291023830a9cacd7f35cd4778769bde158782ddad'
m=json.loads(MANIFEST.read_bytes());saved=json.loads(PRIMARY.read_bytes());results=[];powers={}
specs=[('261002867','n5',812666888.812348,-.12803,2),
       ('261002867','n2',812666888.812348,-.12803,2),
       ('261007236','nb',813044394.710608,2.84802,5),
       ('261007236','n8',813044394.710608,2.84802,5)]
for bid,det,trig,start,n in specs:
    name=f'glg_tte_{det}_bn{bid}_v00.fit';path=RAW/name
    assert sha(path)==m['tte_sha256'][name]
    with fits.open(path,memmap=True) as h:
        obj=h[0].header['OBJECT']
        t=np.asarray(h['EVENTS'].data['TIME'])-trig
        channel=np.asarray(h['EVENTS'].data['PHA'])
        eb=h['EBOUNDS'].data
        selected=np.asarray(eb['CHANNEL'])[(eb['E_MIN']>=50)&(eb['E_MAX']<=300)]
        use=(t>=start)&(t<start+4*n)
        index=np.floor((t[use]-start)*64).astype(int)
        c=channel[use]
        full=np.bincount(index,minlength=n*256)
        overflow=np.bincount(index[c==127],minlength=n*256)
        band=np.bincount(index[np.isin(c,selected)],minlength=n*256)
        histogram=np.bincount(c.astype(int),minlength=128)
    got=saved['bursts'][bid][det]
    for values,key in [(full,'fullband_lightcurve'),(band,'band_lightcurve'),(histogram,'pha_histogram')]:
        assert np.array_equal(values,got[key]),(bid,det,key)
    assert int(full.sum())==got['fullband_count'] and int(overflow.sum())==got['overflow_count']
    p=[]
    for w in range(n):
        sl=slice(w*256,(w+1)*256);q=band[sl]
        value=2*np.abs(np.fft.rfft(q-q.mean()))**2/q.sum() if q.sum() else np.zeros(129)
        value[-1]/=2;p.append(value[20:129])
        record=got['windows'][w]
        assert int(full[sl].sum())==record['fullband_count']
        assert int(overflow[sl].sum())==record['overflow_count']
        dead=float((full[sl].sum()-overflow[sl].sum())*2.6e-6+overflow[sl].sum()*1e-5)
        assert abs(dead-record['nominal_deadtime_seconds'])<1e-12
    p=np.asarray(p);assert np.allclose(p,got['fourier_power'],rtol=0,atol=1e-12)
    powers[bid,det]=p
    results.append(dict(burst=bid,detector=det,header_object=obj,full_count=int(full.sum()),
                        band_count=int(band.sum()),overflow_count=int(overflow.sum()),windows=n,
                        selected_channels=len(selected),max_power_numeric_discrepancy=float(np.max(np.abs(p-np.asarray(got['fourier_power']))))))
for bid,a,b in [('261002867','n5','n2'),('261007236','nb','n8')]:
    assert np.allclose(np.minimum(powers[bid,a],powers[bid,b]),saved['bursts'][bid]['joint_power'],rtol=0,atol=1e-12)
out=dict(status='all_counts_and_descriptive_arrays_reproduced',files=results,
         primary_sha256=sha(PRIMARY),code_sha256=sha(Path(__file__)),
         scope='Same four files, separate count/bin calculation. No p-value, ranking, detection or source interpretation.')
with (BASE/'independent_check.json').open('x',encoding='utf-8') as f:
    json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
