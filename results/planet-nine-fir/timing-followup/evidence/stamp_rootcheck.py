"""Independent checks of the fixed, previously authorized DES stamps only."""
import hashlib
import json
from pathlib import Path
import numpy as np

base = Path('D:/Open-Space/campaigns/pr143-des-fixed-feature-ca3d18bf-20261117')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = json.loads((base/'execution-receipt.json').read_bytes())
assert receipt['head']=='ca3d18bff073623ee3eba1ce26a3b3297b02aff2' and receipt['status']=='complete'
assert sha(base/'console.log')==receipt['console_sha256']
for name, expected in receipt['output_sha256'].items():
    assert sha(base/'run'/name)==expected
assert sha(base/'authorization.json')=='80f77d91347a81a49088437a162e909ef391bca07e568bbd9a118bc5e783ee1d'
data = np.load(base/'run/native_stamps.npz', allow_pickle=False)
record = json.loads((base/'run/diagnostic.json').read_bytes())
yy,xx = np.mgrid[-12:13,-12:13]
distance_squared = xx**2+yy**2
aperture = distance_squared<=16
annulus = (distance_squared>=81)&(distance_squared<=144)
assert aperture.sum()==49 and annulus.sum()==192
checked = 0
differences = []
case_results = {}
for case, result in record.items():
    assert result['reconstruction']['status']=='passed'
    lobes = result['fixed_lobes']
    score = (lobes['first']-lobes['second']-lobes['center_adu'])/lobes['noise_adu']
    assert abs(score-lobes['score'])<1e-12
    assert abs(score-result['original_score'])<1e-5
    # A zero-motion score uses the identical difference-map pixel twice.
    fixed = result['hypothesis_lobes']['sky_fixed']
    assert fixed['first']==fixed['second']
    for control in result['fixed_controls']:
        if 'score' in control:
            a,b = control['lobes']
            assert abs((a-b-lobes['center_adu'])/lobes['noise_adu']-control['score'])<1e-12
    measurements = {}
    for name, summary in result['stamps'].items():
        if name not in data:
            assert summary['status']=='outside_image'
            continue
        stamp = data[name].astype(np.float64)
        dq = data[name+'_dq']
        good = (dq==0)&np.isfinite(stamp)
        valid_aperture = aperture&good
        valid_annulus = annulus&good
        assert valid_aperture.sum()==summary['valid_aperture_pixels']
        assert valid_annulus.sum()==summary['valid_sky_pixels']
        if summary['status']!='ok':
            continue
        sky = float(np.median(stamp[valid_annulus]))
        flux = float(stamp[valid_aperture].sum()-valid_aperture.sum()*sky)
        delta = abs(flux-summary['flux_aperture'])
        assert delta<0.002, (name,delta)
        differences.append(delta)
        design = np.column_stack((np.ones(valid_annulus.sum()),xx[valid_annulus],yy[valid_annulus]))
        coef = np.linalg.solve(design.T@design,design.T@stamp[valid_annulus])
        model = coef[0]+coef[1]*xx+coef[2]*yy
        plane_flux = float((stamp-model)[valid_aperture].sum())
        assert abs(plane_flux-summary['background_crop_sensitivity']['radius4_plane_flux'])<1e-8
        measurements[name]={'median_background_flux_float64_adu':flux,'plane_background_flux_adu':plane_flux}
        checked += 1
    case_results[case]={'shifted_score':score,'measurements':measurements,
                        'same_pixel_score_is_algebraic_zero_motion':True}
answer={'status':'passed','native_stamp_measurements':checked,
        'max_float32_vs_float64_aperture_delta_adu':max(differences),
        'cases':case_results,
        'limits':'Same two native exposures and selected positions; no independent confirmation, PSF likelihood, empirical false-alarm probability or moving-source detection. The zero-motion lobe score is algebraically trivial. Detector-fixed comparison is descriptive and correlated.'}
out=Path('D:/Open-Space/campaigns/pr143-fixed-stamp-rootcheck-20261117.json')
with out.open('x',encoding='utf-8',newline='\n') as f:
    json.dump(answer,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','native_stamps_checked':checked,'max_aperture_numeric_delta':max(differences),'report':str(out)}))
