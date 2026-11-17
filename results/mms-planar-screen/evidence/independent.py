"""Independent native-CDF audit of the frozen MMS1 null; no primary reader imports."""
import hashlib
import json
from pathlib import Path

import cdflib
import numpy as np

BASE = Path('D:/Open-Space/campaigns')
RAW = BASE / 'heliophysics-schema-20261117/raw'
MANIFEST = BASE / 'mms-planar-review-6ebb8aa4/manifest.json'
PRIMARY = BASE / 'mms-planar-dev-6ebb8aa4-f27e423-20261117/development.json'
OUTPUT = BASE / 'mms-independent-dev-6ebb8aa4-20261117.json'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_valid(cdf, name, unit, shape):
    attrs = cdf.varattsget(name)
    assert attrs['UNITS'] == unit
    assert list(cdf.varinq(name).Dim_Sizes) == shape
    values = cdf.varget(name)
    valid = (np.isfinite(values) & (values != attrs['FILLVAL']) &
             (values >= attrs['VALIDMIN']) & (values <= attrs['VALIDMAX']))
    return values, valid

manifest = json.loads(MANIFEST.read_text())
assert digest(MANIFEST) == 'f27e423a21d9f2bd38450fe3719a756848c7e2aa6acff71de0ea25315fe0c272'
names = ['mms1_fgm_brst_l2_20260930132803_v5.575.0.cdf',
         'mms1_fpi_brst_l2_dis-moms_20260930132803_v3.4.0.cdf']
for name in names:
    assert digest(RAW / name) == manifest['inputs'][name]
f, d = [cdflib.CDF(str(RAW / name)) for name in names]
ft, dt = f.varget('Epoch'), d.varget('Epoch')
assert np.all(np.diff(ft) > 0) and np.all(np.diff(dt) > 0)
low = int(cdflib.cdfepoch.compute_tt2000([2026, 9, 30, 13, 28, 3, 413, 441, 756]))
high = int(cdflib.cdfepoch.compute_tt2000([2026, 9, 30, 13, 30, 12, 157, 399, 823]))
fi, di = np.flatnonzero((ft >= low) & (ft <= high)), np.flatnonzero((dt >= low) & (dt <= high))
field, bvalid = read_valid(f, 'mms1_fgm_b_gse_brst_l2', 'nT', [4])
velocity, vvalid = read_valid(d, 'mms1_dis_bulkv_gse_brst', 'km/s', [3])
error, evalid = read_valid(d, 'mms1_dis_bulkv_err_brst', 'km/s', [3])
fg, dg = f.varget('mms1_fgm_flag_brst_l2'), d.varget('mms1_dis_errorflags_brst')
# Find nearest DIS record by a separate scalar search, preserving earlier tie choice.
nearest = np.asarray([int(di[np.argmin(np.abs(dt[di] - stamp))]) for stamp in ft[fi]])
aligned = np.abs(dt[nearest] - ft[fi]) <= 80_000_000
good = ((fg[fi] == 0) & bvalid[fi, :3].all(axis=1) &
        (dg[nearest] == 0) & vvalid[nearest].all(axis=1) &
        evalid[nearest].all(axis=1) & aligned)
coverage = dict(native_fgm_rows=len(ft), native_dis_rows=len(dt),
                cropped_fgm_rows=len(fi), cropped_first_fgm_native_index=int(fi[0]),
                cropped_dis_rows=len(di), good_aligned_fgm_rows=int(good.sum()),
                nearest_dis_within_tolerance=int(aligned.sum()))
t = ft[fi] - ft[fi[0]]
b, v = field[fi, :3], velocity[nearest]
counts = dict(tested_1s_windows=0, quality_cadence_windows=0,
              field_endpoint_windows=0, rotation_windows=0,
              mva_failures=0, speed_failures=0, scanned_to_first_selection=False)
maximum_eligible_rotation = None
for i, stamp in enumerate(t):
    j = int(np.searchsorted(t, stamp + 1_000_000_000))
    if j == len(t) or t[j] - stamp > 1_020_000_000:
        continue
    counts['tested_1s_windows'] += 1
    begin = int(np.searchsorted(t, stamp - 5_000_000_000, side='left'))
    end = int(np.searchsorted(t, stamp + 5_000_000_000, side='right'))
    if (stamp < 5_000_000_000 or t[-1] - stamp < 5_000_000_000 or end - begin < 1000
            or not good[begin:end].all() or np.max(np.diff(t[begin:end])) > 20_000_000
            or not np.isfinite(b[begin:end]).all() or not np.isfinite(v[begin:end]).all()):
        continue
    counts['quality_cadence_windows'] += 1
    a, z = b[i], b[j]
    na, nz = np.linalg.norm(a), np.linalg.norm(z)
    if min(na, nz) < 3:
        continue
    counts['field_endpoint_windows'] += 1
    angle = float(np.degrees(np.arccos(np.clip(np.dot(a, z)/(na*nz), -1, 1))))
    maximum_eligible_rotation = angle if maximum_eligible_rotation is None else max(angle, maximum_eligible_rotation)
    if angle >= 45:
        counts['rotation_windows'] += 1
assert counts['rotation_windows'] == 0, 'A rotation would require completing the later model gates'
primary = json.loads(PRIMARY.read_text())
assert coverage == primary['coverage'], (coverage, primary['coverage'])
assert counts == primary['gates'], (counts, primary['gates'])
assert primary['status'] == 'no_qualified_candidate' and primary['selected'] is None
result = dict(source_hashes={n: manifest['inputs'][n] for n in names},
              manifest_sha256=digest(MANIFEST), primary_sha256=digest(PRIMARY),
              independent_code_sha256=digest(Path(__file__)), coverage=coverage, gates=counts,
              first_last_native_indices=[int(fi[0]), int(fi[-1])],
                first_last_utc=list(cdflib.cdfepoch.encode_tt2000(ft[[fi[0], fi[-1]]])),
              maximum_eligible_rotation_deg=maximum_eligible_rotation,
              status='independent_native_counts_match',
              scope='Same MMS1 observations, separate calculation; overlapping windows are correlated. No other spacecraft science arrays read.')
with OUTPUT.open('x', encoding='utf-8') as out:
    json.dump(result, out, indent=2, allow_nan=False)
    out.write('\n')
print(json.dumps(result, indent=2))
