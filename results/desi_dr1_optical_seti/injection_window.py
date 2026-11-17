"""Save lossless native samples surrounding the fixed 5000–5100 Å rest window."""
import json
from pathlib import Path

import numpy as np

REST_MIN = 5000.
REST_MAX = 5100.


def pack_injection_window(records, meta, output_dir):
    """Supply the complete eligible B-arm list and target metadata from desi_shard.targets."""
    tid = int(meta['targetid'])
    out = Path(output_dir) / f"{meta['survey']}-{meta['program']}-{meta['healpix']}"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{tid}.npz"
    payload = {'metadata': np.array(json.dumps(meta, separators=(',',':')))}
    for i,e in enumerate(records):
        # Keep external samples needed by edge columns, local baseline and peak comparison.
        lo, hi = np.searchsorted(e.wave, np.array([REST_MIN,REST_MAX])*e.doppler)
        margin = 2*e.resolution.shape[0]+10
        sel = slice(max(0,lo-margin),min(len(e.wave),hi+margin+1))
        for name in ('wave','flux','template','mask','ivar'):
            payload[f'{i}_{name}'] = np.asarray(getattr(e,name))[sel]
        payload[f'{i}_resolution'] = np.asarray(e.resolution)[:, sel]
        payload[f'{i}_expid'] = np.asarray(e.expid)
        payload[f'{i}_velocity_kms'] = np.asarray(e.velocity_kms)
        payload[f'{i}_velocity_error_kms'] = np.asarray(e.velocity_error_kms)
    # Write a temporary file then rename atomically. ZIP compression is
    # lossless. Per-star files let discarded survey shards be safely regenerated.
    tmp = path.with_suffix('.partial')
    with tmp.open('wb') as f:
        np.savez_compressed(f, **payload)
    tmp.replace(path)
    return path


def read_window(path):
    from laser_search import Exposure
    with np.load(path,allow_pickle=False) as f:
        meta=json.loads(str(f['metadata']))
        records=[]
        for i in range(len(meta['exposures'])):
            values = {key: f[f'{i}_{key}'] for key in
                      ('wave','flux','template','mask','ivar','resolution')}
            values.update(expid=int(f[f'{i}_expid']),velocity_kms=float(f[f'{i}_velocity_kms']),
                          velocity_error_kms=float(f[f'{i}_velocity_error_kms']))
            records.append(Exposure(**values))
    return meta,records
