"""Guarded two-image descriptive mode. Never reads the prospective third image."""
import argparse, hashlib, json, math, os, sys
from pathlib import Path

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def authorize(manifest, authorization, expected_manifest, expected_authorization, input_paths, output):
    if digest(authorization) != expected_authorization: raise ValueError("authorization digest mismatch")
    if digest(manifest) != expected_manifest: raise ValueError('manifest digest mismatch')
    m=json.loads(Path(manifest).read_text()); a=json.loads(Path(authorization).read_text())
    if a.get('phase')!='two-image-development' or a.get('manifest_sha256')!=expected_manifest: raise ValueError('phase unauthorized')
    if len(input_paths)!=2 or len(set(map(lambda x:str(Path(x).resolve()),input_paths)))!=2: raise ValueError('exactly two distinct inputs required')
    if Path(output).exists(): raise ValueError('output already exists')
    if m['python']!=sys.version.split()[0]: raise ValueError('python version mismatch')
    if digest(__file__)!=m['code_sha256']: raise ValueError('code mutation')
    expected=('PUNCH_L2_CTM_20261101000029_v0l.fits','PUNCH_L2_CTM_20261101000829_v0l.fits')
    if tuple(r['name'] for r in m['inputs']) != expected: raise ValueError('unexpected fixed input names')
    for path, row in zip(input_paths,m['inputs']):
        if Path(path).name!=row['name'] or digest(path)!=row['sha256']: raise ValueError('input mismatch')
    if len(m['inputs'])!=2: raise ValueError('invalid manifest input count')
    return m

def characterize(paths):
    import numpy as np
    from astropy.io import fits
    frames=[]; uncertainties=[]; header=[]
    for path in paths:
        with fits.open(path,memmap=False) as h:
            x=h['PRIMARY DATA ARRAY']; u=h['UNCERTAINTY ARRAY']
            if x.data is None or u.data is None or x.data.ndim!=2 or x.data.shape!=u.data.shape: raise ValueError('invalid image/uncertainty schema')
            if not x.header.get('BUNIT') or x.header.get('BUNIT')!=u.header.get('BUNIT'): raise ValueError('uncertainty units unresolved')
            frames.append(np.array(x.data,dtype=float)); uncertainties.append(np.array(u.data,dtype=float)); header.append(x.header.copy())
    keys=('CTYPE1','CTYPE2','CUNIT1','CUNIT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','PC1_1','PC1_2','PC2_1','PC2_2','BUNIT')
    if frames[0].shape!=frames[1].shape or any(k not in header[0] or header[0][k]!=header[1].get(k) for k in keys): raise ValueError('matching WCS/units not established')
    valid=np.isfinite(frames[0]) & np.isfinite(frames[1]) & np.isfinite(uncertainties[0]) & np.isfinite(uncertainties[1]) & (uncertainties[0]>=0) & (uncertainties[1]>=0)
    if not np.any(valid): raise ValueError('no valid common pixels')
    d=frames[1][valid]-frames[0][valid]
    return {'mode':'descriptive','units':str(header[0]['BUNIT']),'shape':list(frames[0].shape),'common_finite_pixels':int(valid.sum()),'total_pixels':int(valid.size),'difference_median':float(np.median(d)),'difference_p05_p95':[float(v) for v in np.percentile(d,[5,95])], 'limits':'No physical significance; shared backgrounds, masks, calibration and registration unqualified.'}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--manifest',required=True); p.add_argument('--manifest-sha256',required=True); p.add_argument('--authorization',required=True); p.add_argument('--authorization-sha256',required=True); p.add_argument('--input',nargs=2,required=True); p.add_argument('--output',required=True); args=p.parse_args()
    m=authorize(args.manifest,args.authorization,args.manifest_sha256,args.authorization_sha256,args.input,args.output)
    import numpy, astropy
    if numpy.__version__!=m['numpy'] or astropy.__version__!=m['astropy']: raise ValueError('dependency version mismatch')
    result=characterize(args.input)
    with open(args.output,'x',encoding='utf-8',newline='\n') as f: json.dump(result,f,indent=2,allow_nan=False); f.write('\n')
if __name__=='__main__': main()
