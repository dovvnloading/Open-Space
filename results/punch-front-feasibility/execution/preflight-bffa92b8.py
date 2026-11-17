"""Guarded, descriptive two-CTM development; no prospective image access."""
import argparse, hashlib, json, sys
from contextlib import ExitStack
from pathlib import Path

EXPECTED=('PUNCH_L2_CTM_20261101000029_v0l.fits','PUNCH_L2_CTM_20261101000829_v0l.fits')
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def resolved(p):
    path=Path(p)
    if not path.is_absolute() or path.drive.upper()!='D:': raise ValueError('absolute D: path required')
    return str(path.resolve())
def authorize(manifest,authorization,expected_manifest,expected_auth,inputs,output):
    if digest(manifest)!=expected_manifest or digest(authorization)!=expected_auth: raise ValueError('digest mismatch')
    m=json.loads(Path(manifest).read_text()); a=json.loads(Path(authorization).read_text())
    out=resolved(output)
    if a.get('phase')!='two-image-development' or a.get('manifest_sha256')!=expected_manifest or a.get('output_path')!=out or m.get('output_path')!=out: raise ValueError('phase/output unauthorized')
    if Path(out).exists(): raise ValueError('output collision')
    if len(inputs)!=2 or len(m['inputs'])!=2 or tuple(x['name'] for x in m['inputs'])!=EXPECTED: raise ValueError('input scope')
    if m['python']!=sys.version.split()[0] or digest(__file__)!=m['code_sha256']: raise ValueError('runtime/code mismatch')
    for path,row in zip(inputs,m['inputs']):
        if resolved(path)!=row['path'] or Path(path).name!=row['name'] or digest(path)!=row['sha256']: raise ValueError('input path/hash mismatch')
    if resolved(inputs[0])==resolved(inputs[1]): raise ValueError('duplicate input')
    return m

def inspect_headers(hdus):
    from astropy.wcs import WCS
    images=[]; signatures=[]
    for h in hdus:
        x=h['PRIMARY DATA ARRAY']; u=h['UNCERTAINTY ARRAY']; hx=x.header; hu=u.header
        shape=tuple(hx.get('NAXIS'+str(i)) for i in (1,2))
        if hx.get('NAXIS')!=2 or hu.get('NAXIS')!=2 or any(not isinstance(v,int) or v<=0 for v in shape) or shape!=tuple(hu.get('NAXIS'+str(i)) for i in (1,2)): raise ValueError('image/uncertainty schema')
        if not hx.get('BUNIT') or hx.get('BUNIT')!=hu.get('BUNIT'): raise ValueError('uncertainty units unresolved')
        # Header identifies an uncertainty array, but does not certify its probability distribution.
        if x.name!='PRIMARY DATA ARRAY' or u.name!='UNCERTAINTY ARRAY': raise ValueError('extension meaning')
        w=WCS(hx,naxis=2).to_header(relax=True)
        if not all(k in w for k in ('CTYPE1','CTYPE2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2')): raise ValueError('incomplete WCS')
        observer=tuple((k,hx.get(k)) for k in ('HGLN_OBS','HGLT_OBS','DSUN_OBS','RSUN_REF','DATE-OBS'))
        # Compare celestial WCS, not timestamp; observer data are retained for review.
        signatures.append((shape,hx['BUNIT'],w.tostring(sep='\n',endcard=False,padding=False)))
        images.append((x,u,observer))
    if signatures[0]!=signatures[1]: raise ValueError('matching WCS/units not established')
    return images

def characterize(paths):
    import numpy as np
    from astropy.io import fits
    with ExitStack() as stack:
        hdus=[stack.enter_context(fits.open(p,memmap=False,lazy_load_hdus=False)) for p in paths]
        pairs=inspect_headers(hdus) # all headers checked before either .data access
        arrays=[]
        for x,u,_ in pairs:
            a=np.asarray(x.data,dtype=float); b=np.asarray(u.data,dtype=float)
            if a.ndim!=2 or a.shape!=b.shape: raise ValueError('decoded shape mismatch')
            arrays.append((a,b))
        (a,ua),(b,ub)=arrays
        finite=np.isfinite(a)&np.isfinite(b)&np.isfinite(ua)&np.isfinite(ub)
        positive=np.isfinite(ua)&np.isfinite(ub)&(ua>0)&(ub>0)
        valid=finite&positive
        if not valid.any(): raise ValueError('no finite pixels with positive uncertainties')
        diff=b[valid]-a[valid]
        return {'mode':'descriptive','units':str(pairs[0][0].header['BUNIT']),'shape':list(a.shape),'total_pixels':int(a.size),'common_finite_pixels':int(finite.sum()),'positive_uncertainty_pixels':int(positive.sum()),'eligible_pixels':int(valid.sum()),'difference_median':float(np.median(diff)),'difference_p05_p95':[float(v) for v in np.percentile(diff,[5,95])], 'observer_headers':[list(pair[2]) for pair in pairs], 'limits':'Finite positive-uncertainty filter is not a validated quality mask. Shared background, calibration and registration remain unqualified; no significance.'}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--manifest',required=True); p.add_argument('--manifest-sha256',required=True); p.add_argument('--authorization',required=True); p.add_argument('--authorization-sha256',required=True); p.add_argument('--input',nargs=2,required=True); p.add_argument('--output',required=True); args=p.parse_args()
    # No receipt is written for failed authorization, which could otherwise write arbitrary paths.
    m=authorize(args.manifest,args.authorization,args.manifest_sha256,args.authorization_sha256,args.input,args.output)
    try:
        import numpy, astropy
        if numpy.__version__!=m['numpy'] or astropy.__version__!=m['astropy']: raise ValueError('dependency version mismatch')
        result={'status':'completed','result':characterize(args.input)}
    except Exception as exc:
        result={'status':'inconclusive','error_type':type(exc).__name__,'reason':str(exc)}
    result.update({'manifest_sha256':args.manifest_sha256,'authorization_sha256':args.authorization_sha256})
    with open(args.output,'x',encoding='utf-8',newline='\n') as f: json.dump(result,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(result,allow_nan=False))
    if result['status']!='completed': raise SystemExit(2)
if __name__=='__main__': main()
