"""Guarded, descriptive two-CTM development; no prospective image access."""
import argparse, hashlib, json, math, sys
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
    images=[]; signatures=[]; wcs=[]
    for h in hdus:
        x=h['PRIMARY DATA ARRAY']; u=h['UNCERTAINTY ARRAY']; hx=x.header; hu=u.header
        shape=tuple(hx.get('NAXIS'+str(i)) for i in (1,2))
        if hx.get('NAXIS')!=2 or hu.get('NAXIS')!=2 or any(not isinstance(v,int) or v<=0 for v in shape) or shape!=tuple(hu.get('NAXIS'+str(i)) for i in (1,2)): raise ValueError('image/uncertainty schema')
        if not hx.get('BUNIT') or hx.get('BUNIT')!=hu.get('BUNIT'): raise ValueError('uncertainty units unresolved')
        # Header identifies an uncertainty array, but does not certify its probability distribution.
        if x.name!='PRIMARY DATA ARRAY' or u.name!='UNCERTAINTY ARRAY': raise ValueError('extension meaning')
        w=WCS(hx,naxis=2).to_header(relax=True)
        if not all(k in w for k in ('CTYPE1','CTYPE2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2')): raise ValueError('incomplete WCS')
        if WCS(hx,naxis=2).has_distortion: raise ValueError('unsupported distortion')
        if (hx.get('CTYPE1'),hx.get('CTYPE2')) != ('HPLN-ARC','HPLT-ARC'): raise ValueError('unsupported coordinate frame')
        observer=tuple((k,hx.get(k)) for k in ('HGLN_OBS','HGLT_OBS','DSUN_OBS','RSUN_REF','DATE-OBS'))
        # Compare celestial WCS, not timestamp; observer data are retained for review.
        signatures.append((shape,hx['BUNIT'],hx))
        wcs.append(WCS(hx,naxis=2))
        images.append((x,u,observer))
    if signatures[0][:2]!=signatures[1][:2]: raise ValueError('shape or units mismatch')
    h0,h1=signatures[0][2],signatures[1][2]
    for key,tol in (('HGLN_OBS',0.001),('HGLT_OBS',0.001),('DSUN_OBS',1e6),('RSUN_REF',0.0)):
        a,b=h0.get(key),h1.get(key)
        if a is None or b is None or not math.isfinite(float(a)) or not math.isfinite(float(b)) or abs(float(a)-float(b))>tol: raise ValueError('observer mismatch: '+key)
    return images,wcs

def project(source, valid, w0, w1, shape):
    """Map target centers to input centers; require all four stencil pixels valid."""
    import numpy as np
    out=np.full(shape,np.nan); height,width=source.shape
    for start in range(0,shape[0],128):
        yy,xx=np.mgrid[start:min(start+128,shape[0]),0:shape[1]]
        lon,lat=w0.all_pix2world(xx,yy,0); xp,yp=w1.all_world2pix(lon,lat,0)
        good=np.isfinite(xp)&np.isfinite(yp)&(xp>=0)&(yp>=0)&(xp<width-1)&(yp<height-1)
        ix=np.floor(np.where(good,xp,0)).astype(int); iy=np.floor(np.where(good,yp,0)).astype(int)
        good &= valid[iy,ix]&valid[iy,ix+1]&valid[iy+1,ix]&valid[iy+1,ix+1]
        dx=xp-ix; dy=yp-iy
        values=(1-dx)*(1-dy)*source[iy,ix]+dx*(1-dy)*source[iy,ix+1]+(1-dx)*dy*source[iy+1,ix]+dx*dy*source[iy+1,ix+1]
        out[start:start+len(yy)]=np.where(good,values,np.nan)
    return out

def characterize(paths):
    import numpy as np
    from astropy.io import fits
    with ExitStack() as stack:
        hdus=[stack.enter_context(fits.open(p,memmap=False,lazy_load_hdus=False)) for p in paths]
        pairs,wcs=inspect_headers(hdus)
        arrays=[]
        for x,u,_ in pairs:
            a=np.asarray(x.data,dtype=float); packed=np.asarray(u.data,dtype=float)
            if a.ndim!=2 or a.shape!=packed.shape: raise ValueError('decoded shape mismatch')
            with np.errstate(divide='ignore',invalid='ignore'): sigma=a/packed
            mask=np.isfinite(a)&np.isfinite(packed)&(packed!=0)&np.isfinite(sigma)&(sigma>0)
            arrays.append((a,mask))
        (a,ma),(b,mb)=arrays
        projected=project(b,mb,wcs[0],wcs[1],a.shape)
        valid=ma&np.isfinite(projected)
        if not valid.any(): raise ValueError('no valid common pixels')
        diff=projected[valid]-a[valid]
        # Fixed 64x64 block grid; entries without valid samples are null.
        grid=[]
        for ys in np.array_split(np.arange(a.shape[0]),64):
            row=[]
            for xs in np.array_split(np.arange(a.shape[1]),64):
                block=np.ix_(ys,xs); mask=valid[block]; values=(projected[block]-a[block])[mask]
                row.append(float(np.mean(values)) if values.size else None)
            grid.append(row)
        return {'mode':'descriptive','units':str(pairs[0][0].header['BUNIT']),'shape':list(a.shape),'total_pixels':int(a.size),'native_valid_pixels':[int(ma.sum()),int(mb.sum())],'eligible_pixels':int(valid.sum()),'difference_grid_64':grid,'difference_median':float(np.median(diff)),'difference_p05_p95':[float(v) for v in np.percentile(diff,[5,95])], 'observer_headers':[list(pair[2]) for pair in pairs], 'limits':'Packed reciprocal fractional uncertainty decoded as data/packed. Four-pixel valid stencil, not validated quality mask. Angular registration only; shared background, observer parallax and correlated errors unqualified; no significance.'}

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
