"""Measure the quiescent sky at ZTF aperture centers in archival Pan-STARRS r deep stacks."""
import gzip, hashlib, io, json, math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from astropy.io import ascii,fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np

OUT=Path(__file__).parent/"image-search"/"static-reference"
OUT.mkdir(exist_ok=True)
receipts=[]
def get(base,params):
    url=base+"?"+urlencode(params)
    rec={"url":url}
    receipts.append(rec)
    with urlopen(Request(url,headers={"User-Agent":"Open-Space-static-reference/1.0 (+https://github.com/dovvnloading/Open-Space)"}),timeout=60) as r:
        raw=r.read()
        rec.update(bytes=len(raw),status=r.status,date=r.headers.get("Date"),sha256=hashlib.sha256(raw).hexdigest())
    return raw

def measure(ra,dec,zp):
    first="https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    raw=get(first,{"ra":ra,"dec":dec,"filters":"r","type":"stack"})
    tab=ascii.read(raw.decode(),format="basic")
    if not len(tab):raise ValueError("No covering stack")
    filename=str(tab[0]["filename"])
    raw=get("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi",{"ra":ra,"dec":dec,"format":"fits","red":filename,"size":720})
    # MAST's fitscut.cgi converts asinh skycell data to linear flux.
    with fits.open(io.BytesIO(raw)) as h:
        hdr=h[0].header
        dat=np.asarray(h[0].data,dtype=float)
        w=WCS(hdr).celestial
        exposure=float(hdr["EXPTIME"])
        scale=float(np.mean(proj_plane_pixel_scales(w))*3600)
        x,y=map(float,w.all_world2pix(ra,dec,0))
        yy,xx=np.mgrid[:dat.shape[0],:dat.shape[1]]
        def ap(x0,y0,radius):
            d=np.hypot(xx-x0,yy-y0)*scale
            aa=dat[d<=radius];bg=dat[(d>=11.2)&(d<=15.3)]
            if np.mean(np.isfinite(aa))<.95 or np.mean(np.isfinite(bg))<.95:return None
            return float(np.nansum(aa-np.nanmedian(bg)))/exposure
        aperture=[]
        for r in (3.55,6.1):
            value=ap(x,y,r)
            controls=[ap(x+dx/scale,y+dy/scale,r) for dx in (-45,-30,-15,0,15,30,45) for dy in (-45,-30,-15,0,15,30,45) if dx or dy]
            controls=np.array([v for v in controls if v is not None])
            if value is None or len(controls)<20:raise ValueError("Insufficient unmasked reference pixels")
            ctr=np.median(controls);sigma=1.4826*np.median(abs(controls-ctr))
            aperture.append({"radius_arcsec":r,"counts_per_second":value,"estimated_PanSTARRS_r_mag":float(25-2.5*np.log10(value)) if value>0 else None,
                             "rough_ZTF_equivalent_DN_same_r":value*10**(.4*(zp-25)),"local_control_sigma_per_second":float(sigma),
                             "controls_ge":int(sum(controls>=value)),"ncontrols":len(controls)})
        return {"filename":filename,"ra":ra,"dec":dec,"exptime":exposure,"pixscale":scale,
                "shape":list(dat.shape),"ztf_exposure_zeropoint":zp,"apertures":aperture}
known=json.loads((OUT.parent/"stack"/"stack-results.json").read_text())["frames"]
points=[("20251025541736",167.839372253,-11.366983117,26.048378)]+[(f["file"],f["ra"],f["dec"],f["zp"]) for f in known if f["file"] in ("20251027542164","20251027547164","20251029546134")]
results={}
for k,a,b,z in points:
    try:
        results[k]=measure(a,b,z)
    except Exception as exc:
        err=str(exc)
        receipts[-1]["error"]=err
        results[k]={"error":err}
    print("PS1 SKY",k,results[k],flush=True)
(OUT/"results.json").write_text(json.dumps({"positions":results,"receipts":receipts},indent=2,allow_nan=False)+"\n")
