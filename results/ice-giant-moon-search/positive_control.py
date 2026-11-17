"""Ephemeris-directed recovery of S/2025 U1 in the ten discovery exposures.

Use an excerpt of URA184 spanning emission times (moon_spk.py) and a cache
directory. Products are public; only the needed byte ranges are downloaded.
Output includes the byte ranges and SHA-256 digests needed to replay the test.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
from pathlib import Path
import time
from urllib.request import Request, urlopen
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates
import warnings
from geometry import query, direction, tangent_basis
from moon_spk import offsets

ROOTS = [f"jw06379001001_021{activity:02d}_{i:05d}_nrcb1"
         for activity in (1, 3) for i in range(1, 6)]
RADIUS = 24
HALF_FIELD = 218
USED = set()


def ranges(cache, filename, start, length):
    path = cache / f"{filename}.{start}.{length}.bin"
    USED.add(path)
    if not path.exists():
        end = start + length - 1
        url = (f"https://stpubdata.s3.amazonaws.com/jwst/public/{filename[:7]}/"
               + filename.split("_")[0] + "/" + filename + f"?bytes={start}-{end}")
        for attempt in range(4):
            try:
                with urlopen(Request(url, headers={"Range":f"bytes={start}-{end}"}), timeout=75) as r:
                    if r.status != 206 or not r.headers.get("Content-Range","").startswith(f"bytes {start}-{end}/"):
                        raise IOError("range not honored")
                    data = r.read(length+1)
                if len(data) != length:
                    raise IOError("truncated range")
                temp = path.with_suffix(".part")
                temp.write_bytes(data)
                temp.replace(path)
                break
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
    data = path.read_bytes()
    if len(data) != length:
        raise IOError(f"Invalid cache entry {path.name}")
    return data


def metadata(cache, root):
    rfile = root + "_rateints.fits"
    cfile = root + "_cal.fits"
    part = io.BytesIO(ranges(cache, rfile, 0, 65536))
    primary, sci = fits.Header.fromfile(part), fits.Header.fromfile(part)
    nio, ny, nx = (sci[f"NAXIS{i}"] for i in (3,2,1))
    start = part.tell()
    step = ((sci.data_size+2879)//2880)*2880
    # Walk successive FITS headers. All ancillary image headers fit one block.
    p = start + step
    ext = {}
    while "INT_TIMES" not in ext or "DQ" not in ext:
        head = fits.Header.fromfile(io.BytesIO(ranges(cache, rfile, p, 2880)))
        ext[head["EXTNAME"]] = (p, head)
        p += 2880 + ((head.data_size+2879)//2880)*2880
    p, h = ext["INT_TIMES"]
    raw = ranges(cache, rfile, p, 2880 + ((h.data_size+2879)//2880)*2880)
    table = fits.BinTableHDU.readfrom(io.BytesIO(raw)).data
    header = io.BytesIO(ranges(cache, cfile, 0, 65536))
    fits.Header.fromfile(header)
    cs = fits.Header.fromfile(header)
    dq_off, dq = ext["DQ"]
    assert dq["BITPIX"]==32 and sci["BITPIX"]==-32 and dq["NAXIS3"]==nio
    return dict(file=rfile, sci=start, dq=dq_off+2880, nx=nx, ny=ny,
                times=table['int_mid_MJD_UTC'].copy(), wcs=WCS(cs, naxis=2),
                mt=(cs["MT_RA"], cs["MT_DEC"]), sec=primary['EFFINTTM'])


def predictions(meta, kernel, horizons):
    mjd = np.concatenate([m['times'] for m in meta])
    planet = query(799, "JWST", mjd + 2400000.5, horizons)
    xyz = offsets(kernel, Time(mjd,format="mjd",scale="utc").tdb.jd, planet)
    los,east,north = tangent_basis(planet[:,0],planet[:,1])
    # Differential tangent offsets refer to the actual observer at each epoch.
    xi = (xyz*east).sum(axis=1)/(xyz*los).sum(axis=1)
    eta = (xyz*north).sum(axis=1)/(xyz*los).sum(axis=1)
    jobs = []
    for m in meta:
        l,e,n = tangent_basis(*m["mt"])
        for idx,t in enumerate(m['times']):
            j=len(jobs)
            v = l + e*xi[j] + n*eta[j]
            coords = np.rad2deg([np.arctan2(v[1],v[0]), np.arctan2(v[2],np.hypot(v[1],v[0]))])
            x,y = m['wcs'].all_world2pix(*coords,0)
            jobs.append((m,idx,t,float(x),float(y)))
    return jobs


def wide_images(cache, m):
    x,y = m['wcs'].all_world2pix(*m['mt'],0)
    x0,y0 = int(x)-HALF_FIELD,int(y)-HALF_FIELD
    width = 2*HALF_FIELD+1
    n = width*m['nx']*4
    def read(i):
        layers=[]
        for ext,dtype in (('sci','>f4'),('dq','>i4')):
            offset = m[ext]+(i*m['ny']+y0)*m['nx']*4
            data = ranges(cache,m['file'],offset,n)
            a = np.frombuffer(data,dtype).reshape(width,-1)[:,x0:x0+width]
            layers.append(a)
        sci,dq = layers
        return np.where((dq & (1|2|4|1024))==0,sci,np.nan)
    with ThreadPoolExecutor(10) as pool:
        images=np.array(list(pool.map(read,range(len(m['times'])))))
    # All integrations of one exposure track the planet at the same dither.
    # The rings are fixed in detector coordinates while each moon moves.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",message="All-NaN slice encountered")
        template=np.nanmedian(images,axis=0)
    return images-template,x0,y0


def patches_for_positions(residuals, x0, y0, positions):
    patches={key:[] for key in positions}
    pp=np.arange(-RADIUS,RADIUS+1)
    for i,frame in enumerate(residuals):
        good=np.isfinite(frame).astype(float)
        data=np.nan_to_num(frame)
        def smooth(sigma):
            return gaussian_filter(data,sigma)/np.maximum(gaussian_filter(good,sigma),0.01)
        high_pass=smooth(0.9)-smooth(4.)
        for body,tracks in positions.items():
            x,y=tracks[i]
            yy,xx=np.meshgrid(pp+y-y0,pp+x-x0,indexing='ij')
            image=map_coordinates(high_pass,[yy,xx],order=1,cval=np.nan)
            valid=map_coordinates(good,[yy,xx],order=1,mode='constant',cval=0.)
            image[valid<0.9]=np.nan
            patches[body].append(image)
    return patches


def other_positions(jobs, body, horizons):
    jd=np.array([j[2] for j in jobs])+2400000.5
    planet=query(799,"JWST",jd,horizons);sat=query(body,"JWST",jd,horizons)
    los,east,north=tangent_basis(planet[:,0],planet[:,1])
    v=direction(sat[:,0],sat[:,1])
    xi=(v*east).sum(axis=1)/(v*los).sum(axis=1)
    eta=(v*north).sum(axis=1)/(v*los).sum(axis=1)
    result=[]
    for (m,*_), dx,dy in zip(jobs,xi,eta):
        l,e,n=tangent_basis(*m["mt"])
        v=l+e*dx+n*dy
        ra,dec=np.rad2deg([np.arctan2(v[1],v[0]),np.arctan2(v[2],np.hypot(v[0],v[1]))])
        result.append(m["wcs"].all_world2pix(ra,dec,0))
    return np.array(result)


def summary(img):
    yy,xx=np.indices(img.shape)-RADIUS
    cen=xx*xx+yy*yy<=3**2
    ann=(xx*xx+yy*yy>=8**2)&(xx*xx+yy*yy<=20**2)
    b=np.nanmedian(img[ann]);mad=np.nanmedian(np.abs(img[ann]-b))*1.4826
    peak=np.nanargmax(np.where(cen,img,np.nan))
    py,px=np.unravel_index(peak,img.shape)
    return [(img[RADIUS,RADIUS]-b)/mad,(img[py,px]-b)/mad,px-RADIUS,py-RADIUS,mad]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache",type=Path,help="small FITS headers and timing tables")
    p.add_argument("pixels",type=Path,help="original FITS SCI/DQ byte ranges")
    p.add_argument("kernel",type=Path)
    p.add_argument("horizons",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    for directory in (a.cache,a.pixels,a.out):directory.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(10) as pool:
        meta=list(pool.map(lambda r:metadata(a.cache,r),ROOTS))
    jobs=predictions(meta,a.kernel,a.horizons)
    pos={"S2025U1":np.array([j[3:5] for j in jobs]),
         "Bianca":other_positions(jobs,708,a.horizons),"Portia":other_positions(jobs,712,a.horizons)}
    all_images={key:[] for key in pos}
    for k,m in enumerate(meta):
        data,x0,y0=wide_images(a.pixels,m)
        coords={key:v[10*k:10*k+10] for key,v in pos.items()}
        images=patches_for_positions(data,x0,y0,coords)
        for key,values in images.items():all_images[key].extend(values)
    with (a.out/"control_ephemeris.csv").open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(["body","RATEINTS","integration","mid_mjd_utc","x_pred","y_pred","central_filtered_DN_s"])
        for key,images in all_images.items():
            for (m,i,t,*_), (x,y),image in zip(jobs,pos[key],images):
                w.writerow([key,m["file"],i+1,t,round(x,4),round(y,4),image[RADIUS,RADIUS]])
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    fig,axes=plt.subplots(1,3,figsize=(7.8,3.0),constrained_layout=True)
    with (a.out/"control_summary.csv").open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(["body","subset","n_integrations","predicted_contrast_MAD","peak_within_3px_MAD","peak_dx","peak_dy","local_MAD_DN_s","contrast_at_full_peak_MAD"])
        for key,images in all_images.items():
            arr=np.array(images)
            splits={"first_five_exposures":arr[:50],"last_five_exposures":arr[50:],"all":arr}
            stacks={name:np.nanmedian(z,axis=0) for name,z in splits.items()}
            dx,dy=summary(stacks["all"])[2:4]
            for i,(name,im) in enumerate(stacks.items()):
                vals=summary(im)
                yy,xx=np.indices(im.shape)-RADIUS
                ann=(xx*xx+yy*yy>=8**2)&(xx*xx+yy*yy<=20**2)
                bg=np.nanmedian(im[ann]); z=(im[RADIUS+dy,RADIUS+dx]-bg)/vals[-1]
                w.writerow([key,name,len(splits[name]),*vals,z])
                if key=="S2025U1":
                    scale=vals[-1]
                    axes[i].imshow((im-bg)/scale,origin="lower",extent=(-24.5,24.5,-24.5,24.5),
                                   norm=Normalize(-3,8),cmap="gray",interpolation="nearest")
                    axes[i].plot(0,0,"+",color="#ef5e4c",ms=7)
                    axes[i].set_title(["First 50 integrations","Last 50 integrations","All 100 integrations"][i])
                    axes[i].set_xlabel("Detector x offset (pixels)")
                    if i==0:axes[i].set_ylabel("Detector y offset (pixels)")
    fig.suptitle("S/2025 U1 at its predicted Uranus-relative position")
    fig.savefig(a.out/"s2025u1_control.png",dpi=160)
    plt.close(fig)
    with (a.out/"control_byte_ranges.csv").open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n");w.writerow(["productFilename","offset","length","sha256"])
        for path in sorted(USED):
            name,start,length,_=path.name.rsplit(".",3)
            w.writerow([name,start,length,hashlib.sha256(path.read_bytes()).hexdigest()])


if __name__ == "__main__":
    main()
