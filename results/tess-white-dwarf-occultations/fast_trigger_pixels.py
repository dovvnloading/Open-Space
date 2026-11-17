"""Check every fast-cadence discovery box against native 20-s target pixels and off-time windows."""
import argparse
import csv
import gzip
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from ffi_acquire import propagated, tsv_rows
from ffi_images import event_masks
from replication import direct
from roster_pixels import dilution, field
from short_pixels import acquire
from wd_search import mad


def relative_levels(t, values, epoch, width, indices):
    middle, left, right = indices
    window = max(.12,5*width)
    lo, a, b, hi = np.searchsorted(t,[epoch-window,epoch-2*width,epoch+2*width,epoch+window])
    local = np.median(np.vstack((values[lo:a],values[b:hi])),axis=0)
    if not np.isfinite(local).all() or np.any(local<=0):
        return None
    center, before, after = (np.mean(values[i],axis=0)/local for i in (middle,left,right))
    shoulders = .5*(before+after)
    return np.stack((1-center,shoulders-1,shoulders-center,before-after),axis=1)


def measure(t, numbers, values, trigger):
    epoch = trigger["epoch_btjd"]
    width = trigger["duration_min"]/1440
    indices = [np.flatnonzero(mask) for mask in event_masks(t,epoch,width,None)]
    result = dict(trigger=trigger,counts_center_left_right=[len(i) for i in indices])
    if not all(len(i) for i in indices):
        return dict(result,status="incomplete fast coverage")
    levels = relative_levels(t,values,epoch,width,indices)
    if levels is None:
        return dict(result,status="invalid reference")
    # Shift whole observed cadence-number triplets, retaining exact missing-sample patterns.
    orig = numbers[indices[0][0]]
    offsets = [numbers[i]-orig for i in indices]
    step = int(np.ceil(4*width*86400/20))+1
    grid = np.arange((numbers[0]-orig)//step,(numbers[-1]-orig)//step+1)
    grid = grid[grid!=0]
    if len(grid)>512:
        grid=grid[np.linspace(0,len(grid)-1,512,dtype=int)]
    control = []
    for displacement in grid*step:
        targets = [orig+displacement+v for v in offsets]
        ixs = [np.minimum(np.searchsorted(numbers,v),len(t)-1) for v in targets]
        if not all(np.array_equal(numbers[i],v) for i,v in zip(ixs,targets)):
            continue
        center = t[ixs[0][0]] + epoch-t[indices[0][0]]
        if center-t.min()<max(.12,5*width) or t.max()-center<max(.12,5*width):
            continue
        contrasts = relative_levels(t,values,center,width,ixs)
        if contrasts is not None:control.append(contrasts)
    sampling = "same cadence offsets"
    if len(control)<25:
        # Requiring every exact cadence survives poorly for long boxes with a few flags.
        # Keep the observed per-window counts and sample from each corresponding box.
        counts = [len(i) for i in indices]
        rng = np.random.default_rng(20419)
        fallback = []
        for displacement in grid*step:
            center = epoch+displacement*20/86400
            if center-t[0]<max(.12,5*width) or t[-1]-center<max(.12,5*width):continue
            bounds=[(-.5,.5),(-1.8,-.8),(.8,1.8)]
            windows=[np.arange(np.searchsorted(t,center+lo*width),np.searchsorted(t,center+hi*width,side="right")) for lo,hi in bounds]
            if any(len(w)<n for w,n in zip(windows,counts)): continue
            sampled=[rng.choice(w,n,replace=False) for w,n in zip(windows,counts)]
            contrasts=relative_levels(t,values,center,width,sampled)
            if contrasts is not None:fallback.append(contrasts)
        if len(fallback)>len(control):
            control=fallback
            sampling="same counts per box and flank"
    control = np.asarray(control).reshape(-1,2,4)
    metrics={}
    for i,name in enumerate(("aperture","wd_image")):
        center,excess,dip,lr=levels[i]
        vals=control[:,i,0]
        metrics[name]=dict(center_loss_vs_context=float(center),flanks_excess_vs_context=float(excess),
                           flank_minus_center=float(dip),left_minus_right=float(lr),controls=len(vals),
                           control_median_center_loss=float(np.median(vals)) if len(vals) else None,
                           control_mad_sigma=float(mad(vals)) if len(vals) else None,
                           off_time_equal_or_deeper=int(np.sum(vals>=center)) if len(vals) else None)
    return dict(result,status="measured",control_sampling=sampling,**metrics)


def visit(row, star, neighbors, receipt, cache):
    ev=[item for item in row["events"] if item["trigger"]["lc"].endswith("_fast-lc.fits")]
    if not ev:return None
    names={item["trigger"]["lc"].replace("-lc.fits","-tp.fits") for item in ev}
    name,=names
    path=cache/name
    url=direct(name)
    digest=acquire(path,url)
    with fits.open(path) as h:
        if (int(h[0].header["TICID"])!=int(row["tic"]) or int(h[0].header["SECTOR"])!=row["sector"]
            or abs(h[1].header["TIMEDEL"]*86400-20)>.1):
            raise ValueError("Fast-product identity mismatch")
        d=h[1].data
        t=d["TIME"].astype(float)
        good=(d["QUALITY"]==0)&np.isfinite(t)
        year=2000+(np.median(t[good])+2457000-2451545)/365.25
        wcs=WCS(h[2].header,naxis=2)
        sx,sy=map(float,wcs.world_to_pixel_values(*propagated(star,year)))
        y,x=np.indices(h[2].data.shape)
        dist=np.hypot(x-sx,y-sy)
        ap=dist<=1.5
        if min(sx,sy,x.max()-sx,y.max()-sy)<2.5:
            raise ValueError("Truncated science disk")
        out_of_window=good.copy()
        for item in ev:
            e=item["trigger"]
            out_of_window &= abs(t-e["epoch_btjd"])>2*e["duration_min"]/1440
        typical=np.nanmedian(d["FLUX"][np.flatnonzero(out_of_window)[::100]].astype(float),axis=0)
        disk=(dist<=2.5)&np.isfinite(typical)
        mat=np.stack((typical[disk],np.ones(disk.sum()),x[disk]-sx,y[disk]-sy),axis=1)
        weight=np.linalg.pinv(mat)[0]
        cube=d["FLUX"]
        aperture=cube[:,ap].sum(axis=1,dtype=float)
        wd=cube[:,disk].astype(float)@weight
        good &= np.isfinite(aperture+wd)
        t=t[good]
        numbers=d["CADENCENO"][good].astype(int)
        values=np.stack((aperture[good],wd[good]),axis=1)
        results=[]
        for item in ev:
            result=measure(t,numbers,values,item["trigger"])
            result["two_minute_status"]=item["pixels"]["status"]
            results.append(result)
    return dict(tic=row["tic"],sector=row["sector"],source_id=row["source_id"],url=url,
                sha256=digest,bytes=path.stat().st_size,target_xy=[sx,sy],quality0_finite=len(t),
                two_minute_reference_sha256=row["sha256"],gaia_query=receipt,
                **dilution(star,neighbors,wcs,year,x[ap],y[ap],sx,sy),events=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("two_minute_results",type=Path);p.add_argument("cache",type=Path)
    p.add_argument("out",type=Path);p.add_argument("--workers",type=int,default=2)
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True)
    rows=([json.loads(f.read_text()) for f in sorted(a.two_minute_results.glob("*.json"))]
          if a.two_minute_results.is_dir() else json.load(gzip.open(a.two_minute_results)))
    parent={r["GaiaEDR3"]:r for r in tsv_rows(Path("parent.tsv"))}
    stars={r["tic"]:dict(parent[r["gaia"]],tmag=r["tmag"]) for r in csv.DictReader(open("tic_match.csv"))
           if r["status"]=="source_id_and_position"}
    snapshots=json.load(gzip.open("extension/extension_catalogue_responses.json.gz"))
    hosts={r["tic"] for r in rows}
    fields={tic:field(stars[tic],a.cache,snapshots) for tic in hosts}
    def run(row):
        return visit(row,stars[row["tic"]],*fields[row["tic"]],a.cache)
    with ThreadPoolExecutor(a.workers) as pool:
        results=[r for r in pool.map(run,rows) if r is not None]
    a.out.write_bytes(gzip.compress(json.dumps(results,allow_nan=False,separators=(",",":")).encode(),mtime=0))
