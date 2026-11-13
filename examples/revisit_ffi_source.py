"""Test a frozen sky location and frequency against independent full-frame observations."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from open_space.ffi import _wave, condition, detect
from open_space.ffi_archive import acquire
from open_space.tess_archive import atomic_bytes
from open_space.tess_eb_catalog import separation


def inspect(cube, blind_sources, hypothesis):
    time, flux, error, wcs = cube
    ra, dec, p = [hypothesis[k] for k in ("ra_deg", "dec_deg", "fixed_period_days")]
    x,y = wcs.world_to_pixel_values(ra, dec)
    col, row = int(np.floor(x+0.5)), int(np.floor(y+0.5))
    if not (0 <= col < flux.shape[2] and 0 <= row < flux.shape[1]):
        return {"outcome":"off_detector"}
    series, _, nuisance, width = condition(time,flux,error,excluded=(col,row))
    design = _wave(time,p**-1)
    fit = np.linalg.lstsq(design,series,rcond=None)[0][:2]
    # Reference phase chosen from one predetermined detector pixel nearest the supplied sky position.
    k = row*width+col
    direction = fit[:,k]/max(1e-15,np.linalg.norm(fit[:,k]))
    image = (direction @ fit).reshape(flux.shape[1:])
    yy,xx = np.indices(image.shape)
    patch = (abs(xx-col)<=2)&(abs(yy-row)<=2)
    ring = (abs(xx-col)<=3)&(abs(yy-row)<=3)&~((abs(xx-col)<=1)&(abs(yy-row)<=1))
    image = image - np.median(image[ring])
    weights = np.maximum(image,0)*patch
    cx,cy = (weights*xx).sum()/weights.sum(),(weights*yy).sum()/weights.sum()
    sky = tuple(map(float,wcs.pixel_to_world_values(cx,cy)))
    signal = design[:,:2] @ fit[:,k]
    trace = series[:,k]
    nuisance_corr = np.max(abs(trace@nuisance)/np.maximum(1e-12,np.linalg.norm(trace)*np.linalg.norm(nuisance,axis=0)))
    fractions = ((time-1325)/p)%1
    phase = []
    for left in np.arange(20)/20:
        hits = (fractions >= left) & (fractions < left+.05)
        phase.append(float(np.mean(trace[hits])) if hits.any() else None)
    near = sorted(blind_sources,key=lambda d:separation(ra,dec,d["ra_deg"],d["dec_deg"]))
    closest = near[0] if near else None
    offset = separation(ra,dec,closest["ra_deg"],closest["dec_deg"]) if closest else None
    stamp = image[max(0,row-2):row+3, max(0,col-2):col+3]
    return {"outcome":"measured","nearest_blind_peak":closest,"nearest_blind_offset_arcsec":offset,
            "peak_in_matching_radius":bool(closest and offset<=31.5),
            "fixed_frequency_power":float(np.var(signal)/np.var(trace)),
            "fixed_frequency_semiamplitude_e_per_s":float(np.linalg.norm(fit[:,k])),
            "fixed_frequency_spectral_window":float(abs(np.exp(2j*np.pi*(time-time[0])/p).mean())**2),
            "field_max_abs_nuisance_correlation":float(nuisance_corr),
            "negative_to_positive_flux":float(np.maximum(-image[patch],0).sum()/weights.sum()),
            "fixed_frequency_center_icrs_deg":sky,
            "fixed_frequency_center_offset_arcsec":separation(ra,dec,*sky),
            "fixed_frequency_signal_map_e_per_s":stamp.tolist(), "x_pixel":col, "y_pixel":row,
            "phase_bins_e_per_s":phase,"phase_reference_btjd":1325.0}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan",type=Path,default=Path("data/ffi-source-revisit.json"))
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    raw=args.plan.read_bytes()
    plan=json.loads(raw)
    outcomes=[]
    for product in plan["products"]:
        cube,meta=acquire(plan["field"],product,args.cache)
        peaks,counts=detect(cube)
        outcomes.append({"product":meta,"blind_search":counts,"cases":[
            {**h,**inspect(cube,peaks,h)} for h in plan["hypotheses"]]})
    out={"schema_version":1,"input_plan_sha256":sha256(raw).hexdigest(),"sectors":outcomes}
    atomic_bytes(args.output,(json.dumps(out,indent=2,allow_nan=False)+"\n").encode())


if __name__=="__main__":
    main()
