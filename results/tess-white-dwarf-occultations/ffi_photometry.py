"""Measure small apertures and a Gaia-based crowding sensitivity range on archived TESS cutouts."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.special import ndtr

from ffi_acquire import propagated, tsv_rows

def load_pixel_tools(path):
    spec=importlib.util.spec_from_file_location("existing_tess_pixels",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def match_tic(star, rows):
    ra,dec=propagated(star,2000.)
    matches=[]
    for row in rows:
        sep=3600*np.hypot((float(row["RAJ2000"])-ra)*np.cos(np.deg2rad(dec)),float(row["DEJ2000"])-dec)
        if row["GAIA"]==star["GaiaEDR3"] and sep<2 and not row.get("Disp"):
            matches.append((row,sep))
    return matches[0] if len(matches)==1 else (None,None)

def deposit(x,y,sx,sy,sigma):
    def axis(pixel,center):
        return ndtr((pixel+.5-center)/sigma)-ndtr((pixel-.5-center)/sigma)
    return np.sum(axis(x,sx)*axis(y,sy))

def gaia_t(row,worst=False):
    g=float(row["Gmag"])
    c=(float(row["BPmag"])-float(row["RPmag"])) if row["BPmag"] and row["RPmag"] else (3.5 if worst else 1.)
    c=np.clip(c,-.2,3.5)
    return g-.00522555*c**3+.0891337*c**2-.633923*c+.0324473

def run(star,root,out,px):
    source=star["GaiaEDR3"]
    sector=int(star["sector"])
    coordinate_epoch=float(star.get("coordinate_epoch") or 2026.8)
    base=f"{source}-s{sector}"
    result=dict(source_id=source,wd=star["WDJname"],sector=sector)
    path=root/(base+".fits")
    if not path.is_file():
        return dict(result,status="no cutout")
    result["fits_sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
    assert result["fits_sha256"]==json.loads(path.with_suffix(".fits.json").read_text())["sha256"]
    tic,sep=match_tic(star,tsv_rows(root/f"{source}-tic82.tsv"))
    if tic is None:
        return dict(result,status="excluded: no unique Gaia-ID and <2 arcsec TIC match")
    result.update(tic=tic["TIC"],tic_sep_arcsec=sep,tmag=float(tic["Tmag"]))
    neighbours=tsv_rows(root/f"{source}-gaia.tsv")
    if sum(r["Source"]==source for r in neighbours)!=1:
        raise ValueError("WD Gaia ID missing or duplicate in cone")
    with fits.open(path,memmap=True) as fitsfile:
        primary,ext=fitsfile[0].header,fitsfile[1].header
        if (int(primary["SECTOR"]),int(primary["CAMERA"]),int(primary["CCD"]))!=tuple(int(star[k]) for k in ("sector","camera","ccd")):
            raise ValueError("unexpected detector/sector")
        data=fitsfile[1].data
        result["arrays_sha256"]={name:hashlib.sha256(data[name].tobytes()).hexdigest() for name in ("TIME","FLUX","QUALITY")}
        result["date_obs"]=primary["DATE-OBS"]
        result["date_end"]=primary["DATE-END"]
        ny,nx=data["FLUX"].shape[1:]
        yy,xx=np.indices((ny,nx))
        tx,ty=px.target_pixel(ext,*propagated(star,coordinate_epoch))
        aperture=(xx-tx)**2+(yy-ty)**2<=1.5**2
        if nx!=9 or ny!=9 or not (2<tx<6 and 2<ty<6):
            raise ValueError("WD lies outside central cutout")
        neighbours_xy=[]
        for r in neighbours:
            if not r.get("Gmag"):
                continue
            x,y=px.target_pixel(ext,*propagated(r,coordinate_epoch))
            neighbours_xy.append((r,x,y))
        # Pixel-integrated circular Gaussian at three widths. Test a 0.3 mag
        # faintward target colour error and red colour for unidentified neighbours.
        fractions=[]
        contribution=[]
        for sigma in (.65,.85,1.05):
            throughput=deposit(xx[aperture],yy[aperture],tx,ty,sigma)
            for worst in (False,True):
                ratios=[]
                for r,x,y in neighbours_xy:
                    if r["Source"]==source:
                        continue
                    tmag=gaia_t(r,worst=worst)
                    relative=10**(-.4*(tmag-float(tic["Tmag"])-(.3 if worst else 0.)))
                    integral=deposit(xx[aperture],yy[aperture],x,y,sigma)
                    ratio=relative*integral/throughput
                    ratios.append(ratio)
                    if sigma==.85 and not worst:
                        contribution.append((r["Source"],x,y,tmag,ratio))
                value=1/(1+sum(ratios))
                fractions.append(value)
                if sigma==.85 and not worst:
                    nominal=value
        peak=sorted(contribution,key=lambda row: -row[-1])[:5]
        result.update(camera=int(primary["CAMERA"]),ccd=int(primary["CCD"]),target_xy=[tx,ty],
                      target_fraction=float(nominal),min_target_fraction=float(min(fractions)),
                      max_target_fraction=float(max(fractions)),n_gaia=len(neighbours),
                      gaia_without_G=len(neighbours)-len(neighbours_xy),neighbours_top5=peak)
        if min(fractions)<.8 or result["gaia_without_G"]:
            return dict(result,status="excluded: Gaia crowding cut")
        edge=(xx==0)|(xx==nx-1)|(yy==0)|(yy==ny-1)
        inside3=np.hypot(xx-tx,yy-ty)<3.0
        time=np.asarray(data["TIME"],dtype=float)
        cube=np.asarray(data["FLUX"],dtype=float)
        sky=np.nanmedian(cube[:,edge],axis=1)
        sub=cube-sky[:,None,None]
        aperture_flux=np.sum(sub[:,aperture],axis=1)
        errs=np.asarray(data["FLUX_ERR"][:,aperture],dtype=float)
        sky_error=1.4826*np.nanmedian(np.abs(cube[:,edge]-sky[:,None]),axis=1)*1.253/np.sqrt(edge.sum())
        flux_err=np.sqrt(np.sum(errs**2,axis=1)+(aperture.sum()*sky_error)**2)
        keep=(data["QUALITY"]==0)&np.isfinite(time+aperture_flux)
        if keep.sum()<100:
            return dict(result,status="excluded: insufficient quality=0 finite exposures")
        median=float(np.median(aperture_flux[keep]))
        if median<=0:
            return dict(result,status="excluded: nonpositive aperture baseline")
        weights=np.maximum(0,sub[:,inside3])
        denom=weights.sum(axis=1)
        cx=np.divide(weights@xx[inside3],denom,out=np.full(len(time),np.nan),where=denom>0)
        cy=np.divide(weights@yy[inside3],denom,out=np.full(len(time),np.nan),where=denom>0)
        dt=float(ext["TIMEDEL"])*86400
        if abs(dt-200)>1: raise ValueError("not 200-second FFI")
        lc=out/f"tic{tic['TIC']}-s{sector}.npz"
        np.savez_compressed(lc,time=time[keep],flux=aperture_flux[keep],flux_err=flux_err[keep],
                            cadence_s=dt,target_fraction=nominal,min_target_fraction=min(fractions),
                            centroid_x=cx[keep],centroid_y=cy[keep],background=sky[keep],
                            pixel_row=np.flatnonzero(keep),aperture=aperture,gaia_id=source,tic=tic["TIC"],sector=sector)
        return dict(result,status="included",quality0_points=int(sum(keep)),total_points=len(time),start_btjd=float(time[keep][0]),
                    end_btjd=float(time[keep][-1]),cadence_s=dt,median_flux_e_s=median,
                    lc_file=lc.name,lc_sha256=hashlib.sha256(lc.read_bytes()).hexdigest())

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("out",type=Path)
    parser.add_argument("--roster",default="s*.csv",help="input roster glob (history.csv for auxiliary visits)")
    parser.add_argument("--pixel-tools",type=Path,default=Path(__file__).resolve().parent.parent/"tess-asymmetric-dips/pixels.py")
    args=parser.parse_args()
    args.out.mkdir(exist_ok=True,parents=True)
    px=load_pixel_tools(args.pixel_tools)
    records=[]
    visits=[]
    for roster in sorted(args.input.glob(args.roster)):
        for star in csv.DictReader(roster.open()):
            result=run(star,args.input,args.out,px)
            print(star["WDJname"],result.get("tic"),result["status"],result.get("min_target_fraction"),flush=True)
            records.append(result)
            if result["status"]=="included":
                # IAU nominal solar mass parameter [m^3/s^2], radius [m], logg in cgs.
                mass=float(star["MassH"])
                radius=np.sqrt(1.3271244e20*mass/(10**float(star["loggH"])*.01))/6.957e8
                visits.append(dict(target=star["GaiaEDR3"],path=result["lc_file"],
                                   cadence_s=result["cadence_s"],light_fraction=result["target_fraction"],
                                   mass_solar=mass,radius_solar=radius))
    with (args.out/"ffi_visits.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=["target","path","cadence_s","light_fraction","mass_solar","radius_solar"])
        writer.writeheader()
        writer.writerows(visits)
    (args.out/"ffi_measurements.json").write_text(json.dumps(records,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
