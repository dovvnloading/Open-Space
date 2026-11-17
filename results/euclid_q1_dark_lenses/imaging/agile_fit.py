"""Independent pixel fit of the near-IR-faintest new AgileLens A/B target (Q1)."""
import argparse
import hashlib
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares,lsq_linear
from scipy.signal import fftconvolve

from cutouts import CutoutArchive,image_extensions
from measure import HERE,blank_measurement,ab,table,stamp

ID="EUCLJ035520.57-490933.17"
TILE="102020532"
RA,DEC=58.83570340215,-49.15921368584
# Nearest successful published cutout on this tile, selected by separation below.
BANDS=("VIS","Y","J","H")


def obtain(band,cache,ra,dec,tile,prefix_name):
    key="VIS" if band=="VIS" else "NIR-"+band
    base=f"https://irsa.ipac.caltech.edu/ibe/data/euclid/q1/MER/{tile}/"+("VIS/" if band=="VIS" else "NISP/")
    listing=urlopen(base,timeout=60).read().decode()
    arrays,manifest={},[]
    for name,prefix in [("flux",f"EUC_MER_BGSUB-MOSAIC-{key}_TILE"),("rms",f"EUC_MER_MOSAIC-{key}-RMS_TILE")]:
        files=re.findall(r'href="('+re.escape(prefix)+r'[^"]+\.fits)"',listing)
        if len(files)!=1:raise ValueError("Ambiguous Q1 mosaic selection")
        url=base+files[0]+f"?center={ra},{dec}deg&size=30arcsec&gzip=false"
        path=cache/f"{prefix_name}_{band}_{name}.fits"
        if not path.exists():
            for attempt in range(3):
                try:
                    raw=urlopen(url,timeout=90).read()
                    image_extensions(raw)
                    path.write_bytes(raw)
                    break
                except OSError:
                    if attempt==2:raise
                    time.sleep(3*(attempt+1))
        data=path.read_bytes()
        manifest.append(dict(band=band,product=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),url=url))
        arrays[name]=next(iter(image_extensions(data).values()))
    return band,arrays,manifest


def sky_position(h,x,y):
    cr=np.deg2rad([float(h["CRVAL1"]),float(h["CRVAL2"])])
    diff=np.array([x+1-float(h["CRPIX1"]),y+1-float(h["CRPIX2"])])
    cd=np.array([[float(h["CD1_1"]),float(h["CD1_2"])],[float(h["CD2_1"]),float(h["CD2_2"])]])
    xi,eta=np.deg2rad(cd@diff)
    denom=np.cos(cr[1])-eta*np.sin(cr[1])
    ra=cr[0]+np.arctan2(xi,denom)
    dec=np.arctan2(np.sin(cr[1])+eta*np.cos(cr[1]),np.hypot(xi,denom))
    return np.rad2deg([ra,dec])


def structures(p,model,x,y):
    # Locations and widths in pixels on the native 0.1 arcsec MER grid.
    dx,dy=x-p[0],y-p[1]
    foreground=np.exp(-(dx**2+dy**2)/(2*p[3]**2))
    if model=="SIS+shear":
        rad=np.maximum(np.hypot(dx,dy),1e-7)
        bx=dx-p[2]*dx/rad-p[6]*dx-p[7]*dy-p[8]
        by=dy-p[2]*dy/rad-p[7]*dx+p[6]*dy-p[9]
        u=bx*np.cos(p[12])+by*np.sin(p[12]);v=-bx*np.sin(p[12])+by*np.cos(p[12])
        source=np.exp(-.5*((u/p[10])**2+(v/p[11])**2))
    else:
        u=dx*np.cos(p[12])+dy*np.sin(p[12]);v=-dx*np.sin(p[12])+dy*np.cos(p[12])
        rad=np.hypot(u,v/p[6]);ang=np.arctan2(v/p[6],u)
        source=np.exp(-.5*((rad-p[2])/p[10])**2+p[7]*np.cos(ang)+p[8]*np.sin(ang)+p[9]*np.cos(2*ang)+p[11]*np.sin(2*ang))
    return foreground,source


def fit_vis(im,err,psf,rguess):
    # Crop is centred on the arc seen in the released VIS stamp, ~2.2" from
    # the catalogued position. Both model families use the same pixels.
    yy,xx=np.indices(im.shape)
    cy,cx=(np.array(im.shape)-1)/2
    fit_radius=max(12,rguess*2)
    valid=(np.hypot(xx-cx,yy-cy)<=fit_radius)&np.isfinite(im)&np.isfinite(err)&(err>0)
    amp=float(np.max(im[valid]));sigma=float(np.median(err[valid]))
    results={}
    # Evaluate the surface brightness on a 3x3 subpixel grid before convolution.
    def averaged(p,model):
        fg,src=0.,0.
        for ox in (-1/3,0.,1/3):
            for oy in (-1/3,0.,1/3):
                a,b=structures(p,model,xx+ox,yy+oy);fg+=a/9;src+=b/9
        return fg,src
    for model in ("SIS+shear","intrinsic_ring"):
        if model=="SIS+shear":
            low=[11,11,2,.3,0,0,-.25,-.25,-4,-4,.20,.20,-3.15,-4*sigma]
            high=[19,19,9,4,10*amp,10*amp,.25,.25,4,4,4,4,3.15,4*sigma]
            start=[15,15,5.5,1,.3*amp,.8*amp,.01,-.02,.3,.3,1,.6,.3,0]
        else:
            low=[11,11,2,.3,0,0,.5,-3,-3,-3,.35,-3,-3.15,-4*sigma]
            high=[19,19,10,4,10*amp,10*amp,1,3,3,3,3,3,3.15,4*sigma]
            start=[15,15,6,1,.3*amp,.7*amp,.9,.1,.1,.1,1,.1,.1,0]
        low[:3]=[cx-5,cy-5,.45*rguess]
        high[:3]=[cx+5,cy+5,1.75*rguess]
        start[:3]=[cx,cy,rguess]
        def predict(p):
            fg,src=averaged(p,model)
            return fftconvolve(p[4]*fg+p[5]*src,psf,mode="same")+p[-1]
        def residual(p):return ((predict(p)-im)/err)[valid]
        candidates=[]
        for seed in (-.8,.8):
            init=np.array(start);init[8]=seed;init[12]=seed
            fit=least_squares(residual,init,bounds=(low,high),max_nfev=240,x_scale="jac")
            candidates.append(fit)
        best=min(candidates,key=lambda fit:fit.cost)
        fg,src=averaged(best.x,model)
        results[model]=(best.x,2*best.cost,int(valid.sum()),predict(best.x),fg,src)
    return results


def main():
    import csv
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("catalog",type=Path)
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--index-cache",type=Path)
    p.add_argument("--faint-vis-a",action="store_true",help="also supported: faint VIS AgileLens A, at its arc centre")
    args=p.parse_args();args.cache.mkdir(exist_ok=True,parents=True)
    tile,ra0,dec0,tag,obj=TILE,RA,DEC,"agile",ID
    reg=(slice(153,184),slice(122,153)); guess=5.5
    if args.faint_vis_a:
        tile,ra0,dec0,tag,obj="102019126",58.77131071304,-50.27601482573,"agile_faint_a","EUCLJ035505.11-501633.65"
        reg=(slice(124,177),slice(120,173));guess=10.
    archive=CutoutArchive("lens",args.index_cache)
    rows=[r for r in csv.DictReader(args.catalog.open()) if r["tile_index"]==tile and f'lens/{r["id_str"]}/{r["id_str"]}.fits' in archive.entries]
    near=min(rows,key=lambda r:(np.cos(np.deg2rad(dec0))*(float(r["right_ascension"])-ra0))**2+(float(r["declination"])-dec0)**2)
    id_=near["id_str"]
    pd=archive.members(id_,[id_+".fits"])[id_+".fits"];psf_products=image_extensions(pd)
    with ThreadPoolExecutor(4) as pool: inputs=list(pool.map(lambda b:obtain(b,args.cache,ra0,dec0,tile,tag),BANDS))
    manifest=[r for _,_,m in inputs for r in m]
    manifest.append(dict(band="all",product="same-tile PSFs",bytes=len(pd),sha256=hashlib.sha256(pd).hexdigest(),
                         url=f"https://zenodo.org/records/15025832/files/lens.zip#lens/{id_}/{id_}.fits"))
    table(HERE/(tag+"_inputs.csv"),manifest)
    data={b:a for b,a,_ in inputs}
    def psf(b):
        name="VIS" if b=="VIS" else "NIR_"+b
        q=psf_products[name+"_PSF"][1].astype(float);return q/q.sum()
    vh,vi=data["VIS"]["flux"];ve=data["VIS"]["rms"][1]
    results=fit_vis(vi[reg],ve[reg],psf("VIS"),guess)
    models,phot=[],[]
    ae=[]
    for model,(v,chi,n,pred,g,s) in results.items():
        ra,dec=sky_position(vh,reg[1].start+v[0],reg[0].start+v[1])
        models.append(dict(model=model,ra_fit=ra,dec_fit=dec,radius_arcsec=v[2]*.1,foreground_width_arcsec=v[3]*.1,
                      residual_diagonal_chi2=chi,usable_pixels=n,n_fit_parameters=len(v),parameters=" ".join(f"{x:.8g}" for x in v)))
    v,_,_,pred,fg,src=results["SIS+shear"]
    draw=[]
    for b in BANDS:
        h,im=data[b]["flux"];err=data[b]["rms"][1]
        pixx,pixy=reg[1].start+v[0],reg[0].start+v[1]
        yy,xx=np.indices(fg.shape)
        valid=(np.hypot(xx-v[0],yy-v[1])<=max(12,2*guess))&np.isfinite(im[reg])&np.isfinite(err[reg])&(err[reg]>0)
        fgconv=fftconvolve(fg,psf(b),mode="same");sc=fftconvolve(src,psf(b),mode="same")
        design=np.column_stack((fgconv[valid],sc[valid],np.ones(valid.sum())))
        w=1/err[reg][valid]
        fit=lsq_linear(design*w[:,None],im[reg][valid]*w,bounds=([0,0,-np.inf],np.inf))
        dy,dx,bg,scatter,nsky=blank_measurement(im,err,np.isfinite(im)&np.isfinite(err)&(err>0),pixy,pixx,3)
        ap=(round(pixy)+dy,round(pixx)+dx);loc=(round(v[1])+dy,round(v[0])+dx)
        raw=im[ap].sum()-bg;ss=fit.x[1]*sc[loc].sum()
        noise=max(scatter,np.sqrt(np.sum(err[ap].astype(float)**2)));zp=float(h["MAGZERO"])
        phot.append(dict(band=b,blended_core_snr=raw/noise,source_subtracted_core_snr=(raw-ss)/noise,
                   core_source_subtracted_mag=ab(zp,raw-ss),noise_only_5sigma_mag=ab(zp,5*noise),
                   core_source_fraction=ss/raw,fit_foreground_mag=ab(zp,fit.x[0]*fgconv.sum()),sky_apertures=nsky))
        dy,dx,bg,scatter,nsky=blank_measurement(im,err,np.isfinite(im)&np.isfinite(err)&(err>0)&(err<1e15),pixy,pixx,v[2])
        ap=(round(pixy)+dy,round(pixx)+dx);loc=(round(v[1])+dy,round(v[0])+dx)
        net=im[ap].sum()-bg-fit.x[1]*sc[loc].sum()
        sig=max(scatter,np.sqrt(np.sum(err[ap].astype(float)**2)))
        ae.append(dict(band=b,theta_arcsec=v[2]*.1,ap_source_sub_counts=net,sigma=sig,
                       aperture_mag=ab(zp,net),noise_only_5sigma_mag=ab(zp,5*sig)))
        if b in ("VIS","H"):
            draw.extend([(b+" data",im[reg]),(b+" minus foreground",im[reg]-fit.x[0]*fgconv)])
    table(HERE/(tag+"_models.csv"),models);table(HERE/(tag+"_core.csv"),phot)
    table(HERE/(tag+"_einstein_light.csv"),ae)
    stamp(draw,obj,HERE)
    # Same solar SED and blueshifted Euclid passband convention as inference/.
    import sys
    sys.path.insert(0,str(HERE.parent/"inference"))
    from infer import chi,mass_floor,luminosity
    sun=np.genfromtxt(HERE.parent/"inference"/"sun_shifted_band_ab.csv",delimiter=",",names=True)
    sensitivity=[]
    for zl in (.15,.3,.6,1.,1.5,2.95):
        for zs in (2.,4.,6.):
            if zs<=zl:continue
            beta=1-chi(zl)/chi(zs)
            mass=mass_floor(v[2]*.1,zl)/beta
            for r in ae:
                lum=luminosity(r["aperture_mag"],zl,sun,r["band"].lower())
                sensitivity.append(dict(band=r["band"],assumed_zl=zl,assumed_zs=zs,
                    theta_arcsec=v[2]*.1,mass_msun=mass,lum_lsun=lum,ml_msun_lsun=mass/lum))
    table(HERE/(tag+"_geometry_sensitivity.csv"),sensitivity)


if __name__=="__main__":main()
