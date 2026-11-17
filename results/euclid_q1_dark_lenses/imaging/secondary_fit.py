"""Independent image models and a falsifiable compact-knot SIS prediction for the H M/L lead."""
import argparse
import csv
import hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares,nnls
from scipy.signal import fftconvolve
from cutouts import CutoutArchive,image_extensions
from measure import HERE,table,stamp,blank_measurement,ab

ID="102158888_2659753895653871052"


def profiles(t,model,x,y):
    dx,dy=x-t[0],y-t[1]
    u=dx*np.cos(t[5])+dy*np.sin(t[5]);v=-dx*np.sin(t[5])+dy*np.cos(t[5])
    rr=np.sqrt(u*u+(v/t[4])**2+1e-6)
    bn=2*t[3]-1/3+4/(405*t[3])
    fg=np.exp(np.clip(-bn*((rr/t[2])**(1/t[3])-1),-80,15))
    dx=x-t[6];dy=y-t[7];rad=np.maximum(np.hypot(dx,dy),1e-6)
    if model=="SIS+shear":
        bx=dx-t[8]*dx/rad-t[14]*dx-t[15]*dy-t[9]
        by=dy-t[8]*dy/rad-t[15]*dx+t[14]*dy-t[10]
        u=bx*np.cos(t[13])+by*np.sin(t[13]);v=-bx*np.sin(t[13])+by*np.cos(t[13])
        src=np.exp(-.5*((u/t[11])**2+(v/t[12])**2))
    else:
        # Compact knot or a curved arc depending on radial and azimuthal widths.
        angle=np.arctan2(dy,dx)-t[9]
        ang=np.arctan2(np.sin(angle),np.cos(angle))
        src=np.exp(-.5*(((rad-t[8])/t[10])**2+(ang/t[11])**2))
    return fg,src


def fit_models(im,noise,psf):
    yy,xx=np.indices(im.shape);xx=xx-(im.shape[1]-1)/2;yy=yy-(im.shape[0]-1)/2
    use=(np.hypot(yy,xx)<21)&np.isfinite(im)&np.isfinite(noise)&(noise>0)&(noise<1e15)
    dat=im[use];w=1/noise[use]
    def fit(model,initial):
        lower=[-7,-7,1,.4,.18,-3.15,-8,-8,2];upper=[7,7,25,5,1,3.15,8,8,20]
        if model=="SIS+shear":
            lower += [-18,-18,.4,.4,-3.15,-.3,-.3];upper += [18,18,8,8,3.15,.3,.3]
        else:
            lower += [-6.3,.6,.08];upper += [6.3,6,2.2]
        def render(t):
            lens,src=0.,0.
            for sx in (-.25,.25):
                for sy in (-.25,.25):
                    a,b=profiles(t,model,xx+sx,yy+sy);lens+=a/4;src+=b/4
            f=fftconvolve(lens,psf,mode="same");s=fftconvolve(src,psf,mode="same")
            m=np.column_stack((f[use],s[use]));centre=np.average(m,axis=0,weights=w*w);dmean=np.average(dat,weights=w*w)
            amps,_=nnls((m-centre)*w[:,None],(dat-dmean)*w)
            sky=dmean-centre@amps
            return f*amps[0],s*amps[1],sky
        def err(t):
            f,s,b=render(t)
            return ((f+s+b)[use]-dat)*w
        best=least_squares(err,initial,bounds=(lower,upper),x_scale="jac",max_nfev=140)
        f,s,b=render(best.x)
        return best.x,2*best.cost,f,s,b,best.nfev,best.success
    output={}
    for model in ("SIS+shear","intrinsic_arc"):
        starts=[]
        fg=[1,0,8,1.2,.5,-.3]
        if model=="SIS+shear":
            for rad,offset in [(9,-4),(7,-7),(11,-2)]:
                starts.append(fg+[2,0,rad,offset,1,2,1,.6,.05,.05])
        else:
            for phi in (-2.7,-3.4):starts.append(fg+[2,0,10,phi,1.5,.55])
        output[model]=min((fit(model,t) for t in starts),key=lambda row:row[1])
    return output,use,xx,yy


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("mass_csv",type=Path)
    p.add_argument("--index-cache",type=Path);args=p.parse_args()
    rec=next(r for r in csv.DictReader(args.mass_csv.open()) if r["id_str"]==ID)
    arc=CutoutArchive("lens",args.index_cache)
    names=[ID+".fits","result/mge_lens_light.fits","result/mge_source_light.fits","result/sersic_lens_light.fits"]
    files=arc.members(ID,names)
    table(HERE/"secondary_inputs.csv",[dict(member=f"lens/{ID}/{k}",sha256=hashlib.sha256(v).hexdigest(),bytes=len(v)) for k,v in files.items()])
    raw,light,source,sersic=[image_extensions(files[n]) for n in names]
    h,im=raw["VIS_FLUX"];noise=raw["VIS_RMS"][1];q=raw["VIS_PSF"][1];q=q/q.sum()
    crop=(slice(125,175),slice(125,175));data=im[crop].astype(float)
    models,valid,xx,yy=fit_models(data,noise[crop],q)
    fits=[];draw=[("VIS data",data)]
    for model,(t,chi,fg,src,sky,iterations,success) in models.items():
        region=(yy-t[7])*( -.62)+(xx-t[6])*(-1.06)<0
        local=region&(np.hypot(xx-t[6],yy-t[7])>2)&valid
        peak=np.unravel_index(np.argmax(np.where(local,src,-np.inf)),src.shape)
        fits.append(dict(model=model,radius_arcsec=.1*t[8],centre_x_arcsec=.1*t[6],centre_y_arcsec=.1*t[7],
                         light_mass_offset_arcsec=.1*np.hypot(t[6]-t[0],t[7]-t[1]),residual_diagonal_chi2=chi,
                         usable_pixels=valid.sum(),n_free=len(t)+3,nfev=iterations,converged=bool(success),
                         opposite_peak_vs_noise=src[peak]/noise[crop][peak],opposite_peak_frac=src[peak]/src.max(),
                         params=" ".join(f"{z:.8g}" for z in t)))
        draw+= [(model+" model",fg+src+sky),(model+" residual",data-fg-src-sky)]
    table(HERE/"secondary_models.csv",fits)
    stamp(draw,ID+"_models",HERE)
    # Benchmark uses only the published effective radius, mass centre, and the
    # measured position of the brightest VIS knot, then predicts the opposite side.
    centre=np.array([(im.shape[0]-1)/2+10*float(rec["centre_0_median_pdf"]),(im.shape[1]-1)/2+10*float(rec["centre_1_median_pdf"])])
    theta=10*float(rec["einstein_radius_effective_median_pdf"]);full_y,full_x=np.indices(im.shape)
    def pad(l):
        out=np.zeros(im.shape,dtype=float);oy,ox=(np.array(im.shape)-l.shape)//2
        out[oy:oy+l.shape[0],ox:ox+l.shape[1]]=l
        return out
    r_im=im-pad(light["VIS"][1]);rr=np.hypot(full_y-centre[0],full_x-centre[1])
    smooth=gaussian_filter(r_im,.6)
    pri=np.array(np.unravel_index(np.argmax(np.where((rr>.65*theta)&(rr<1.42*theta),smooth,-np.inf)),im.shape))
    offset=pri-centre;rplus=np.linalg.norm(offset);rminus=2*theta-rplus
    counter=centre-offset*rminus/rplus
    out=[]
    for band in ["VIS","NIR_H"]:
        h,im=raw[band+"_FLUX"];err=raw[band+"_RMS"][1]
        weights=np.isfinite(err)&(err>0)&(err<1e15)&np.isfinite(im)
        for label,l in [('MGE',light[band][1]),('Sersic',sersic.get(band,(None,None))[1])]:
            if l is None:continue
            res=im-pad(l)
            for radius in (.20,.30,.40):
                vals=[]
                for y,x in (pri,counter):
                    dy,dx,bg,sd,_=blank_measurement(im,err,weights,y,x,radius*10)
                    ind=(round(y)+dy,round(x)+dx)
                    vals.append((res[ind].sum()-bg,max(sd,np.sqrt((err[ind].astype(float)**2).sum()))))
                (pflux,sp),(cflux,sc)=vals;ratio=rminus/rplus
                out.append(dict(band=band,light_model=label,aperture_arcsec=radius,theta_E_arcsec=.1*theta,
                            outer_radius_arcsec=.1*rplus,inner_radius_arcsec=.1*rminus,flux_ratio_sis=ratio,
                            primary_counts=pflux,primary_snr=pflux/sp,counter_counts=cflux,counter_snr=cflux/sc,
                            predicted_counter_snr=ratio*pflux/sc,difference_sigma=(ratio*pflux-cflux)/np.hypot(sc,ratio*sp)))
    table(HERE/"secondary_counterimage.csv",out)
    h,im=raw["NIR_H_FLUX"];err=raw["NIR_H_RMS"][1]; y,x=centre
    dy,dx,bg,sd,_=blank_measurement(im,err,np.isfinite(im)&(err>0)&(err<1e15),y,x,theta)
    ind=(round(y)+dy,round(x)+dx);all_light=float(im[ind].sum()-bg)
    net=all_light-pad(source["NIR_H"][1])[ind].sum()
    results=[]
    for name,f in [('all light',all_light),('raw minus MGE source',net),('MGE light map',pad(light['NIR_H'][1])[ind].sum()),('Sersic light map',pad(sersic['NIR_H'][1])[ind].sum())]:
        results.append(dict(model=name,H_counts=f,H_ap_mag=ab(float(h["MAGZERO"]),f),ML_relative_to_subtracted=net/f))
    table(HERE/"secondary_H_deblend.csv",results)


if __name__=="__main__":main()
