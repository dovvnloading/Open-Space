"""PSF-matched arc-knot colours for the two conditional M/L outliers."""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter,maximum_filter
from scipy.signal import fftconvolve
from cutouts import CutoutArchive,image_extensions
from measure import HERE,aperture,ab,table

IDS=["102158887_2655306065654633189","102158888_2659753895653871052"]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("mass_csv",type=Path);p.add_argument("--index-cache",type=Path)
    args=p.parse_args();mass={r["id_str"]:r for r in csv.DictReader(args.mass_csv.open())}
    arc=CutoutArchive("lens",args.index_cache);rows=[]
    for obj in IDS:
        record=mass[obj];rad=float(record["einstein_radius_effective_median_pdf"])*10
        d=arc.members(obj,[obj+".fits","result/mge_lens_light.fits"]);raw=image_extensions(d[obj+".fits"])
        light=image_extensions(d["result/mge_lens_light.fits"]);imgs,psfs={},{}
        for band in ("VIS","NIR_Y","NIR_J","NIR_H"):
            head,im=raw[band+"_FLUX"]
            _,l=light[band];oy=(im.shape[0]-l.shape[0])//2;ox=(im.shape[1]-l.shape[1])//2
            residual=im.astype(float);residual[oy:oy+l.shape[0],ox:ox+l.shape[1]]-=l
            imgs[band]=residual
            kernel=raw[band+"_PSF"][1].astype(float);psfs[band]=kernel/kernel.sum()
        y,x=(np.array(im.shape)-1)/2+[float(record["centre_0_median_pdf"])*10,float(record["centre_1_median_pdf"])*10]
        yy,xx=np.indices(im.shape); rr=np.hypot(xx-x,yy-y)
        image=gaussian_filter(imgs["VIS"],.6)
        mask=(rr>.65*rad)&(rr<1.42*rad)&(image==maximum_filter(image,5))
        pts=sorted(zip(*np.where(mask)),key=lambda p:-image[p])
        if not pts:continue
        pts=[p for p in pts if image[p]>.20*image[pts[0]]]
        for band in ("NIR_Y","NIR_J","NIR_H"):
            # Cross-convolution puts VIS and each NIR image on exactly the
            # same effective PSF without deconvolution or a Gaussian PSF guess.
            v=fftconvolve(np.nan_to_num(imgs["VIS"]),psfs[band],mode="same")
            n=fftconvolve(np.nan_to_num(imgs[band]),psfs["VIS"],mode="same")
            dy,dx=aperture(2.3,0.,0.)
            blank=[(i,j) for i in range(25,im.shape[0]-25,10) for j in range(25,im.shape[1]-25,10) if 70<np.hypot(i-y,j-x)<125]
            for k,(py,px) in enumerate(pts,1):
                vals=[]
                for name,arr in [("VIS",v),(band,n)]:
                    sky=np.array([arr[i+dy,j+dx].sum() for i,j in blank])
                    bg=np.median(sky);sd=1.4826022*np.median(np.abs(sky-bg))
                    flux=arr[py+dy,px+dx].sum()-bg
                    vals.append((flux,sd,ab(float(raw[name+"_FLUX"][0]["MAGZERO"]),flux)))
                f1,s1,m1=vals[0];f2,s2,m2=vals[1]
                dist=np.hypot(py-y,px-x)*.1
                rows.append(dict(id_str=obj,knot=k,band=band,dx_arcsec=(px-x)*.1,dy_arcsec=(py-y)*.1,
                    radius_arcsec=dist,vis_snr=f1/s1,nir_snr=f2/s2,vis_minus_nir=m1-m2,
                    sky_only_colour_sigma=(2.5/np.log(10))*np.hypot(s1/f1,s2/f2) if min(f1,f2)>0 else np.nan))
    table(HERE/"knot_colours.csv",rows)


if __name__=="__main__":main()
