"""Apply the same forced beam measurement along translated tracks in a 10-arcminute field."""
import argparse
from pathlib import Path
import numpy as np
from scipy.signal import fftconvolve
from replay import Image, predicted, csv_out


def matched_filter(im, sky):
    # All control offsets are integer pixel translations, giving identical design matrices.
    x0,y0=im.wcs.world_to_pixel_values(*sky)
    cx,cy=round(float(x0)),round(float(y0))
    yy,xx=np.mgrid[-9:10,-9:10]
    dx,dy=xx+cx-x0,yy+cy-y0
    use=(np.abs(dx)<9)&(np.abs(dy)<9)
    X=np.column_stack([im.kernel(dx[use],dy[use]),np.ones(use.sum()),dx[use],dy[use]])
    k=np.zeros(xx.shape)
    k[use]=np.linalg.pinv(X)[0]
    filt=fftconvolve(np.nan_to_num(im.image),k[::-1,::-1],mode="same")
    return cx,cy,filt


def main(directory,out):
    imgs=[Image(p) for p in sorted(directory.glob("*VLASS*.fits"))]
    if len(imgs)!=4:
        raise ValueError("Four 601-pixel FITS images required")
    offsets=[(x,y) for x in range(-270,271,6) for y in range(-270,271,6) if np.hypot(x,y)>15]
    cols=[]
    for im in imgs:
        sky=predicted(im.date)
        cx,cy,filt=matched_filter(im,sky)
        amplitudes=np.array([filt[cy+y,cx+x] for x,y in offsets])
        assert abs(filt[cy,cx]-im.fit(sky))<1e-9
        cols.append(amplitudes)
    controls=np.array(cols)
    target=np.array([im.fit(predicted(im.date)) for im in imgs])
    med=np.median(controls,axis=1)
    mad=1.4826*np.median(abs(controls-med[:,None]),axis=1)
    weights=1/mad[:3]**2
    stat=np.average(controls[:3],axis=0,weights=weights)
    a=np.average(target[:3],weights=weights)
    rows=[]
    for i in range(4):
        rows.append(dict(test=f"epoch{i+1}",n_controls=len(offsets),target_mjy=target[i],median_mjy=med[i],
            mad_mjy=mad[i],n_ge=int(sum(controls[i]>=target[i]))))
    rows.append(dict(test="prior_epochs_inverse_variance",n_controls=len(offsets),target_mjy=a,
        median_mjy=np.median(stat),mad_mjy=1.4826*np.median(abs(stat-np.median(stat))),n_ge=int(sum(stat>=a))))
    joint=(controls[0]>=target[0]) & (controls[2]>=target[2]) & (controls[3]>=target[3])
    print(rows,"joint",int(joint.sum()),"n",len(offsets))
    csv_out(out/"field-ranks.csv",rows)
    np.savetxt(out/"translated-tracks.csv.gz",np.column_stack([offsets,controls.T]),delimiter=",",
        fmt=["%d","%d","%.6f","%.6f","%.6f","%.6f"],header="dx_pix,dy_pix,epoch1_mjy,epoch2_mjy,epoch3_mjy,epoch4_mjy",comments="")
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    main(args.directory,args.output)
