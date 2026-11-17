"""Independent image-position lens fit on the strongest SLDE H-band priority candidate."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "imaging"))
from cutouts import CutoutArchive, image_extensions

ID = "102158891_2700857611652676803"
BANDS = ("VIS", "NIR_Y", "NIR_J", "NIR_H")
PIXSCALE = 0.1  # arcsec on common grid, checked against FITS


def inputs(cache):
    cache.mkdir(parents=True, exist_ok=True)
    raw, model = cache/(ID + ".fits"), cache/"mge_lens_light.fits"
    if not raw.exists() or not model.exists():
        data = CutoutArchive("lens", cache).members(ID, [ID + ".fits", "result/mge_lens_light.fits"])
        raw.write_bytes(data[ID + ".fits"])
        model.write_bytes(data["result/mge_lens_light.fits"])
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw, model]}
    if hashes != {ID + ".fits": "40f9b755b5d9f369caaae9df3200f9ab54faaa6562bec7fabaea02a84c04ca02",
                  "mge_lens_light.fits": "374b0b025f19668ac8c9dbe4d975902a43d00b8c91358e94e40e6ff12e92c455"}:
        raise ValueError("Unexpected input cutout or foreground snapshot")
    return image_extensions(raw.read_bytes()), image_extensions(model.read_bytes()), hashes


def centroid(images):
    data = images["VIS_FLUX"][1].astype(float)[125:175, 125:175]
    y, x = np.indices(data.shape, dtype=float)
    in_core = (x-24.5)**2 + (y-24.5)**2 < 5**2
    weights = np.maximum(data - np.median(data[in_core]), 0) * in_core
    return np.array([(weights*x).sum(), (weights*y).sum()]) / weights.sum()


def radial_foreground(data, center):
    # Independent smooth radial foreground in narrow shells. A 30th percentile
    # rejects bright arc sectors. Constant sky is fitted by the radial estimate.
    y, x = np.indices(data.shape)
    radius = np.hypot(x-center[0], y-center[1])
    bins = np.arange(0.5, radius.max()+2)
    level = [np.percentile(data[np.abs(radius-r)<0.75], 30) for r in bins
             if np.any(np.abs(radius-r)<0.75)]
    return np.interp(radius, bins[:len(level)], gaussian_filter(np.array(level), 0.6))


def peak_positions(residual, center, smoothing=0.95):
    # Identify peaks without any input lens parameters or catalogue radius.
    y, x = np.indices(residual.shape)
    radius = np.hypot(x-center[0], y-center[1]) * PIXSCALE
    smoothed = gaussian_filter(residual, smoothing)
    peaks = np.where((smoothed == maximum_filter(smoothed, 7)) & (radius>0.55) & (radius<2.15))
    peaks = sorted(zip(*peaks), key=lambda yx: -smoothed[yx])[:4]
    positions = []
    for py, px in peaks:
        patch = smoothed[py-1:py+2, px-1:px+2]
        yy, xx = np.mgrid[-1:2, -1:2]
        basis = np.column_stack([xx.flatten()**2, (xx*yy).flatten(), yy.flatten()**2,
                                 xx.flatten(), yy.flatten(), np.ones(9)])
        a,b,c,d,e,_ = np.linalg.lstsq(basis,patch.ravel(),rcond=None)[0]
        offset = np.linalg.solve([[2*a,b],[b,2*c]], [-d,-e])
        if np.max(np.abs(offset)) > 1:
            offset = np.zeros(2)
        positions.append(np.array([px,py])+offset)
    return np.array(sorted(positions, key=lambda p: np.arctan2(p[1]-center[1], p[0]-center[0])))


def fit(positions, center):
    """Solve SIS+constant shear lens equation with one shared source and fixed light centre."""
    x,y = ((positions-center)*PIXSCALE).T
    r = np.hypot(x,y)
    # theta - b theta/r - Gamma theta = beta. Linear in b, g1, g2, beta_x, beta_y.
    A = np.vstack([np.column_stack([x/r, x, y, np.ones(4), np.zeros(4)]),
                   np.column_stack([y/r, -y, x, np.zeros(4), np.ones(4)])])
    b, g1, g2, bx, by = np.linalg.lstsq(A, np.r_[x,y], rcond=None)[0]
    shear = np.array([[g1, g2],[g2, -g1]])
    theta = np.column_stack([x,y])
    beta = theta - b*theta/r[:,None] - theta @ shear.T
    det = []
    corrections = []
    for i,p in enumerate(theta):
        jac = np.eye(2)-shear - b*(np.eye(2)/r[i] - np.outer(p,p)/r[i]**3)
        det.append(np.linalg.det(jac))
        corrections.append(np.linalg.solve(jac,beta[i]-[bx,by]))
    return np.array([b,g1,g2,bx,by]), np.asarray(det), np.array(corrections), beta


def colours(im, model, positions, center):
    # Apply the other three PSFs to each plane to give the same effective PSF.
    y,x = np.indices((150,150))
    rng = np.random.default_rng(312)
    xy = rng.uniform(12,137,(1200,2))
    blank = xy[np.hypot(xy[:,0]-(center[0]+50),xy[:,1]-(center[1]+50))>31]
    out = {}
    for band in BANDS:
        head,data = im[band+"_FLUX"]
        residual = data[75:225,75:225].astype(float)-model[band][1]
        for other in BANDS:
            if other != band:
                psf = im[other+"_PSF"][1].astype(float)
                residual = fftconvolve(residual, psf/psf.sum(), mode="same")
        # Sum aperture on common grid, no aperture correction; colours require equal effective PSF.
        def aperture(pos):
            px,py = pos
            iy,ix=int(round(py)),int(round(px))
            dx,dy=x[iy-5:iy+6,ix-5:ix+6]-px, y[iy-5:iy+6,ix-5:ix+6]-py
            return residual[iy-5:iy+6,ix-5:ix+6][dx*dx+dy*dy<=3.5**2].sum()
        sky = np.array([aperture(p) for p in blank])
        factor = 10**((23.9-float(head["MAGZERO"]))/2.5)
        med = np.median(sky)
        noise = 1.4826*np.median(abs(sky-med)) * factor
        values = [(aperture(p+50)-med)*factor for p in positions]
        out[band] = {"uJy": values, "empirical_sky_sigma_uJy": noise,
                     "flags": ["below_3sigma" if v < 3*noise else "measured" for v in values]}
    return out


def compare_colours(phot):
    # Scale one common SED by a separate image magnification. The first scale fixes the degeneracy.
    flux = np.array([phot[b]["uJy"] for b in BANDS]).T
    noise = np.array([phot[b]["empirical_sky_sigma_uJy"] for b in BANDS])
    def residual(parameters):
        sed = parameters[:4]
        scales = np.r_[1,parameters[4:]]
        return ((scales[:,None]*sed-flux)/noise).ravel()
    result = least_squares(residual, np.r_[flux[0], 1,1,1], bounds=(0,np.inf))
    return {"VIS_minus_H_AB": (-2.5*np.log10(flux[:,0]/flux[:,-1])).tolist(),
            "common_SED_chi2_diag": float(np.sum(result.fun**2)), "nominal_dof": 9,
            "relative_brightness": np.r_[1,result.x[4:]].tolist()}


def main(cache, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    im, model, hashes = inputs(cache)
    if abs(abs(float(im["VIS_FLUX"][0]["CD1_1"]))*3600-PIXSCALE) > 1e-5:
        raise ValueError("Pixel scale changed")
    data = im["VIS_FLUX"][1][125:175,125:175].astype(float)
    center = centroid(im)
    residual = data-radial_foreground(data, center)
    pos = peak_positions(residual,center)
    p,det,delta,beta = fit(pos,center)
    sigma_c, sigma_peak = 0.05, 0.1  # explicit registration & morphology sensitivity in arcsec
    rng = np.random.default_rng(442)
    sim = np.array([fit(pos+rng.normal(0,sigma_peak/PIXSCALE,pos.shape),
                        center+rng.normal(0,sigma_c/PIXSCALE,2))[0] for _ in range(2000)])
    cross = []
    for foreground, sub in [("radial", residual), ("released_MGE",data-model["VIS"][1][50:100,50:100])]:
        for smoothing in [0.7,0.95,1.3]:
            xy = peak_positions(sub,center,smoothing)
            cross.append({"foreground": foreground, "smoothing_pix": smoothing,
                          "b_arcsec": float(fit(xy,center)[0][0])})
    result = {"object":ID, "input_sha256":hashes, "pixel_origin": "0-based x,y in raw[125:175,125:175]",
              "lens_center_pix":center.tolist(), "image_peaks_pix":pos.tolist(),
              "b_sis_arcsec":p[0],"gamma1":p[1],"gamma2":p[2],"source_xy_arcsec":p[3:].tolist(),
              "effective_critical_area_radius_arcsec":p[0]*np.sqrt(1+(p[1]**2+p[2]**2)/2)/(1-p[1]**2-p[2]**2),
              "source_plane_rms_arcsec":float(np.sqrt(np.mean(np.sum((beta-p[3:])**2,axis=1)))),
              "signed_inverse_magnifications":det.tolist(),"linearized_image_offsets_arcsec":delta.tolist(),
              "image_offset_rms_arcsec":float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),
              "centroid_sensitivity": {"peak_sigma_arcsec":sigma_peak,"center_sigma_arcsec":sigma_c,
                                       "b_p16_p50_p84":np.percentile(sim[:,0],[16,50,84]).tolist()},
              "subtraction_sensitivity":cross,"arc_aperture_photometry":colours(im,model,pos,center)}
    result["arc_colour_consistency"] = compare_colours(result["arc_aperture_photometry"])
    output.joinpath("refit_control.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    fig,axes=plt.subplots(1,3,figsize=(9.4,3.65),layout="constrained")
    for ax,data0,title in zip(axes[:2],[data,residual],["Observed VIS","Radial foreground removed"]):
        ax.imshow(data0,origin="lower",cmap="gray",vmin=-0.005,vmax=0.040)
        ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
    for i,(x,y) in enumerate(pos): axes[1].text(x+2,y+1,str(i+1),color="#54f3ea",fontsize=10)
    t=np.linspace(0,2*np.pi,360)
    r=p[0]*(1-p[1]*np.cos(2*t)-p[2]*np.sin(2*t))/(1-p[1]**2-p[2]**2)
    axes[1].plot(center[0]+r*np.cos(t)/PIXSCALE,center[1]+r*np.sin(t)/PIXSCALE,
                 color="#f8bb47",alpha=.9,lw=.9)
    axes[1].plot(*center,"+",color="#ff763d")
    ax=axes[2];ax.set_xlim(p[3]-0.031,p[3]+0.031);ax.set_ylim(p[4]-0.031,p[4]+0.031)
    ax.scatter(*beta.T,color="#187383")
    for i,(x,y) in enumerate(beta):ax.annotate(str(i+1),(x,y),xytext=(3,3),textcoords="offset points")
    ax.scatter(*p[3:],marker="+",color="black",s=70)
    ax.set_aspect("equal");ax.set_title("Back-projected peaks");ax.set_xlabel("Source x (arcsec)");ax.set_ylabel("Source y (arcsec)")
    axes[0].plot([3,13],[3,3],color="white",lw=2);axes[0].text(3,4,"1″",color="white")
    fig.savefig(output/"refit_control.png",dpi=150);plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args()
    main(args.cache,args.output)
