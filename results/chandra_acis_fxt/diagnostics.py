"""Generate blinded-selection follow-up cutouts from original level-2 rows."""
import csv
import gzip
import hashlib
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from astropy.io import fits

from archive import dump
from search import _gti_map, _orientation


SOURCE = "#0856a0"
BACKGROUND = "#a74500"
PEAK = "#efc15f"
TRACE = "#087c74"
DTYPE = [(name, "f8") for name in ("t", "x", "y", "e")] + [
    ("row", "i8"), ("cx", "i4"), ("cy", "i4"), ("grade", "i4")
]


def _get_data(path, trigger):
    ccd = int(trigger["ccd_id"])
    start, stop = trigger["start_met_s"], trigger["stop_met_s"]
    with fits.open(path, memmap=False) as hdul:
        evt = hdul["EVENTS"]
        if int(evt.header["OBS_ID"]) != int(trigger["obsid"]):
            raise ValueError("Diagnostic input and trigger identify different observations")
        matches = [(a,b) for a,b in _gti_map(hdul)[ccd] if a <= start < stop <= b]
        if len(matches) != 1:
            raise ValueError("Trigger does not lie within exactly one CCD good-time interval")
        lo, hi = matches[0]
        if not np.allclose([start-lo, hi-stop], [trigger["before_s"], trigger["after_s"]],
                           atol=1e-6, rtol=0):
            raise ValueError("Trigger baseline differs from input GTI")
        names = {name.upper():name for name in evt.columns.names}
        chunks = []
        for i in range(0,len(evt.data),150000):
            block = evt.data[i:i+150000]
            col = lambda k: block[names[k]]
            status = col("STATUS")
            mask = (col("CCD_ID")==ccd) & (col("TIME")>=lo) & (col("TIME")<hi)
            mask &= (col("ENERGY")>=500) & (col("ENERGY")<=7000)
            mask &= (status == 0).reshape(len(block),-1).all(axis=1)
            mask &= np.isin(col("GRADE"), [0,2,3,4,6]) & np.isfinite(col("X")) & np.isfinite(col("Y"))
            picks = np.flatnonzero(mask)
            table = np.empty(len(picks),dtype=DTYPE)
            table["row"] = i+picks
            for key, field in (("t","TIME"),("x","X"),("y","Y"),("e","ENERGY"),("cx","CHIPX"),
                               ("cy","CHIPY"),("grade","GRADE")):
                table[key] = col(field)[picks]
            chunks.append(table)
    if not chunks:
        raise ValueError("No events for the diagnostic image")
    return np.concatenate(chunks), float(lo), float(hi)


def _events_csv(path, data, aperture, annulus, start, stop, p0, p1):
    with path.open("wb") as file:
        with gzip.GzipFile(filename="",mode="wb",fileobj=file,mtime=0) as gz:
            with io.TextIOWrapper(gz,encoding="ascii",newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["events_row_zero_based","region","phase","in_peak","time_met_s",
                                 "x","y","chipx","chipy","energy_eV","grade"])
                for idx in np.flatnonzero(aperture|annulus):
                    event = data[idx]
                    t = float(event["t"])
                    region = "aperture" if aperture[idx] else "annulus"
                    phase = "before" if t < start else "trigger" if t < stop else "after"
                    writer.writerow([int(event["row"]),region,phase,int(p0<=t<p1),t,
                                     float(event["x"]),float(event["y"]),int(event["cx"]),
                                     int(event["cy"]),float(event["e"]),int(event["grade"])])


def _density(ax,x,y,bins,limits):
    image, xs, ys = np.histogram2d(x,y,bins=bins,range=limits)
    vals = np.log1p(image)
    high = max(1., np.quantile(vals[image>0],.995)) if image.any() else 1.
    ax.imshow(vals.T,extent=[xs[0],xs[-1],ys[0],ys[-1]],origin="lower",cmap="Greys",
              vmin=0,vmax=high,interpolation="nearest",aspect="equal")


def _emit_fig(figure_path, curves_path, data, trig, gti, aperture, annulus, far):
    t0, t1 = trig["start_met_s"], trig["stop_met_s"]
    w = t1-t0
    p0 = trig["peak_start_met_s"]-t0
    p1 = p0+trig["peak_window_s"]
    lo, hi = np.array(gti)-t0
    r, x, y = (trig[k] for k in ("aperture_radius_px","aperture_x","aperture_y"))
    t = data["t"]-t0
    on, peak = (t>=0)&(t<w), (t>=p0)&(t<p1)
    sx, sy = data["x"]-x, data["y"]-y
    matrix = _orientation(data)
    fig, axs = plt.subplots(3,3,figsize=(14,11),layout="constrained")
    try:
        fig.suptitle(f"Chandra ObsID {trig['obsid']}  •  CCD {trig['ccd_id']}  •  "
                     f"ICRS {trig['ra_deg']:.6f}°, {trig['dec_deg']:+.6f}°\n"
                     f"0.5–7 keV  |  trigger starts at mission time {t0:.3f} s  |  "
                     f"window {w:.0f} s  |  amber time band marks the 60 s peak",size=14)
        field, sky, chip = axs[0]
        _density(field,sx,sy,(220,220),[[sx.min()-1,sx.max()+1],[sy.min()-1,sy.max()+1]])
        if matrix is not None:
            axis = matrix[1]/np.linalg.norm(matrix[1])
            line = np.array([-800.,800.])
            xlim,ylim=field.get_xlim(),field.get_ylim()
            field.plot(line*axis[0],line*axis[1],lw=.8,color=TRACE,ls=":",label="readout axis")
            field.set(xlim=xlim,ylim=ylim)
        field.add_patch(Circle((0,0),4*r,fc="none",ec=SOURCE,lw=1.3))
        if matrix is not None:
            field.legend(fontsize=7,loc="upper right")
        field.set(title="CCD field over the full good-time interval",xlabel="Sky X offset (pixel)",ylabel="Sky Y offset (pixel)")
        _density(sky,sx,sy,(80,80),[[-5*r,5*r],[-5*r,5*r]])
        close_on = on & (np.abs(sx)<=5*r) & (np.abs(sy)<=5*r)
        sky.scatter(sx[close_on],sy[close_on],s=7,c=SOURCE,alpha=.7)
        for radius,color in [(r,SOURCE),(2*r,BACKGROUND),(4*r,BACKGROUND)]:
            sky.add_patch(Circle((0,0),radius,fc="none",ec=color,lw=.9))
        sky.set(xlim=(-5*r,5*r),ylim=(-5*r,5*r),title="Trigger photons and extraction apertures",
                xlabel="Sky X offset (pixel)",ylabel="Sky Y offset (pixel)")
        src_on = aperture & on
        cx0,cy0 = (float(np.median(data[k][src_on])) for k in ("cx","cy"))
        d = max(30,4*r)
        _density(chip,data["cx"],data["cy"],(90,90),[[cx0-d,cx0+d],[cy0-d,cy0+d]])
        for mask,color,label in ((src_on & ~peak,SOURCE,"trigger"),(aperture & peak,BACKGROUND,"peak")):
            chip.scatter(data["cx"][mask],data["cy"][mask],s=12,c=color,alpha=.7,label=label)
        chip.legend(fontsize=7)
        chip.set(xlim=(cx0-d,cx0+d),ylim=(cy0-d,cy0+d),title="Aperture events on the detector",
                 xlabel="CHIPX column (pixel)",ylabel="CHIPY row (pixel)")
        for ax in axs[1:].ravel():
            ax.axvspan(0,w,color="#d2e4f2",alpha=.6,zorder=0)
            ax.axvspan(p0,p1,color=PEAK,alpha=.45,zorder=0)
            ax.set_xlabel("Time from trigger start (s)")
        zoom_edges = np.linspace(-w,2*w,37)
        full_edges = np.linspace(lo,hi,101)
        with curves_path.open("w",encoding="ascii",newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["series","start_met_s","stop_met_s","live_s","aperture_counts",
                             "annulus_counts","ccd_counts_outside_4r"])
            for label,edges,ax in [("full_gti",full_edges,axs[1,0]),("local",zoom_edges,axs[1,1])]:
                # Exclude the rightmost boundary to match FITS half-open selection.
                counts = [np.histogram(t[mask & (t<edges[-1])],bins=edges)[0]
                          for mask in (aperture,annulus,far)]
                live = np.diff(edges)*trig["live_fraction"]
                for j in range(len(live)):
                    writer.writerow([label,edges[j]+t0,edges[j+1]+t0,live[j],
                                     *[int(c[j]) for c in counts]])
                ax.stairs(counts[0]/live,edges,color=SOURCE,label="aperture")
                ax.stairs(counts[1]/(trig["annulus_area_ratio"]*live),edges,color=BACKGROUND,
                          label="scaled annulus")
                ax.set(xlim=(edges[0],edges[-1]),ylabel="Counts / live second")
                ax.legend(fontsize=7)
                if label=="local":
                    ccd_ax=axs[1,2]
                    ccd_ax.stairs(counts[2]/live,edges,color=TRACE)
                    ccd_ax.set(xlim=(edges[0],edges[-1]),ylabel="Counts / live second")
        axs[1,0].set_title("Source and background over the full GTI")
        axs[1,1].set_title("Source and background near the flare")
        axs[1,2].set_title("CCD background outside the outer annulus")
        local=(t>=-w)&(t<2*w)
        energy=axs[2,0]
        for mask,color,label in ((annulus,BACKGROUND,"annulus"),(aperture,SOURCE,"aperture")):
            selected=mask & local
            energy.scatter(t[selected],data["e"][selected]/1000,s=10,c=color,alpha=.65,label=label)
        energy.set(xlim=(-w,2*w),ylim=(.4,7.1),ylabel="Energy (keV)",title="Energy and arrival time")
        energy.legend(fontsize=7)
        for i,field in enumerate(("cx","cy")):
            ax=axs[2,i+1]
            ax.scatter(t[aperture & local],data[field][aperture & local],s=10,c=SOURCE,label="aperture")
            independent = far & local
            if matrix is not None and independent.any():
                # Predict the detector row/column of the fixed aperture center
                # from independent photons. Bin medians trace telescope dither.
                inv = np.linalg.inv(matrix)
                pred = data[field][independent] - sx[independent]*inv[0,i] - sy[independent]*inv[1,i]
                times = t[independent]
                groups = np.digitize(times,zoom_edges)-1
                med = [np.median(pred[groups==j]) if (groups==j).any() else np.nan
                       for j in range(len(zoom_edges)-1)]
                ax.plot((zoom_edges[:-1]+zoom_edges[1:])/2,med,color=TRACE,marker=".",lw=1,
                        label="field-derived center")
            ax.set(xlim=(-w,2*w),title=("Detector column versus time","Detector row versus time")[i],
                   ylabel=("CHIPX (pixel)","CHIPY (pixel)")[i])
            ax.legend(fontsize=7)
        fig.savefig(figure_path,dpi=140)
    finally:
        plt.close(fig)


def save_diagnostics(event_path, trigger, output_dir):
    """Write a diagnostic figure, original EVENTS row index, and light curves."""
    event_path, out = Path(event_path),Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    data, lo, hi = _get_data(event_path,trigger)
    start,stop=(trigger[k] for k in ("start_met_s","stop_met_s"))
    p0=trigger["peak_start_met_s"]
    p1=p0+trigger["peak_window_s"]
    t=data["t"]
    on,peak=(t>=start)&(t<stop),(t>=p0)&(t<p1)
    dr=np.hypot(data["x"]-trigger["aperture_x"],data["y"]-trigger["aperture_y"])
    r=trigger["aperture_radius_px"]
    source,annulus,far=dr<=r,(dr>=2*r)&(dr<=4*r),dr>4*r
    masks = {"on_counts":source&on,"before_counts":source&(t<start),"after_counts":source&(t>=stop),
             "annulus_counts":annulus&on,"peak_counts":source&peak,"peak_annulus_counts":annulus&peak,
             "ccd_on_counts":far&on,"ccd_off_counts":far&~on}
    counts = {key:int(mask.sum()) for key,mask in masks.items()}
    if any(counts[k] != trigger[k] for k in counts):
        raise ValueError(f"Diagnostic selection disagrees with trigger counts: {counts}")
    rows = out/"event_rows.csv.gz"
    fig = out/"diagnostics.png"
    curves = out/"lightcurves.csv"
    _events_csv(rows,data,source,annulus,start,stop,p0,p1)
    _emit_fig(fig,curves,data,trigger,(lo,hi),source,annulus,far)
    def file_hash(path):
        with path.open("rb") as stream:
            return hashlib.file_digest(stream,"sha256").hexdigest()
    receipt = dict(input_name=event_path.name,input_sha256=file_hash(event_path),
                   ccd_id=int(trigger["ccd_id"]),gti_met_s=[lo,hi],counts=counts,
                   apertures_sky_px=[r,2*r,4*r],sky_center_px=[trigger["aperture_x"],trigger["aperture_y"]],
                   event_table_rows=int(np.count_nonzero(source|annulus)),
                   files={p.name:file_hash(p) for p in (rows,fig,curves)},
                   interpretation="Images show log-scaled counts per spatial bin over the full CCD GTI. "
                                  "Colored markers identify selected trigger photons. Detector traces "
                                  "use independent field photons and a fitted chip-to-sky rotation. "
                                  "GTI columns and exposure-normalized curves use the triggered CCD only.")
    dump(out/"diagnostics.json",receipt)
    return receipt
