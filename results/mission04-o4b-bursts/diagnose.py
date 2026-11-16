"""Descriptive time-local detector checks for a ranked H/L foreground trigger.

python diagnose.py CACHE OUTPUT GPS DURATION F_LOW F_HIGH
Optional dependency: lalsuite, for reading the public 64-second auxiliary frame.
No diagnostic here defines vetoes or contributes to the false-alarm estimate.
"""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import h5py
import numpy as np
from scipy import signal

from acquire import FS, DURATION, DQ_REQUIRED, INJ_REQUIRED, fetch, require_flags
from detector_context import acquire_context, cw_lines, spectral_bands


def local_power(values, fs, frame_start, centers, duration, band, remove_bands=()):
    size = round(duration * fs)
    freq = np.fft.rfftfreq(size, 1/fs)
    bins = (freq >= band[0]) & (freq < min(band[1],fs*0.45))
    if not np.any(bins):
        return {"available":False, "reason":"no usable Fourier bins"}
    # A high-pass attenuates leakage from large low-frequency / DC levels.
    # Remove a least-squares linear trend from each prefiltered Hann tile.
    if not np.isfinite(values).all():
        return {"available":False, "reason":"nonfinite frame"}
    hi = signal.butter(4, min(24,0.2*fs), btype="highpass", fs=fs, output="sos")
    values = signal.sosfiltfilt(hi, np.asarray(values, dtype=float))
    removed_fraction = 0.0
    if remove_bands:
        # Resolve narrow persistent lines over the full local uninterrupted
        # record. A half-second end taper is inside the one-second guard.
        # Pad each interval by 0.05 Hz (>=3 frequency bins for a 64-s record).
        highres = np.fft.rfftfreq(len(values),1/fs)
        remove = np.zeros(len(highres),dtype=bool)
        for lo,hi,_ in remove_bands:
            if hi+0.05 >= band[0] and lo-0.05 < band[1]:
                left,right = np.searchsorted(highres,[lo-0.05,hi+0.05],side="left")
                remove[left:right] = True
        relevant = (highres>=band[0])&(highres<band[1])
        removed_fraction = float(np.mean(remove[relevant]))
        taper_edges = signal.windows.tukey(len(values),alpha=min(1,fs/len(values)))
        fourier = np.fft.rfft(values*taper_edges)
        fourier[remove] = 0
        values = np.fft.irfft(fourier,n=len(values))
    starts = np.rint((centers-frame_start)*fs).astype(int) - size//2
    windows = np.lib.stride_tricks.sliding_window_view(values, size)[starts]
    taper = signal.windows.hann(size,sym=False)
    power = np.abs(np.fft.rfft(signal.detrend(windows,axis=1)*taper,axis=1)[:,bins])**2
    reference = np.median(power[1:],axis=0)
    valid = np.isfinite(reference) & (reference > 0)
    if not np.any(valid):
        return {"available":False,"reason":"off-source power vanishes"}
    ratios = power[:,valid]/reference[valid]
    energy = np.mean(ratios,axis=1)
    peak = np.argmax(ratios[0])
    return {"available":True,"fraction_band_notched":removed_fraction,
            "average_to_offsource_median":float(energy[0]/np.median(energy[1:])),
            "offsource_ge":int(np.sum(energy[1:] >= energy[0])), "offsource_n":len(energy)-1,
            "peak_bin_hz":float(freq[bins][valid][peak]), "bin_width_hz":float(fs/size),
            "peak_share_of_band_power":float(ratios[0,peak]/np.sum(ratios[0])),
            "all_norm_energies":energy/np.median(energy[1:])}


def centers_near(gps, frame_start, frame_end, duration, good=None, mask_start=None):
    # Whole-second offsets preserve the target's position within a tile.
    offsets = np.arange(-64,65)
    centers = gps+offsets
    keep = (centers >= frame_start+1+duration/2)&(centers < frame_end-1-duration/2)
    if good is not None:
        for i,c in enumerate(centers):
            if not keep[i]:
                continue
            low = int(np.floor(c-duration/2)) - mask_start
            high = int(np.ceil(c+duration/2)) - mask_start
            keep[i] = np.all(good[low:high])
    if not np.any(keep & (offsets == 0)):
        return None
    return np.r_[gps, centers[keep & (np.abs(offsets)>2)]]


def strain_data(cache, block, gps, duration, band):
    out = {}
    for ifo,name in block["files"].items():
        mask_gps = block["gps"]
        with h5py.File(cache/name) as f:
            dq = require_flags(f["quality/simple"], "DQmask", "DQShortnames", DQ_REQUIRED, mask_gps)
            good = dq & require_flags(f["quality/injections"], "Injmask","InjShortnames",INJ_REQUIRED,mask_gps)
            cw_free = require_flags(f["quality/injections"], "Injmask","InjShortnames",("NO_CW_HW_INJ",),mask_gps)
            i = int(gps-mask_gps)
            if not good[i]:
                out[ifo] = {"available":False,"reason":"DQ/injection mask excludes target"}
                continue
            # Use the uninterrupted DQ-good span around target and retain 1-s
            # conditioning guards as with the 64-s auxiliary frame below.
            left, right = max(0,i-32),min(DURATION,i+33)
            bad_left = np.flatnonzero(~good[left:i]); bad_right = np.flatnonzero(~good[i:right])
            if bad_left.size:
                left += int(bad_left[-1])+1
            if bad_right.size:
                right = i+int(bad_right[0])
            start,end = mask_gps+left,mask_gps+right
            centers = centers_near(gps,start,end,duration)
            if centers is None or len(centers)<8:
                out[ifo] = {"available":False,"reason":"insufficient uninterrupted reference"}
                continue
            raw = f["strain/Strain"][left*FS:right*FS]
            info = local_power(raw,FS,start,centers,duration,band)
            info.pop("all_norm_energies",None)
            if info["available"] and ifo in ("H1","L1"):
                lines = list(spectral_bands(cache,ifo,gps))
                if not cw_free[i]:
                    lines += list(cw_lines(cache,gps))
                # Flag a catalog line whose physical interval intersects the
                # dominant Fourier resolution bin, without assigning causation.
                notched = local_power(raw,FS,start,centers,duration,band,lines)
                notched.pop("all_norm_energies",None)
                info["after_catalog_line_removal"] = notched
                peak = info["peak_bin_hz"]
                half = info["bin_width_hz"]/2
                info["dominant_bin_overlaps"] = [
                    {"low_hz":float(lo),"high_hz":float(hi),"label":label}
                    for lo,hi,label in lines if hi >= peak-half and lo < peak+half]
            out[ifo] = info
    return out


def auxiliary(cache,gps,duration,band):
    import lal
    import lalframe as lf

    url = f"https://datafind.gwosc.org/api/v1/gwf/H/H1_AUX_AR1/{int(gps)},{int(gps)+1}/osdf.json"
    with urlopen(url,timeout=30) as response:
        objects = json.load(response)
    if len(objects)!=1:
        return {"available":False,"reason":"no unique public H1 auxiliary frame","query":url}
    object_url = objects[0]
    path = cache / object_url.rsplit("/",1)[1]
    # The GWOSC OSDF namespace publishes its authoritative direct-read origin.
    fetch("https://origin.ligo.caltech.edu:8443/"+object_url.removeprefix("osdf:///"),path)
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream,"sha256").hexdigest()
    frame = lf.FrFileOpenURL(str(path))
    reader = lf.FrameUFrFileOpen(str(path),"r")
    toc = lf.FrameUFrTOCRead(reader)
    out, scores = {}, []
    for i in range(lf.FrameUFrTOCQueryProcN(toc)):
        channel = lf.FrameUFrTOCQueryProcName(toc,i)
        typ = lf.FrFileQueryChanType(frame,channel,0)
        loaders = {lal.S_TYPE_CODE: lf.FrFileReadREAL4TimeSeries, lal.D_TYPE_CODE: lf.FrFileReadREAL8TimeSeries}
        if typ not in loaders:
            out[channel] = {"available":False,"reason":"unsupported channel sample type"}
            continue
        series = loaders[typ](frame,channel,0)
        fs = round(1/series.deltaT)
        start = float(series.epoch)
        centers = centers_near(gps,start,start+len(series.data.data)/fs,duration)
        if centers is None or len(centers)<8:
            out[channel] = {"available":False,"reason":"edge or insufficient reference"}
            continue
        info = local_power(series.data.data,fs,start,centers,duration,band)
        if info["available"]:
            scores.append(info.pop("all_norm_energies"))
        out[channel] = info
    result = {"query":url,"object_url":object_url,"sha256":digest,"channels":out}
    if scores:
        maxima = np.max(np.vstack(scores),axis=0)
        result["max_over_channels"] = {"on_source":float(maxima[0]),"offsource_ge":int(np.sum(maxima[1:]>=maxima[0])),
                                      "offsource_n":len(maxima)-1,"channels":len(scores)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("--selection",type=Path, help="Use frozen holdout manifest and fetch the target block")
    parser.add_argument("output",type=Path)
    for label in ("gps","duration","low","high"):
        parser.add_argument(label,type=float)
    args = parser.parse_args()
    acquire_context(args.cache)
    selection = args.selection or args.cache/"search_selection.json"
    b = next(b for b in json.loads(selection.read_text())["blocks"]
             if b["gps"] <= args.gps < b["gps"]+DURATION)
    if args.selection is not None:
        from holdout_data import acquire_holdout
        acquire_holdout(args.cache,b)
        if "virgo_source" in b:
            rec = b["virgo_source"]
            fetch(rec["url"],args.cache/rec["name"],rec["md5"])
            b["files"]["V1"] = rec["name"]
    band = (args.low,args.high)
    result = {"gps":args.gps, "duration":args.duration, "band_hz":band,
              "strain":strain_data(args.cache,b,args.gps,args.duration,band),
              "auxiliary":auxiliary(args.cache,args.gps,args.duration,band)}
    args.output.write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
