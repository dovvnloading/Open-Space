"""Rebuild the MAXI/TESS flare measurements. Requires numpy and pandas.

    python analyze.py --out rebuilt

The complete light curves are retrieved from the mission archive. Only a short
byte range of each much larger target-pixel file is needed.
"""
import argparse
import hashlib
import io
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

EVENTS = [
    ("AT Mic", "2018-08-16T19:22:36", "8346898070", "J2041-324",
     "tess2018206045859-s0001-0000000441398770-0120-s_lc.fits",
     "d541fa0cf2bca7520b461ecfca4127463846aa1a4a16ef0ce9f6fbefdf6276fb"),
    ("UX Ari", "2021-10-22T21:52:51", "9509734289", "J0326+287",
     "tess2021284114741-s0044-0000000029701408-0215-s_lc.fits",
     "f00c11fb9f90d587e4f4919489c93d1bd4e7c842e1e8716cb1784fd096d20867"),
]

def digest(data):
    return hashlib.sha256(data).hexdigest()

def bytes_at(url, start=None, stop=None):
    headers = {} if start is None else {"Range": f"bytes={start}-{stop}"}
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=100) as response:
                if start is None:
                    return response.read()
                if response.status == 206:
                    if not response.headers["Content-Range"].startswith(f"bytes {start}-"):
                        raise ValueError("Unexpected HTTP byte range")
                else:
                    for position in range(0, start, 65536):
                        response.read(min(65536, start - position))
                data = response.read(stop - start + 1)
                if len(data) != stop - start + 1:
                    raise ValueError("Incomplete archive byte range")
                return data
        except OSError:
            if attempt == 2:
                raise

def header(buf, pos):
    values = {}
    while True:
        card = buf[pos:pos + 80].decode("ascii")
        pos += 80
        if len(card) != 80:
            raise ValueError("Incomplete FITS header")
        if card.startswith("END     "):
            return values, (pos + 2879) // 2880 * 2880
        if card[8:10] == "= ":
            values[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()

def table_layout(data):
    primary, offset = header(data, 0)
    ext, offset = header(data, offset)
    numeric = {"D": ">f8", "E": ">f4", "K": ">i8", "J": ">i4", "I": ">i2", "B": "u1"}
    columns = []
    for i in range(1, int(ext["TFIELDS"]) + 1):
        count, code = re.fullmatch(r"(\d*)([A-Z])", ext[f"TFORM{i}"]).groups()
        shape = (int(count),) if count and int(count) > 1 else ()
        columns.append((ext[f"TTYPE{i}"], numeric[code], shape))
    dtype = np.dtype(columns)
    assert dtype.itemsize == int(ext["NAXIS1"])
    return ext, offset, dtype

def iso(tjd):
    return datetime.fromtimestamp((tjd + 2457000 - 2440587.5) * 86400, timezone.utc).isoformat(timespec="seconds")

def median_window(t, y, hours, minimum):
    series = pd.Series(y, index=pd.to_timedelta(t, unit="d"))
    return series.rolling(pd.Timedelta(hours=hours), center=True, min_periods=minimum).median().values

def optical(d, ref, compact, out, tag):
    good = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
    obs = d[good]
    t = obs["TIME"] - obs["TIMECORR"] - 69.184 / 86400
    hours = (t - ref) * 24
    inner, outer, integ = (1.5, 4., (-.45, 1.35)) if compact else (6., 12., (-5., 3.))
    search = (-.5, .5) if compact else (-2., 1.)
    variations = []
    for key in ["SAP_FLUX", "PDCSAP_FLUX"]:
        flux = obs[key].astype(float)
        for degree in [1, 2]:
            for scale in [1., 1.33]:
                bg = (abs(hours) >= inner * scale) & (abs(hours) <= outer * scale)
                model = np.polyfit(hours[bg], flux[bg], degree)
                excess = flux / np.polyval(model, hours) - 1
                candidates = np.where((hours > search[0]) & (hours < search[1]))[0]
                peak = candidates[np.argmax(excess[candidates])]
                use = np.where((hours >= integ[0]) & (hours <= integ[1]))[0]
                step = np.diff(t[use]) * 86400
                area = (excess[use[1:]] + excess[use[:-1]]) / 2 * step
                ed = area[step < 246].sum()
                variations.append([excess[peak], ed])
                if (key, degree, scale) == ("PDCSAP_FLUX", 1, 1.):
                    best = dict(peak_utc=iso(t[peak]), peak_excess=excess[peak],
                                seconds_peak_before_alert=(ref - t[peak]) * 86400,
                                excess_at_alert=float(np.interp(0, hours, excess)),
                                equivalent_duration_s=ed, integration_hours=integ,
                                exposure_integrated_s=step[step < 246].sum())
                    span = abs(hours) < (3 if compact else 9)
                    pd.DataFrame(dict(cadence=obs["CADENCENO"][span].astype(int),
                                      utc_tjd=t[span], sap=obs["SAP_FLUX"][span].astype(float),
                                      pdcsap=flux[span], fractional_excess=excess[span])).to_csv(
                                          out / (tag + "_tess.csv"), index=False, float_format="%.9f")
    best["peak_excess_sensitivity_range"] = [min(v[0] for v in variations), max(v[0] for v in variations)]
    best["equivalent_duration_sensitivity_s"] = [min(v[1] for v in variations), max(v[1] for v in variations)]
    flux = obs["PDCSAP_FLUX"].astype(float)
    baseline = median_window(t, flux, 12, 200)
    score = median_window(t, flux / baseline - 1, 20/60, 7)
    for dt in [-.20, .20]:
        count = np.searchsorted(t, t + dt + .04) - np.searchsorted(t, t + dt - .04)
        score[count < 8] = np.nan
    value = score[np.argmin(abs(t - ref))]
    assert np.isfinite(value)
    best.update(timing_control_score=value,
                control_at_least_as_high=int(np.sum(score >= value)),
                control_eligible=int(np.sum(np.isfinite(score))))
    return best

def target_pixels(d, ref, lc_url, compact, out, tag):
    url = lc_url.replace("_lc.fits", "_tp.fits")
    h, offset, dtype = table_layout(bytes_at(url, 0, 34559))
    hours = (d["TIME"] - d["TIMECORR"] - 69.184/86400 - ref) * 24
    near = np.where(abs(hours) < (1.25 if compact else 3.5))[0]
    lo, hi = int(near[0]), int(near[-1])
    byte_range = [offset + lo*dtype.itemsize, offset + (hi+1)*dtype.itemsize - 1]
    raw = bytes_at(url, *byte_range)
    pix = np.frombuffer(raw, dtype=dtype)
    assert np.allclose(pix["TIME"], d["TIME"][lo:hi+1], rtol=0, atol=1e-7)
    dim = tuple(map(int, h["TDIM5"].strip("()").split(",")))[::-1]
    cube = pix["FLUX"].astype(float).reshape(-1, *dim)
    hours = hours[lo:hi+1]
    good = pix["QUALITY"] == 0
    is_peak = (abs(hours + (.166 if compact else 1.215)) < (.065 if compact else .20)) & good
    pre = (hours < (-.5 if compact else -3.1)) & good
    base = np.median(cube[pre], axis=0)
    diff = np.median(cube[is_peak], axis=0) - base
    def centroid(image):
        y, x = np.indices(image.shape)
        weights = np.where(image > .025*image.max(), image, 0)
        return np.array([(x*weights).sum(), (y*weights).sum()]) / weights.sum()
    sky, flare = centroid(base), centroid(diff)
    np.savetxt(out / (tag+"_difference.csv"), diff, delimiter=",", fmt="%.4f")
    return dict(url=url, bytes_inclusive=byte_range, raw_sha256=digest(raw),
                quiescent_xy=sky.tolist(), difference_xy=flare.tolist(),
                displacement_pixels=float(np.linalg.norm(sky-flare)),
                pre_cadences=int(pre.sum()), peak_cadences=int(is_peak.sum()))

def xray(src, ref, out, tag, compact):
    url = f"https://maxi.riken.jp/pubdata/v7.7l/{src}/glcscan_lcbg_gsc0124578.csv"
    data = bytes_at(url)
    clipped = b"".join(row for row in data.splitlines(keepends=True)
                       if abs(float(row.split(b",", 1)[0]) - ref - 56999.5) < 16)
    (out / (tag+"_maxi_scans.csv")).write_bytes(clipped)
    scan = np.loadtxt(io.BytesIO(clipped), delimiter=",")
    scan = scan[(scan[:, 11] == 0) & (scan[:, 6] > 0)]
    dt = (scan[:, 0] - ref - 56999.5)*24
    j = np.argmin(abs(dt))
    assert abs(dt[j])*60 < 4
    null = abs(dt) > 12
    score = scan[:, 5] / scan[:, 6]
    result = dict(url=url, subset_sha256=digest(clipped), scan_utc=iso(scan[j,0]-56999.5),
                  exposure_s=(scan[j,2]-scan[j,1])*86400,
                  flux_2_4=scan[j,5], error_2_4=scan[j,6],
                  flux_4_10=scan[j,7], error_4_10=scan[j,8],
                  null_at_least_as_high=int(np.sum(score[null]>=score[j])),
                  null_scans=int(null.sum()))
    if not compact:
        dt_pair = np.diff(scan[:,0])*24
        pair = null[1:] & null[:-1] & (dt_pair > 1.1) & (dt_pair < 2)
        flux, err = scan[:,5], scan[:,6]
        z = (flux[:-1]/err[:-1]**2 + flux[1:]/err[1:]**2)/np.sqrt(err[:-1]**-2 + err[1:]**-2)
        result.update(following_flux_2_4=flux[j+1], following_error_2_4=err[j+1],
                      two_scan_score=z[j], pairs_at_least_as_high=int(np.sum(z[pair]>=z[j])),
                      eligible_pairs=int(pair.sum()))
    return result

def main(out):
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, dt, alert, src, filename, sha in EVENTS:
        tag = name.replace(" ", "").lower()
        compact = name == "AT Mic"
        _, sector, tic, *_ = filename.split("-")
        path = "/".join(tic[i:i+4] for i in range(0,16,4))
        url = f"https://s3.amazonaws.com/stpubdata/tess/public/tid/{sector}/{path}/{filename}"
        raw = bytes_at(url)
        assert digest(raw) == sha
        h, offset, dtype = table_layout(raw)
        d = np.frombuffer(raw, dtype=dtype, count=int(h["NAXIS2"]), offset=offset)
        ref = datetime.fromisoformat(dt).replace(tzinfo=timezone.utc).timestamp()/86400+2440587.5-2457000
        results[name] = dict(alert_utc=dt, alert_url=f"https://maxi.riken.jp/alert/novae/{alert}/{alert}.htm",
                             tess_url=url, tess_sha256=sha, optical=optical(d, ref, compact, out, tag),
                             pixel=target_pixels(d, ref, url, compact, out, tag), xray=xray(src, ref, out, tag, compact))
    (out / "measurements.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args().out)
