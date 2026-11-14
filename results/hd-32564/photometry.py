"""Test circular conjunctions in two TESS sectors and calibrate box-depth recovery."""
import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parent
PRODUCTS = {
    5: ("tess2018319095959-s0005-0000000213078996-0125-s",
        "407d80dd1d4e3b20f4c0f771204085643c166973b70c1e849fe85ddd2509ff0f"),
    32: ("tess2020324010417-s0032-0000000213078996-0200-s",
         "ac7d1b397c980892c8f690104bbadcb2c1b79daa839d1c8f4d6781e3bbec83ea")
}
DURATION_H = (2., 3.5)


def detrend(flux):
    return (flux / median_filter(flux, size=721, mode="nearest") - 1) * 1e6


def depth(series, phase):
    samples = [y[np.abs((ph-phase+.5) % 1-.5) < width/2]
               for ph, y, width in series]
    return -np.mean(np.concatenate(samples))


def run(raw_dir):
    model = min(json.loads((ROOT/"harps_predictions.json").read_text()), key=lambda m: m["rms"])
    tests = [[], []]
    injected = [[], []]
    coverage = {}
    for sector, (name, expected_sha) in PRODUCTS.items():
        filename = name + "_lc.fits"
        url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/" + filename
        if raw_dir:
            raw = (raw_dir / filename).read_bytes()
        else:
            with urlopen(url, timeout=90) as stream:
                raw = stream.read()
        if hashlib.sha256(raw).hexdigest() != expected_sha:
            raise ValueError("Archive digest differs for " + name)
        with fits.open(io.BytesIO(raw)) as hdus:
            if hdus[0].header["TICID"] != 213078996:
                raise ValueError("Unexpected target")
            rows = hdus[1].data
            t, flux = rows["TIME"].astype(float) + 2457000., rows["PDCSAP_FLUX"].astype(float)
            ok = (rows["QUALITY"] == 0) & np.isfinite(t) & np.isfinite(flux) & (flux > 0)
            t, flux = t[ok], flux[ok]
        y = detrend(flux)
        counts = []
        for i, duration in enumerate(DURATION_H):
            p = model["P"][i]
            cosine, sine = model["params"][2*i:2*i+2]
            conjunction = .25 + np.arctan2(sine, cosine)/(2*np.pi)
            phase = ((t - 2456000.) / p - conjunction + .5) % 1 - .5
            width = duration / (24*p)
            inside = abs(phase) < width/2
            tests[i].append((phase, y, width))
            # Re-run the same flattening after injection into the archived PDC flux.
            altered = flux * (1 - .0002 * inside)
            injected[i].append((phase, detrend(altered), width))
            counts.append(int(inside.sum()))
        coverage[sector] = {"good_cadences":len(t), "box_cadences":counts,
                            "file":filename, "sha256":expected_sha, "url":url}
    results = []
    for i, series in enumerate(tests):
        shifts = np.arange(200)/200
        offsets = [p for p in shifts if min(p, 1-p) > 1.5*series[0][2]]
        controls = [depth(series, p) for p in offsets]
        d = depth(series, 0)
        results.append({"period_d":round(model["P"][i],6), "duration_h":DURATION_H[i],
                        "conjunction_depth_ppm":round(float(d),2),
                        "off_phase_05_95_percentiles_ppm":np.round(np.percentile(controls,[5,95]),2).tolist(),
                        "injected_200ppm_depth_increment_ppm":round(float(depth(injected[i],0)-d),2)})
    return dict(sectors=coverage, transit_tests=results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", type=Path, help="Use local files with the original MAST filenames")
    print(json.dumps(run(p.parse_args().raw_dir), indent=2))
