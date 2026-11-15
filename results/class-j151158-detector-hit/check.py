"""Compare the raw morphology of a CLASS line with arc lines in the same SDSS fiber."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from numpy.polynomial.legendre import legval
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent
FIBER = 565
ROW = (FIBER - 1) % 320
RAW = "https://data.sdss.org/sas/dr17/sdss/spectro/data/53091/"
REDUX = "https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/1385/"
INPUTS = {"rr": RAW + "sdR-r2-00025898.fit.gz",
          "arc": RAW + "sdR-r2-00025900.fit.gz",
          "flat": REDUX + "spFlat-r2-00025899.fits.gz",
          "cfr-full": REDUX + "spCFrame-r2-00025898.fits"}


def raw_image(path):
    # Older SDSS writes SIMPLE=F and unsigned CCD counts without a BZERO card.
    with gzip.open(path, "rb") as handle:
        data = handle.read()
    end = next(i for i in range(0, len(data), 80) if data[i:i+8] == b"END     ")
    offset = ((end + 80 + 2879) // 2880) * 2880
    header = fits.Header.fromstring(data[:end+80].decode("ascii"))
    width, height = header["NAXIS1"], header["NAXIS2"]
    return np.frombuffer(data, dtype=">u2", count=width*height, offset=offset).reshape(height, width)


def stamp(raw, x, y):
    # opConfig-50000.par for the SDSS cameras: 40 prescan columns, 21 rows.
    cut = raw[y+21-9:y+21+10, x+40-4:x+40+5].astype(float)
    background = np.median(cut[[0, 1, 2, -3, -2, -1]], axis=0)
    return cut[5:14] - background


def shape(image):
    # Report positive residual fractions within 5x5, centered on local peak.
    iy, ix = np.unravel_index(np.argmax(image[2:7, 2:7]), (5, 5))
    core = np.maximum(image[iy:iy+5, ix:ix+5], 0)
    total = core.sum()
    return float(core.max()/total), float(np.sort(core.ravel())[-2:].sum()/total)


def main(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    expected = json.loads((ROOT / "detector.json").read_text())["inputs"]
    for name, url in INPUTS.items():
        path = cache / name
        if not path.exists():
            with urlopen(url, timeout=120) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != expected[name]["sha256"]:
                raise ValueError(f"Archive bytes changed for {url}")
            path.write_bytes(content)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected[name]["sha256"]:
            raise ValueError(f"Unexpected cached bytes for {name}")
    raw = raw_image(cache/"rr")
    arc = raw_image(cache/"arc")
    with fits.open(cache/"flat") as hdus:
        tr = hdus[1].data[0]
        assert tr["FUNC"].strip().lower() == "legendre"
        trace = legval(2*(np.arange(2048)-tr["XMIN"])/(tr["XMAX"]-tr["XMIN"])-1,
                       tr["COEFF"][ROW])
    with fits.open(cache/"cfr-full") as hdus:
        wave = 10.0 ** hdus[3].data[ROW].astype(float)
        # This wavelength is the brightest coadded pixel at the cataloged line.
        pos = int(np.argmin(abs(wave-6446.1)))
    x = np.rint(trace).astype(int)
    target = stamp(raw, x[pos], pos)
    # Trace-summed arc spectrum from the same fiber and readout geometry.
    y = np.arange(2048)
    flux = sum(arc[y+21, x+40+dx].astype(float) for dx in (-1, 0, 1))
    peaks = find_peaks(flux, prominence=650, distance=6)[0]
    peaks = [p for p in peaks if abs(p-pos) <= 250 and 10 <= p < 2038]
    rows = []
    valid = []
    for p in peaks:
        patch = stamp(arc, x[p], p)
        if np.max(arc[p+21-5:p+21+6, x[p]+40-5:x[p]+40+6]) >= 60000:
            continue
        # Exclude blended lamp peaks or an extracted neighboring fiber.
        cy,cx = np.unravel_index(patch[2:7,2:7].argmax(), (5,5))
        if abs(cx-2)>1 or abs(cy-2)>1:
            continue
        valid.append((p, patch))
        rows.append((p, wave[p], *shape(patch)))
    with (out/"lamp-controls.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["trace_row", "wavelength_vac_A", "peak_pixel_fraction", "top_two_fraction"])
        writer.writerows(rows)
    assert valid
    nearest = min(valid, key=lambda t:abs(t[0]-pos))
    vmax = 1
    fig, axes = plt.subplots(1, 2, figsize=(7.3,3.4), constrained_layout=True)
    for ax, patch, title in zip(axes, [target, nearest[1]],
                                [f"Galaxy exposure, {wave[pos]:.1f} Å",
                                 f"Arc lamp, {wave[nearest[0]]:.1f} Å"]):
        im=ax.imshow(patch/patch.max(),vmin=0,vmax=vmax,origin="lower",cmap="magma",
                     extent=(-4.5,4.5,-4.5,4.5),interpolation="nearest")
        ax.set(xlabel="Across fiber (pixels)",ylabel="Dispersion (pixels)",title=title)
    fig.colorbar(im, ax=axes, label="Fraction of peak signal",shrink=.86)
    fig.savefig(out/"detector.png", dpi=150)
    yc,xc=np.unravel_index(np.argmax(target[2:7,2:7]), (5,5))
    result=dict(fiber=FIBER, extracted_detector_row=pos, raw_predicted_xy0=[int(x[pos]+40),pos+21],
                raw_peak_xy0=[int(x[pos]+40+xc-2),int(pos+21+yc-2)],
                vacuum_wavelength_A=float(wave[pos]),target_shape=shape(target),
                controls=len(rows), largest_control_fraction=np.max(np.asarray(rows)[:,2:],axis=0).tolist(),
                target_patch_adu=target.tolist(),
                inputs={k:{"url":v,"sha256":hashlib.sha256((cache/k).read_bytes()).hexdigest()} for k,v in INPUTS.items()})
    (out/"detector.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["target_shape"], result["controls"], result["largest_control_fraction"])


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=ROOT/"cache")
    parser.add_argument("--out",type=Path,default=ROOT)
    args=parser.parse_args()
    main(args.cache, args.out)
