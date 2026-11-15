"""Check whether the two catalogued SDSS spectra provide independent noise."""
import hashlib
import io
import json
import urllib.request
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter


def main():
    root = Path(__file__).resolve().parent
    inputs = json.loads((root / "sources.json").read_text())
    spectra = []
    for source in inputs:
        with urllib.request.urlopen(source["url"], timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != source["sha256"]:
            raise ValueError("Input has changed")
        spectra.append(np.loadtxt(io.BytesIO(data), skiprows=2, delimiter=",").T)
    (wavelength, brighter), (other_wave, fainter) = spectra
    # The brighter-star export has per-frequency flux. Its scale cancels in r.
    brighter = brighter * 2.99792458e18 / wavelength**2
    fainter = np.interp(wavelength, other_wave, fainter)
    keep = (wavelength > 4500) & (wavelength < 5500)
    a, b = brighter[keep], fainter[keep]
    a, b = a - median_filter(a, 21), b - median_filter(b, 21)
    result = {"wavelength_A": [4500, 5500], "samples": len(a),
              "residual_r": float(np.corrcoef(a, b)[0,1]),
              "shift_control_r": {str(lag): float(np.corrcoef(a[lag:], b[:-lag])[0,1])
                                  for lag in [5, 10, 20, 40]}}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
