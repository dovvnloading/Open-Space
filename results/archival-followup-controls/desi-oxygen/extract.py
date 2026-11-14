"""Rebuild compact pixel windows from public DESI DR1 products, verifying exact input bytes."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
from astropy.io import fits
from scipy.sparse import dia_matrix

HERE = Path(__file__).resolve().parent
C = 299792.458
Z = 0.8416150156358213
LINES = {"strong": 5008.239, "weak": 4960.295}
PRODUCT_ROOT = "https://data.desi.lbl.gov/public/dr1/"
SPECTRA = PRODUCT_ROOT + "spectro/redux/iron/healpix/main/"
PRODUCTS = {
    "desi-dark.fits": (SPECTRA + "dark/135/13510/coadd-main-dark-13510.fits",
                       "78a6e975249207c00e8f91ee773b933657238e4952153249a1f32567bdb6133f"),
    "desi-bright.fits": (SPECTRA + "bright/135/13510/coadd-main-bright-13510.fits",
                         "d4593263729af05837cd9ec733567b30c76576f84eb6dfb5f5beedf094bb3e12"),
    "cframe-z6.fits.gz": (PRODUCT_ROOT + "spectro/redux/iron/exposures/20211111/00108405/cframe-z6-00108405.fits.gz",
                          "252ab4082c2431bb3ebe75e8acaa22adc87954c267aa8f4e08b7d8db70a5ff5a"),
    "fastspec.fits.gz": (PRODUCT_ROOT + "vac/dr1/fastspecfit/iron/v3.0/healpix/main/dark/135/13510/fastspec-main-dark-13510.fits.gz",
                         "601a5c869e0c8c5c20f02796514c70af5c104b3075113745c0550be3ea55a744"),
}
TARGETS = {"offset": ("desi-dark.fits", 39628534716566336),
           "central": ("desi-bright.fits", 39628534716566366)}


def verify(path, expected):
    with path.open("rb") as stream:
        actual = hashlib.file_digest(stream, "sha256").hexdigest()
    if actual != expected:
        raise ValueError(f"SHA256 mismatch for {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    for filename, (url, digest) in PRODUCTS.items():
        path = args.cache / filename
        if not path.exists():
            part = path.with_suffix(".partial")
            urlretrieve(url, part)
            verify(part, digest)
            part.replace(path)
        verify(path, digest)

    packet = {}
    metadata = {"products": PRODUCTS, "targets": {}, "redshift_reference": Z}
    for aperture, (filename, target) in TARGETS.items():
        with fits.open(args.cache / filename, memmap=False) as hdus:
            table = hdus["FIBERMAP"].data
            row = np.flatnonzero(table["TARGETID"] == target).item()
            exposures = hdus["EXP_FIBERMAP"].data
            exp = exposures[exposures["TARGETID"] == target]
            metadata["targets"][aperture] = {
                "targetid": str(target), "coadd_fiberstatus": int(table["COADD_FIBERSTATUS"][row]),
                "exposures": [{k: int(e[k]) for k in ("NIGHT", "EXPID", "TILEID", "FIBER")} for e in exp]
            }
            rd = hdus["Z_RESOLUTION"].data[row].astype(float)
            wave = hdus["Z_WAVELENGTH"].data
            offsets = rd.shape[0] // 2 - np.arange(rd.shape[0])
            response = dia_matrix((rd, offsets), shape=(len(wave), len(wave))).tocsr()
            for name, line in LINES.items():
                key = f"{aperture}_{name}_"
                select = abs(wave / (line * (1 + Z)) - 1) * C < 1100
                packet[key + "wave"] = wave[select].copy()
                for column in ("FLUX", "IVAR", "MASK"):
                    packet[key + column.lower()] = hdus["Z_" + column].data[row, select].copy()
                packet[key + "R"] = response[select, :][:, select].toarray()

    with fits.open(args.cache / "cframe-z6.fits.gz", memmap=False) as hdus:
        table = hdus["FIBERMAP"].data
        sky = table["OBJTYPE"] == "SKY"
        good = sky & (table["FIBERSTATUS"] == 0)
        metadata["sky"] = {"designated": int(sky.sum()), "zero_status": int(good.sum()),
                           "fiber_ids": table["FIBER"][good].astype(int).tolist()}
        wave = hdus["WAVELENGTH"].data
        for name in LINES:
            sample = packet["offset_" + name + "_wave"]
            idx = np.searchsorted(wave, sample)
            np.testing.assert_allclose(wave[idx], sample)
            for column in ("FLUX", "IVAR", "MASK"):
                packet[f"sky_{name}_{column.lower()}"] = hdus[column].data[good][:, idx].copy()

    with fits.open(args.cache / "fastspec.fits.gz") as hdus:
        i = np.flatnonzero(hdus["METADATA"].data["TARGETID"] == int(metadata["targets"]["offset"]["targetid"])).item()
        cols = ("OII_DOUBLET_RATIO", "OIII_5007_SIGMA", "OIII_5007_FLUX", "OIII_4959_FLUX",
                "OII_3726_FLUX", "OII_3729_FLUX", "OIII_5007_CHI2")
        metadata["fastspecfit"] = {k: float(hdus["FASTSPEC"].data[k][i]) for k in cols}
    np.savez_compressed(HERE / "windows.npz", **packet)
    (HERE / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
