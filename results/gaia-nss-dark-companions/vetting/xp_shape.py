"""Compare reddened optical continuum shapes of HD 264291 with fixed companion scenarios."""

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "epochs"))
from catalog_io import read_header, table_dtype


SOURCE = "3378588057203660160"
XP_URL = "https://gea.esac.esa.int/data-server/data?" + urllib.parse.urlencode(
    dict(
        RETRIEVAL_TYPE="XP_SAMPLED",
        ID=SOURCE,
        FORMAT="CSV",
        DATA_STRUCTURE="DATAMODEL_STANDARD",
        RELEASE="Gaia DR3",
    )
)
CK_URL = "https://ssb.stsci.edu/trds/grid/ck04models/ckp00/"
PRIMARY_TEMPS = np.arange(7000, 11251, 250)


def acquire(url, path):
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            raw = response.read()
        path.write_bytes(raw)
    raw = path.read_bytes()
    return raw, dict(url=url, sha256=hashlib.sha256(raw).hexdigest(), nbytes=len(raw))


def ccm89(wave, r_v):
    """CCM 1989 A(lambda)/A(V), optical and NIR polynomial branches, wavelength in nm."""
    x = 1000 / wave
    y = x - 1.82
    a = np.polynomial.polynomial.polyval(
        y, [1, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999]
    )
    b = np.polynomial.polynomial.polyval(
        y, [0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002]
    )
    return np.where(x < 1.1, (0.574 - 0.527 / r_v) * x**1.61, a + b / r_v)


def main(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    xp, xp_record = acquire(XP_URL, cache / "xp.csv")
    obs = np.genfromtxt(
        io.BytesIO(xp), delimiter=",", names=True, dtype=None, encoding="ascii"
    )
    assert set(obs["source_id"]) == {int(SOURCE)}
    wave, measured = obs["wavelength"], obs["flux"]

    # Broad optical bins reduce sensitivity to imperfect line spectra and XP resolution.
    edges = np.arange(370, 971, 25)
    selections = [(wave >= lo) & (wave < hi) for lo, hi in zip(edges[:-1], edges[1:])]
    centers = (edges[1:] + edges[:-1]) / 2
    y = np.array([np.mean(measured[sel]) for sel in selections])
    temps = sorted(set(PRIMARY_TEMPS.tolist()) | {5250, 5500, 5750})

    def model(temp):
        name = f"ckp00_{temp}.fits"
        raw, record = acquire(CK_URL + name, cache / name)
        f = io.BytesIO(raw)
        read_header(f)
        h = read_header(f)
        table = np.frombuffer(
            raw, dtype=table_dtype(h), count=int(h["NAXIS2"]), offset=f.tell()
        )
        model_nm = table["WAVELENGTH"].astype(float) / 10
        grid = {g: np.interp(wave, model_nm, table[g]) for g in ("g35", "g40", "g45")}
        return temp, grid, record

    with ThreadPoolExecutor(6) as pool:
        all_models = list(pool.map(model, temps))
    models = {temp: grid for temp, grid, _ in all_models}
    inputs = [xp_record] + [record for _, _, record in all_models]

    # 1.8 Rsun primary and the comparison radii only fix the optical flux ratio.
    # Absolute angular size is profiled out in all models.
    scenarios = {
        "single": [],
        "plus_two_5500K_0p9R": [(5500, "g45", 2 * (0.90 / 1.8) ** 2)],
        "plus_one_7750K_1p65R": [(7750, "g40", (1.65 / 1.8) ** 2)],
        "plus_two_5250K_0p85R": [(5250, "g45", 2 * (0.85 / 1.8) ** 2)],
        "plus_two_5750K_0p95R": [(5750, "g45", 2 * (0.95 / 1.8) ** 2)],
    }
    curves = {}
    best = {}
    for label, others in scenarios.items():
        best_score = np.inf
        for t in PRIMARY_TEMPS:
            for logg in ("g35", "g40", "g45"):
                primary = models[t][logg]
                companion = sum(area * models[t2][g2] for t2, g2, area in others)
                spectrum = primary + companion
                for rv in (2.5, 3.1, 4.0):
                    # Project extinction at each XP wavelength before binning.
                    k = ccm89(wave, rv)

                    def residual(av):
                        reddened = spectrum * 10 ** (-0.4 * av * k)
                        bands = np.array([np.mean(reddened[sel]) for sel in selections])
                        dv = np.log(bands / y)
                        return dv - np.mean(dv)

                    def score(av):
                        return float(np.mean(residual(av) ** 2))

                    opt = minimize_scalar(
                        score, method="bounded", bounds=(0, 2), options={"xatol": 0.001}
                    )
                    if opt.fun >= best_score:
                        continue
                    best_score = opt.fun
                    resid = residual(opt.x)
                    fraction550 = (
                        np.interp(550, wave, companion / spectrum) if others else 0.0
                    )
                    best[label] = dict(
                        temp_primary_K=int(t),
                        logg_primary=float(logg[1:]) / 10,
                        Av=float(opt.x),
                        Rv=rv,
                        companion_flux_fraction_550nm=float(fraction550),
                        rms_logflux_percent=100 * np.sqrt(opt.fun),
                        max_abs_fractional_residual_percent=float(
                            100 * np.max(abs(np.expm1(resid)))
                        ),
                    )
                    curves[label] = y * np.exp(resid)
    with (out / "xp_shape_bands.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(["bin_center_nm", "xp_mean_W_m2_nm", *scenarios])
        writer.writerows(zip(centers, y, *[curves[s] for s in scenarios]))
    (out / "xp_shape.json").write_text(
        json.dumps(
            dict(
                source_id=SOURCE,
                inputs=inputs,
                n_xp_samples=len(wave),
                n_bands=len(centers),
                bin_width_nm=25,
                scenarios=best,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    main(a.cache, a.out)
