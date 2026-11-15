"""Replay Earth-view closest approaches using independent infrared positions."""
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

HERE = Path(__file__).resolve().parent
ORIGIN = np.array([300.209094, 36.497796])
EPOCH = Time(57255.66, format="mjd", scale="tdb").jyear
SPITZER = np.array([0.0, 0.0, 6.1, 372.8, 133.4])


def fit(epochs, reference, seasonal):
    """Equal weight per sky visit. Seasonal fit adopts Spitzer parallax."""
    registered = epochs.join(reference.groupby("epoch")[["dx", "dy"]].median(), on="epoch", validate="many_to_one")
    x = registered.x - registered.dx
    y = registered.y - registered.dy
    dt = epochs.year - EPOCH
    one, zero = np.ones(len(dt)), np.zeros(len(dt))
    if seasonal:
        sign = np.sign(epochs.fx)
        east = np.stack([one, zero, dt, zero, sign, zero], axis=1)
        north = np.stack([zero, one, zero, dt, zero, sign], axis=1)
        x = x - SPITZER[4] * epochs.fx
        y = y - SPITZER[4] * epochs.fy
    else:
        east = np.stack([one, zero, dt, zero, epochs.fx], axis=1)
        north = np.stack([zero, one, zero, dt, epochs.fy], axis=1)
    design = np.vstack([east, north])
    data = np.r_[x, y]
    coefficients = np.linalg.lstsq(design, data, rcond=None)[0]
    residual = data - design @ coefficients
    error = np.sqrt(np.diag(np.linalg.inv(design.T @ design)) * (residual @ residual)
                    / (len(data) - len(coefficients)))
    track = np.r_[coefficients[:4], SPITZER[4]] if seasonal else coefficients
    return track, error


def main():
    sources = pd.read_csv(HERE / "sources.csv", dtype={"source_id": str})
    refs = pd.read_csv(HERE / "reference-offsets.csv", dtype={"source_id": str})
    models = {"published Spitzer": SPITZER}
    for product in ["W2", "L1b"]:
        epochs = pd.read_csv(HERE / f"{product}-epochs.csv")
        for seasonal in [False, True]:
            name = product + (" fixed parallax + scan offsets" if seasonal else " free parallax")
            track, error = fit(epochs, refs[refs["product"] == product], seasonal)
            models[name] = track
            print(name, "offsets, motion, parallax", np.round(track, 2), "formal errors", np.round(error, 2))
    years = np.arange(2032.0, 2042.0, .002)
    time = Time(years, format="jyear", scale="tdb")
    X, Y, Z = get_body_barycentric("earth", time).xyz.value
    ra, dec = np.deg2rad(ORIGIN)
    east = np.sin(ra)*X - np.cos(ra)*Y
    north = np.cos(ra)*np.sin(dec)*X + np.sin(ra)*np.sin(dec)*Y - np.cos(dec)*Z
    for name, lens in models.items():
        for s in sources.itertuples(index=False):
            relative_plx = lens[4] - max(0., s.parallax)
            dx = ((s.ra - ORIGIN[0])*3.6e6*np.cos(dec) + s.pmra*(years-2016)
                  - lens[0] - lens[2]*(years-EPOCH) - relative_plx*east)
            dy = ((s.dec - ORIGIN[1])*3.6e6 + s.pmdec*(years-2016)
                  - lens[1] - lens[3]*(years-EPOCH) - relative_plx*north)
            separation = np.hypot(dx, dy)
            i = np.argmin(separation)
            print(name, s.source_id, time[i].iso[:10], round(separation[i], 1), "mas")


if __name__ == "__main__":
    main()
