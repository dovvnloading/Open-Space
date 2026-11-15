"""Differential Pan-STARRS astrometry using identical foreground/background exposures."""
import numpy as np
import pandas as pd

from forecast import HERE, J2016, MAS_PER_DEG, local_axes, parallax_factors


def fit(data):
    design = np.column_stack([np.ones(len(data)), data[:, 2]])
    return np.linalg.lstsq(design, data[:, :2], rcond=None)[0]


def main():
    data = pd.read_csv(HERE / "ps1-detections-2109.csv")
    good = ((data.psfFlux > 5*data.psfFluxErr) & (data.psfFluxErr > 0)
            & (data.psfQfPerfect > 0.85) & data.infoFlag2.isin([0, 128]))
    data = data[good]
    star = data[data.objID == 167922731544072358].set_index("imageID")
    quasar = data[data.objID == 167922731561582287].set_index("imageID")
    pairs = star.join(quasar, how="inner", lsuffix="_l", rsuffix="_q", validate="one_to_one")
    t = pairs.obsTime_l.to_numpy()
    lens = pd.read_csv(HERE / "gaia-2109.csv", index_col="source_id").loc[2124576942440058752]
    offsets = np.column_stack([(pairs.ra_l-pairs.ra_q)*np.cos(np.radians(lens.dec)),
                               pairs.dec_l-pairs.dec_q])*MAS_PER_DEG
    offsets -= lens.parallax*parallax_factors(t, local_axes(lens.ra, lens.dec))
    nights = pd.DataFrame({"east": offsets[:, 0], "north": offsets[:, 1],
                           "dt": (t-J2016)/365.25}).groupby(t.astype(int)).median().to_numpy()
    origin, velocity = fit(nights)
    rng = np.random.default_rng(210905)
    draws = [fit(nights[rng.integers(len(nights), size=len(nights))])[1] for _ in range(10000)]
    dt = -(origin @ velocity) / (velocity @ velocity)
    print("Paired frames, nights:", len(t), len(nights))
    print("Relative east/north motion and formal 95% night bootstrap (mas/yr):",
          velocity, np.percentile(draws, [2.5, 97.5], axis=0))
    print("PS1-only linear minimum (year, mas):", 2016+dt, np.linalg.norm(origin+velocity*dt))


if __name__ == "__main__":
    main()
