"""Replay the CGCG 476-019 nuclear light curves from frozen public measurements."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "inputs"
TARGET = "J232141.38+262909.4"
DISCOVERY_MJD = 59422
PERIOD = 365.25

def angular_distance(ra, dec, ra_center, dec_center):
    """Great-circle distance in degrees, vectorized."""
    dra = np.radians(ra - ra_center)
    ddec = np.radians(dec - dec_center)
    hav = np.sin(ddec/2)**2 + np.cos(np.radians(dec))*np.cos(np.radians(dec_center))*np.sin(dra/2)**2
    return np.degrees(2*np.arcsin(np.sqrt(np.clip(hav, 0, 1))))


def infrared_visits():
    data = pd.read_csv(DATA / "neowise.csv.gz", dtype={"cc_flags": str, "moon_masked": str})
    good = data[(data.cc_flags == "0000") & (data.ph_qual == "AA") &
                (data.moon_masked == "00") & (data.qi_fact >= 0.5)].copy()
    good = good.sort_values(["source", "mjd"])
    good["visit"] = good.groupby("source")["mjd"].transform(lambda x: (x.diff() > 45).cumsum())
    visits = good.groupby(["source", "visit"]).agg(
        mjd=("mjd", "median"), n=("mjd", "size"),
        w1=("w1mpro", "median"), w2=("w2mpro", "median")).reset_index()
    # Opposite visits have different scan angles. Compare each visit to earlier
    # visits in the same annual season, giving each visit one vote.
    for source, subset in visits.groupby("source"):
        baseline = subset[(subset.mjd < DISCOVERY_MJD) & (subset.n >= 5)]
        for row in subset.itertuples():
            phase = (baseline.mjd - row.mjd + PERIOD/2) % PERIOD - PERIOD/2
            reference = baseline[phase.abs() < 60]
            visits.loc[row.Index, "baseline_visits"] = len(reference)
            for band in ["w1", "w2"]:
                median = reference[band].median()
                visits.loc[row.Index, band+"_excess_mag"] = median - getattr(row, band)
                visits.loc[row.Index, band+"_baseline_mad"] = np.median(abs(reference[band]-median))
    return visits


def direct_photometry():
    data = pd.read_csv(DATA / "ztf-direct.csv.gz")
    data = data[(data.catflags == 0) & (data.magerr < 0.1)].copy()
    data["year"] = pd.to_datetime(data.mjd, unit="D", origin="1858-11-17").dt.year
    return data


def alert_photometry():
    return pd.read_csv(DATA / "ztf-alerts.csv.gz", dtype={"candid": str})


def main():
    visits = infrared_visits()
    visits.to_csv(ROOT / "infrared-visits.csv", index=False, float_format="%.6f")
    target = visits[visits.source == TARGET]
    print("NEOWISE quality measurements/visits:", target.n.sum(), len(target))
    print(target[target.mjd > DISCOVERY_MJD][["mjd","n","w1","w2","w1_excess_mag","w2_excess_mag"]].round(3).to_string(index=False))
    key = visits[visits.mjd.between(60115,60125) | visits.mjd.between(60276,60286)]
    # The narrow-band quantity nu*L_nu at emitted (1+z)*nu_obs equals
    # 4*pi*D_L^2*nu_obs*F_nu,obs. Use the nominal W2 isophotal wavelength.
    peak = key[key.source == TARGET]
    f_jy = 171.787 * 10**(-0.4*peak.w2) * (1 - 10**(-0.4*peak.w2_excess_mag))
    distance_cm = 134.8 * 3.08567758e24
    excess_lum = 4*np.pi*distance_cm**2*(2.99792458e14/4.6)*f_jy*1e-23
    print("W2 excess mJy:", (1000*f_jy).round(3).tolist(),
          "nu*L_nu (1e42 erg/s):", (excess_lum/1e42).round(3).tolist())
    controls = key[(key.source != TARGET) & (key.n >= 5)]
    print("Control visits, largest |W1|/|W2|:", len(controls),
          controls[["w1_excess_mag","w2_excess_mag"]].abs().max().round(3).tolist())
    for cat in ["flaires", "mirong"]:
        rows = pd.read_csv(DATA / (cat+"-positions.csv.gz"))
        sep = angular_distance(rows.ra, rows.dec, 350.422429, 26.485967)
        print(f"{cat}: {len(rows)} rows, nearest {rows.iloc[np.argmin(sep)]['name'].strip()} at {min(sep):.4f} deg")
    direct = direct_photometry()
    annual = direct.groupby(["field","filtercode","year"]).agg(n=("mag","size"),median_mag=("mag","median"))
    annual.to_csv(ROOT / "ztf-annual.csv", float_format="%.6f")
    alerts = alert_photometry()
    print("ZTF direct quality rows / ZTF alerts:", len(direct), len(alerts))
    print("Field 647 median magnitudes:")
    print(annual.loc[647].median_mag.unstack("year").round(3).to_string())
    red = alerts[(alerts.fid == 2) & (alerts.isdiffpos == 1)]
    for start,stop in [(59422,59787),(59945,60310),(60676,61013)]:
        sample = red[red.mjd.between(start,stop)]
        ra, dec = sample[["ra","dec"]].median()
        radius = angular_distance(sample.ra,sample.dec,ra,dec)*3600
        print(f"r-band {start}–{stop}: N={len(sample)}, median RA={ra:.7f}, Dec={dec:.7f}, median radial scatter={np.median(radius):.3f} arcsec")
    assert len(target) == 22 and target.n.sum() == 309
    assert (key[key.source == TARGET].w2_excess_mag > 0.45).all()
    assert (controls.w2_excess_mag.abs() < 0.03).all()
    assert len(alerts)==733 and len(direct)==2724


if __name__ == "__main__":
    main()
