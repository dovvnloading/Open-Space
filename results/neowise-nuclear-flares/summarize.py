"""Combine flare measurements, centroid shifts and catalogue records.

Usage: python summarize.py --flares F.csv --epochs E.csv --positions P.csv --centroid C.csv
                           --vet VET.json --flaires FLAIRES.csv --masterson PUBLISHED.csv --out OUT.csv
                           [--notes NOTES.json]

Adds a distance (cz / 70 km/s/Mpc), the blackbody temperature and luminosity of
the excess emission at peak, the blackbody energy radiated over the observed
visits, and a class:

  published           within 5 arcsec of a flare in Masterson et al. (2024),
                      Jiang et al. (2021) or Necker et al. (2025), SIMBAD lists
                      a WTP transient there, or a hand-entered note says so
  reported transient  a TNS object within 10 arcsec discovered between two
                      years before and one year after the infrared onset
  off-nuclear         centroid moved by more than 4 sigma and the implied
                      transient position is more than 1 arcsec from the nucleus
  active nucleus      SIMBAD classifies the galaxy as an AGN, its quiescent
                      W1 - W2 exceeds 0.5 mag, or a ZTF alert object within
                      5 arcsec has 30 or more detections
  unreported          none of the above

NOTES.json holds hand-entered findings: {id: {"class": ..., "note": ..., "z": ...}}.
A class given there replaces the automatic one; `z` supplies a redshift missing
from 2MRS.
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq

H, C, K, SIGMA = 6.62607e-27, 2.99792e10, 1.38065e-16, 5.67037e-5
NU1, NU2 = C / 3.368e-4, C / 4.618e-4
MPC = 3.0857e24
AGN_TYPES = {"Sy1", "Sy2", "SyG", "AGN", "AG?", "LIN", "Bla", "BLL", "QSO", "rG", "BL?", "Q?"}


def planck(nu, t):
    return 2 * H * nu ** 3 / C ** 2 / np.expm1(H * nu / (K * t))


def blackbody(f1_mjy, f2_mjy, mpc):
    """Temperature (K) and luminosity (erg/s) of a blackbody through two fluxes."""
    if not (f1_mjy > 0 and f2_mjy > 0):
        return np.nan, np.nan
    ratio = f2_mjy / f1_mjy
    try:
        t = brentq(lambda temp: planck(NU2, temp) / planck(NU1, temp) - ratio, 150, 20000)
    except ValueError:
        return np.nan, np.nan
    omega = f2_mjy * 1e-26 / planck(NU2, t)
    return t, 4 * (mpc * MPC) ** 2 * omega * SIGMA * t ** 4


def separation(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = (np.radians(x) for x in (ra1, dec1, ra2, dec2))
    cos = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(cos, -1, 1))) * 3600


def radiated_energy(visits, f, mpc):
    """Blackbody energy (erg) of the excess over the visits from onset on."""
    after = visits[visits.mjd >= f.onset_mjd - 1].sort_values("mjd")
    lum = np.array([blackbody(a - f.level_w1_mjy, b - f.level_w2_mjy, mpc)[1]
                    for a, b in zip(after.f1, after.f2)])
    good = np.isfinite(lum)
    if good.sum() < 2:
        return np.nan
    t = after.mjd.values[good] * 86400.0
    return float(np.sum(0.5 * (lum[good][1:] + lum[good][:-1]) * np.diff(t)))


def main():
    parser = argparse.ArgumentParser()
    for name in ("flares", "epochs", "positions", "centroid", "vet", "flaires", "masterson", "out"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--notes")
    parser.add_argument("--ztf", help="output of ztf_alerts.py")
    args = parser.parse_args()
    flares = pd.read_csv(args.flares, dtype={"id": str})
    positions = pd.read_csv(args.positions, dtype={"id": str}).set_index("id")
    centroid = pd.read_csv(args.centroid, dtype={"id": str}).set_index("id")
    vet = json.load(open(args.vet))
    flaires = pd.read_csv(args.flaires)
    masterson = pd.read_csv(args.masterson)
    notes = json.load(open(args.notes)) if args.notes else {}
    ztf = json.load(open(args.ztf)) if args.ztf else {}
    epochs = pd.read_csv(args.epochs, dtype={"id": str})
    epochs = {name: group for name, group in epochs[epochs.id.isin(flares.id)].groupby("id")}
    rows = []
    for _, f in flares.iterrows():
        ra, dec, cz = positions.loc[f.id, ["ra", "dec", "cz"]]
        note = notes.get(f.id, {})
        if "z" in note:
            cz = note["z"] * 299792.458
        mpc = cz / 70.0 if cz > 300 else np.nan
        t, lum = blackbody(f.peak_excess_w1_mjy, f.peak_excess_w2_mjy, mpc)
        nu_l_nu = 4 * np.pi * (mpc * MPC) ** 2 * NU2 * f.peak_excess_w2_mjy * 1e-26
        record = vet.get(f.id, {"simbad": [], "tns": []})
        simbad = record["simbad"]
        host = next((o for o in simbad if o["type"] not in ("HII", "X", "Rad", "UV", "Cl*", "SN*", "ev", "EB*", "OpC")),
                    simbad[0] if simbad else {"name": "", "type": ""})
        onset_year = 2000 + (f.onset_mjd - 51544.5) / 365.25
        quiet_year = 2000 + (f.last_quiet_mjd - 51544.5) / 365.25
        transients = []
        for o in record["tns"]:
            year = float(o["discovered"][:4]) + (float(o["discovered"][5:7]) - 0.5) / 12 if o["discovered"] else np.nan
            if -2.0 < year - onset_year < 1.0:
                transients.append(f"{o['name']} ({o['type'] or 'unclassified'}, {o['discovered'][:10]})")
        wtp = [o["name"] for o in simbad if o["name"].startswith("WTP")]
        near_m = masterson[separation(ra, dec, masterson.ra.values, masterson.dec.values) < 5]
        near_f = flaires[separation(ra, dec, flaires.ra.values, flaires.dec.values) < 5]
        published = wtp + list(near_m.name) + [f"Flaires {n.strip()}" for n in near_f.Name]
        c = centroid.loc[f.id] if f.id in centroid.index else None
        off = c is not None and c.shift_arcsec > 4 * c.shift_error_arcsec and c.implied_transient_offset_arcsec > 1.0
        alerts = ztf.get(f.id)
        most = max((a["detections"] or 0 for a in alerts), default=0) if alerts is not None else np.nan
        near = [a for a in (alerts or []) if (a["detections"] or 0) >= 3 and a["first_mjd"] > f.onset_mjd - 365
                and a["first_mjd"] < f.onset_mjd + 200 and a["last_mjd"] - a["first_mjd"] < 800]
        agn = host["type"] in AGN_TYPES or f.quiescent_w1_w2 > 0.5 or most >= 30
        if published:
            kind = "published"
        elif transients:
            kind = "reported transient"
        elif off:
            kind = "off-nuclear"
        elif agn:
            kind = "active nucleus"
        else:
            kind = "unreported"
        kind = note.get("class", kind)
        rows.append({
            "id": f.id, "host": " ".join(host["name"].split()), "simbad_type": host["type"], "ra": ra, "dec": dec,
            "cz_kms": cz, "distance_mpc": mpc, "class": kind, "last_quiet_year": quiet_year, "onset_year": onset_year,
            "peak_year": 2000 + (f.peak_mjd - 51544.5) / 365.25, "visits_after_onset": f.visits_after_onset,
            "quiescent_w1_w2": f.quiescent_w1_w2, "baseline_scatter_w2": f.base_scatter_w2,
            "onset_sigma_w2": f.onset_z_w2, "peak_excess_w1_mjy": f.peak_excess_w1_mjy,
            "peak_excess_w2_mjy": f.peak_excess_w2_mjy, "peak_rise_w2_mag": f.peak_dmag_w2,
            "excess_temperature_k": t, "excess_bb_luminosity_erg_s": lum, "nu_l_nu_w2_erg_s": nu_l_nu,
            "excess_bb_energy_erg": radiated_energy(epochs[f.id], f, mpc),
            "last_over_peak_excess": f.last_excess_fraction, "n_peaks": f.n_peaks,
            "centroid_shift_arcsec": c.shift_arcsec if c is not None else np.nan,
            "centroid_shift_error_arcsec": c.shift_error_arcsec if c is not None else np.nan,
            "implied_offset_arcsec": c.implied_transient_offset_arcsec if c is not None else np.nan,
            "published_as": "; ".join(published), "tns_near_onset": "; ".join(transients),
            "tns_any": "; ".join(f"{o['name']} ({o['type'] or 'unclassified'}, {o['discovered'][:10]})"
                                 for o in record["tns"]),
            "ztf_max_detections": most,
            "ztf_burst_near_onset": "; ".join(f"{a['oid']} ({a['detections']} detections, MJD {a['first_mjd']:.0f}-{a['last_mjd']:.0f})"
                                              for a in near),
            "note": note.get("note", ""), "vetted": f.id in vet})
    out = pd.DataFrame(rows).sort_values(["class", "nu_l_nu_w2_erg_s"], ascending=[True, False])
    out.to_csv(args.out, index=False, float_format="%.7g")
    print(out["class"].value_counts().to_dict())


if __name__ == "__main__":
    main()
