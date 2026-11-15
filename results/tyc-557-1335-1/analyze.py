"""Calculate the single-orbit falsification and target/control display-video tests."""
import csv
import json
import math
from pathlib import Path
from statistics import mean

BASE = Path(__file__).resolve().parent / "inputs"


def rows(filename):
    with (BASE / filename).open() as handle:
        return list(csv.DictReader(handle))


def angular_distance(ra1, dec1, ra2, dec2):
    """Haversine distance in arcseconds."""
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    h = math.sin((dec2-dec1)/2)**2 + math.cos(dec1)*math.cos(dec2)*math.sin((ra2-ra1)/2)**2
    return math.degrees(2*math.asin(math.sqrt(h)))*3600


def main():
    source = rows("gaia.csv")[0]
    p, dp, e, de, pl, dpl, k = (float(source[key]) for key in
        ("period", "period_error", "eccentricity", "eccentricity_error", "parallax", "parallax_error", "semi_amplitude_primary"))
    geometry = []
    for row in rows("geometry.csv"):
        rho, drho = float(row["rho_mas"]), float(row["rho_sigma_mas"])
        s = rho / pl
        # Projection only decreases the observed separation; apastron is a(1+e).
        mass_at_e = (s/(1+e))**3/(p/365.25)**2
        weakest_s = (rho-3*drho)/(pl+3*dpl)
        geometry.append(dict(pa_deg=float(row["pa_deg"]), rho_mas=rho, projected_au=s,
                             minimum_mass_msun_at_e=mass_at_e,
                             adverse_three_sigma_mass_msun_any_e=(weakest_s/2)**3/((p+3*dp)/365.25)**2))
    n, ngof, rv_range, gof, snr = (float(source[key]) for key in
        ("rv_n_good_obs_primary", "rv_renormalised_gof", "rv_amplitude_robust", "goodness_of_fit", "rv_expected_sig_to_noise"))
    z = (3.8791 - .6062*math.log10(p) - .041*n - 2.9939*math.log10(ngof)
         + .0153*rv_range/(2*k) + .4038*gof + .9734*math.log10(snr))
    neighbors = [(angular_distance(float(source["ra"]), float(source["dec"]), float(r["ra"]), float(r["dec"])), r)
                 for r in rows("field_60arcsec.csv") if r["source_id"] != source["source_id"]]
    closest, nbr = min(neighbors)
    bright = rows("field_bright_5arcmin.csv")
    bright_comp = [r for r in bright if r["source_id"] != source["source_id"] and float(r["phot_g_mean_mag"]) <= float(source["phot_g_mean_mag"])+2]
    video = rows("video_photometry.csv")
    # Half-frame contacts occur at about frames 51.5, 63.5 and 86.5, 89.
    # Exclude all edge samples from the baseline and use fully in-dip frames only.
    baseline = set(range(120)) - set(range(50,66)) - set(range(85,91))
    primary, secondary = list(range(52,64)), [87,88]
    scan_shift = .146 * max(row["rho_mas"] for row in geometry)
    phot = {}
    for radius in (2,3,4):
        block = {}
        for star in ("target","ref1","ref2","ref3"):
            flux=[float(row[f"{star}_r{radius}"]) for row in video]
            normal=mean(flux[i] for i in sorted(baseline))
            second=mean(flux[i] for i in secondary)
            pairs=[mean(flux[i:i+2]) for i in range(119) if i in baseline and i+1 in baseline]
            block[star] = dict(primary=mean(flux[i] for i in primary)/normal, secondary=second/normal,
                               lower_control_windows=sum(test <= second for test in pairs), control_windows=len(pairs))
        phot[str(radius)] = block
    print(json.dumps(dict(geometry=geometry, bashi_score=1/(1+math.exp(z)),
                          mass_function_msun=1.036149e-7*(1-e**2)**1.5*k**3*p,
                          nearest_neighbor_arcsec=closest,
                          nearest_neighbor_delta_g=float(nbr["phot_g_mean_mag"])-float(source["phot_g_mean_mag"]),
                          bright_neighbors_5arcmin=len(bright_comp),
                          projected_rvs_shift_max_km_s_at_occultation=scan_shift,
                          rvs_k_to_geometric_shift_ratio=k/scan_shift, video_relative_mean_flux=phot), indent=2))


if __name__ == "__main__":
    main()
