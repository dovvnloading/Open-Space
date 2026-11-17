"""Reflected SDSS r and monochromatic greybody fluxes for explicit planet models."""
import csv
from pathlib import Path

import numpy as np
from scipy.constants import h, c, k, astronomical_unit

# (mass [M_Earth], radius [R_Earth], T_eff [K]).
# Fortney et al., arXiv:1604.07424 Table 1; final row Linder & Mordasini nominal model.
MODELS = [(5, 4.12, 36.7), (5, 2.94, 42.2), (5, 2.71, 38.9),
          (8, 4.44, 39.5), (8, 3.44, 46.2), (8, 3.17, 42.6),
          (10, 5.09, 40.3), (10, 3.46, 48.3), (10, 3.16, 45.1),
          (10, 3.7, 47.)]


def main():
    target = Path(__file__).with_name("brightness.csv")
    with target.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mass_Mearth","radius_Rearth","teff_K","distance_au",
                    "geometric_albedo_r","r_AB","F60_Jy","F65_Jy","F90_Jy","F100_Jy"])
        for mass, radius, teff in MODELS:
            r_au = radius * 6371e3 / astronomical_unit
            for d in (500, 600, 700):
                nu = c / (np.array([60,65,90,100])*1e-6)
                thermal = 1e26 * np.pi * (r_au/d)**2 * (2*h*nu**3/c**2) / np.expm1(h*nu/(k*teff))
                for albedo in (0.1, 0.3, 0.75):
                    # Solar apparent SDSS r_AB at 1 au (Willmer 2018). Full phase.
                    rmag = -26.93 - 2.5*np.log10(albedo*(r_au/d**2)**2)
                    w.writerow([mass,radius,teff,d,albedo,f"{rmag:.4f}",*[f"{flux:.6f}" for flux in thermal]])


if __name__ == "__main__":
    main()
