"""Extract foreground spectroscopy, photo-cluster matches, and same-cone velocity controls."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
import numpy as np
from search import C, directions, write_csv

PIXELS = (17996, 17997, 17998, 17999)
BASE = "https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/main/dark/179"


def angular_distance(positions, center):
    return 60 * np.rad2deg(2 * np.arcsin(np.linalg.norm(positions-center, axis=-1)/2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("redrock_dir", type=Path)
    p.add_argument("--out", type=Path, default=Path("."))
    args = p.parse_args()
    out = args.out
    absrows = np.genfromtxt(out / "j0108-absorbers.csv", delimiter=",", names=True, dtype=None)
    qra, qdec = absrows["RA"], absrows["DEC"]
    center = directions(qra.mean(), qdec.mean())
    z0 = np.median(absrows["Z_MGII"])
    cosmology = FlatLambdaCDM(H0=70, Om0=.3)
    scale = cosmology.kpc_proper_per_arcmin(z0).value / 1000
    records = []
    checks = {}
    names = ["targetid", "ra", "dec", "z", "zerr", "deltachi2", "zwarn", "spectype", "healpix"]
    for pixel in PIXELS:
        name = f"redrock-main-dark-{pixel}.fits"
        path = args.redrock_dir / name
        checks[f"{BASE}/{pixel}/{name}"] = hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as hdus:
            spec = hdus["REDSHIFTS"].data
            fibers = hdus["FIBERMAP"].data
            if not np.array_equal(spec["TARGETID"], fibers["TARGETID"]):
                raise ValueError("Unaligned redshift and fibermap rows")
            records.extend(zip(spec["TARGETID"], fibers["TARGET_RA"], fibers["TARGET_DEC"],
                               spec["Z"], spec["ZERR"], spec["DELTACHI2"], spec["ZWARN"],
                               spec["SPECTYPE"], np.full(len(spec), pixel)))
    cat = np.rec.fromrecords(records, names=names)
    sky = directions(cat.ra, cat.dec)
    radius = angular_distance(sky, center)
    impacts = np.column_stack([angular_distance(sky, d) for d in directions(qra, qdec)])
    dv = C * np.log((1+cat.z)/(1+z0))
    secure = (cat.spectype == "GALAXY") & (cat.zwarn == 0) & (cat.deltachi2 >= 25)
    sel = np.flatnonzero(secure & (radius <= 15) & (np.abs(dv) <= 500))
    write_csv(out / "foreground-galaxies.csv",
              [*names, "center_arcmin", "dv_kms", "impact_A_Mpc", "impact_B_Mpc", "impact_C_Mpc"],
              ((*cat[i], radius[i], dv[i], *(impacts[i]*scale)) for i in sel[np.argsort(radius[sel])]))
    controls = []
    for r in (3, 6, 10, 12, 15):
        vs = dv[secure & (radius <= r)]
        for halfwidth in (250, 500, 1000):
            offsets = np.arange(2*halfwidth, 15000, 2*halfwidth)
            for offset in np.r_[-offsets[::-1], 0, offsets]:
                controls.append((r, halfwidth, offset, np.count_nonzero(abs(vs-offset) <= halfwidth)))
    write_csv(out / "galaxy-control-windows.csv",
              ["radius_arcmin", "velocity_halfwidth_kms", "velocity_offset_kms", "galaxies"], controls)
    with (out / "wen-han-cone.tsv").open() as stream:
        lines = [line for line in stream if line.strip() and not line.startswith("#")]
    cluster_rows = csv.DictReader(io.StringIO(lines[0] + "".join(lines[3:])), delimiter="\t")
    matches = []
    for r in cluster_rows:
        zc = float(r["zCl"])
        if abs(zc-z0) > .04:
            continue
        central = directions(float(r["RAJ2000"]), float(r["DEJ2000"]))
        separations = angular_distance(sky, central)*60
        j = np.argmin(separations)
        match = [cat.targetid[j], cat.z[j], cat.zwarn[j], cat.deltachi2[j], separations[j], dv[j]]
        if separations[j] > 2:
            match = [""]*len(match)
        impact = angular_distance(directions(qra, qdec), central)*scale
        matches.append([r["ID"].strip(), r["Name"].strip(), r["RAJ2000"].strip(),
                        r["DEJ2000"].strip(), zc, r["f_zCl"].strip(), r["r500"].strip(),
                        r["Ngal"].strip(), r["Cat"].strip(), *impact, *match])
    write_csv(out / "cluster-crossmatches.csv",
              ["wen_han_id", "name", "ra_deg", "dec_deg", "catalog_z", "catalog_z_is_spec",
               "catalog_R500_Mpc", "catalog_Ngal", "previous_catalog", "impact_A_Mpc",
               "impact_B_Mpc", "impact_C_Mpc", "desi_id", "desi_z", "desi_zwarn",
               "desi_deltachi2", "desi_sep_arcsec", "desi_dv_kms"], matches)
    geometry = []
    qpos = directions(qra, qdec)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        distance = angular_distance(qpos[i], qpos[j])
        geometry.append([i+1, j+1, distance, distance*scale,
                         C*abs(np.log((1+absrows['Z_MGII'][i])/(1+absrows['Z_MGII'][j])))])
    write_csv(out / "sightline-separations.csv",
              ["sightline_i", "sightline_j", "separation_arcmin", "proper_Mpc", "dv_kms"], geometry)
    summary = {"center_ra_deg": qra.mean(), "center_dec_deg": qdec.mean(),
               "z_reference": z0, "proper_Mpc_per_arcmin": scale, "H0_kms_Mpc": 70, "Omega_m": .3,
               "input_urls_sha256": checks}
    (out / "foreground-provenance.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
