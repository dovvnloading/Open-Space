"""Reproduce topocentric pointings and linearized covariance with the orbit fit build."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
BODIES = {"2026 RG39": ("RB15", 27.8), "2026 VV15": ("VE138", 26.0)}


def ephemeris(build, work, state, epoch, designation, magnitude, site, jd):
    source = build / "find_orb"
    output = work / "ephem.txt"
    vector = ",".join(str(value) for value in state)
    command = [
        str(source / "fo"), "-o" + designation.replace(" ", ""),
        f"-v{epoch},{vector},H={magnitude}", f"-x{work}/", f"-e{output}",
        "-E3,5,6,17,20,25", "-C" + site,
        "LINUX_JPL_FILENAME=" + str(build / "de440.bin"), "PERTURBERS=7fe",
        f"EPHEM_START={jd:.12f}", "EPHEM_STEPS=2 1m", "-q",
    ]
    # Explicit -x supplies all configuration. A private child environment omits
    # HOME so Find_Orb does not try to install defaults into the caller's account.
    environment = {k: v for k, v in os.environ.items()
                   if k not in ("HOME", "FO_FORWARD_MODEL")}
    completed = subprocess.run(command, cwd=work, env=environment,
                               capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr[-1000:] or completed.stdout[-1000:])
    values = np.loadtxt(output)
    if values.shape != (2, 14) or abs(values[0, 0] - jd) > 1e-7:
        raise RuntimeError("Unexpected Find_Orb ephemeris format or output epoch")
    return values


def tangent_offset(plus, minus, dec):
    delta_ra = (plus[0] - minus[0] + 180) % 360 - 180
    return np.array([delta_ra * np.cos(np.radians(dec)), plus[1] - minus[1]]) * 3600


def evaluate(record, fit, build):
    designation = record["designation"]
    magnitude = BODIES[designation][1]
    time = datetime.fromisoformat(record["utc"].replace("Z", "+00:00"))
    jd = time.timestamp() / 86400 + 2440587.5
    state = np.asarray(fit["full"]["state"])
    covariance = fit["full"]["covariance"]["site_utc_date_clustered"]
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    with tempfile.TemporaryDirectory() as folder:
        work = Path(folder)
        source = build / "find_orb"
        for data in source.iterdir():
            if data.suffix.lower() in (".txt", ".htm", ".html", ".sof", ".hdr", ".def", ".json", ".pgm"):
                if data.name != "ObsCodes.htm":
                    (work / data.name).symlink_to(data)
        shutil.copyfile(HERE / "stations.txt", work / "ObsCodes.htm")
        shutil.copyfile(source / "environ.def", work / "environ.dat")

        def predict(x):
            return ephemeris(build, work, x, fit["epoch_jd_tdb"], designation,
                             magnitude, record["station"], jd)

        nominal = predict(state)
        row = nominal[0]
        derivatives = []
        for value, direction in zip(eigenvalues, eigenvectors.T):
            perturbation = np.sqrt(max(0, value)) * direction
            plus = predict(state + perturbation)[0, 1:3]
            minus = predict(state - perturbation)[0, 1:3]
            derivatives.append(tangent_offset(plus, minus, row[2]) / 2)
        jacobian = np.array(derivatives)
        variances = np.linalg.eigvalsh(jacobian.T @ jacobian)
        semiaxes = np.sqrt(variances[::-1])

    nominal_offset = tangent_offset(row[1:3], [record["ra_deg"], record["dec_deg"]], row[2])
    expected = record["linearized_cluster_sigma_major_minor_arcsec"]
    if np.linalg.norm(nominal_offset) > 0.15 or not np.allclose(semiaxes, expected, atol=0.08):
        raise ValueError(f"{designation} forecast differs from the archive")
    print(f"{designation} {record['station']} {record['utc']} "
          f"RA={row[1]:.6f} Dec={row[2]:.6f} V={row[7]:.1f} "
          f"rate={row[9]:.3f} arcsec/min, sigma axes={semiaxes[0]:.3f},{semiaxes[1]:.3f} arcsec")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path, help="Directory produced by the orbit build.sh")
    parser.add_argument("--fits", type=Path, default=HERE.parent / "orbits" / "fits")
    args = parser.parse_args()
    packet = json.loads((HERE / "followup.json").read_text())
    for row in packet["predictions"]:
        key = BODIES[row["designation"]][0]
        fit = json.loads((args.fits / (key + ".json")).read_text())
        evaluate(row, fit, args.build.resolve())


if __name__ == "__main__":
    main()
