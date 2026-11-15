"""Reconstruct nominal ring and limb crossings from the frozen public ephemerides."""
import csv
import gzip
import hashlib
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"
AU_KM = 149597870.7
# NASA PDS Ring-Moon Systems Node, Vital Statistics for Saturn's Rings, km.
A_OUTER, A_INNER = 136770., 122050.
ENCKE_OUTER, ENCKE_INNER = 133745., 133423.
# One-bar equatorial and polar reference radii in km, NAIF pck00011.tpc.
EQUATORIAL, POLAR = 60268., 54364.
J2000 = datetime(2000, 1, 1, 12, tzinfo=timezone.utc)


def sexagesimal(value):
    a, b, c = map(float, value.lstrip("+-").split(":"))
    return (-1 if value.startswith("-") else 1) * (a + b / 60 + c / 3600)


def direction(ra, dec):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.stack([np.cos(r) * np.cos(d), np.sin(r) * np.cos(d), np.sin(d)], axis=-1)


def load_ephemeris(entry):
    raw = gzip.decompress((EVIDENCE / entry["file"]).read_bytes())
    if hashlib.sha256(raw).hexdigest() != entry["sha256_uncompressed"]:
        raise ValueError(f'Input digest mismatch: {entry["file"]}')
    if entry["provider"] == "horizons":
        section = raw.decode().split("$$SOE")[1].split("$$EOE")[0].strip()
        rows = list(csv.reader(section.splitlines()))
        topocentric = entry["station"] != "geocenter"
        fields = [0, 3, 4, 8 if topocentric else 6]
        jd, ra, dec, distance = np.array([[float(r[i]) for i in fields] for r in rows]).T
        elevation = np.array([float(r[6]) for r in rows]) if topocentric else None
    else:
        doc = json.loads(raw)
        rows = doc["data"]
        frame = doc["ephemeris"]["reference_frame"]
        if doc["ephemeris"]["planetary_theory"] != "INPOP19" or "astrometric J2000" != frame["type"]:
            raise ValueError("Unexpected Miriade model or frame")
        expected_center = "geocenter" if entry["station"]=="geocenter" else "topocenter"
        if frame["center"] != expected_center:
            raise ValueError("Unexpected Miriade observing site")
        jd = np.array([2451545 + (datetime.fromisoformat(r["Date"]).replace(tzinfo=timezone.utc)
                                  - J2000).total_seconds()/86400 for r in rows])
        ra = np.array([15 * sexagesimal(r["RA"]) for r in rows])
        dec = np.array([sexagesimal(r["DEC"]) for r in rows])
        distance = np.array([r["Dobs"] for r in rows])
        elevation = None
    epoch = J2000 + timedelta(days=float(jd[0] - 2451545))
    seconds = (jd - jd[0]) * 86400
    observer_to_saturn = direction(ra, dec) * (AU_KM * distance[:, None])
    return epoch, seconds, observer_to_saturn, elevation


def surfaces(saturn, source, pole):
    # B is the closest point of the source ray to Saturn's center. The
    # ring intersection is a small number of km along this ray, avoiding
    # subtracting two approximately billion-km vectors to obtain its norm.
    b = (saturn @ source)[:, None] * source - saturn
    along_ring = -(b @ pole) / (source @ pole)
    ring_points = b + along_ring[:, None] * source
    radius = np.linalg.norm(ring_points, axis=1)
    # Minimize the ellipsoid quadratic form over this same source ray.
    c = (EQUATORIAL / POLAR)**2 - 1
    z = source @ pole
    along_limb = -(b @ source + c * (b @ pole) * z) / (1 + c*z*z)
    nearest = b + along_limb[:, None] * source
    limb = np.sqrt(np.sum(nearest**2, axis=1) + c*(nearest @ pole)**2)/EQUATORIAL
    return radius, limb


def crossings(seconds, values, level):
    spline = CubicSpline(seconds, values-level)
    return [brentq(spline, seconds[i], seconds[i+1]) for i in
            np.flatnonzero(np.signbit(values[:-1]-level) != np.signbit(values[1:]-level))]


def local_minimum(seconds, values):
    return minimize_scalar(CubicSpline(seconds, values), bounds=(seconds[0], seconds[-1]),
                           method="bounded", options={"xatol": 1.e-6})


def format_utc(epoch, seconds):
    return (epoch + timedelta(seconds=float(seconds)+0.5)).strftime("%H:%M:%S")


def recover_1989(entry):
    epoch, seconds, saturn, _ = load_ephemeris(entry)
    centuries = (epoch-J2000).total_seconds()/(36525*86400)
    pole = direction(40.589-0.036*centuries, 83.537-0.004*centuries)
    source = direction(15*sexagesimal("18:46:20.58417"),sexagesimal("-22:23:31.87202"))
    radius, _ = surfaces(saturn, source, pole)
    predicted = crossings(seconds, radius, 117530.)
    # Published observed ingress and egress, in seconds after 05:00 UTC.
    observed = [3600 + 16*60 + 26.670, 4*3600 + 10*60 + 42.380]
    return {"year":"1989", "site":"PALOMAR", "nominal_B_edge_km":117530,
            "measured_utc":["06:16:26.670", "09:10:42.380"],
            "predicted_utc":[format_utc(epoch,t) for t in predicted],
            "predicted_minus_measured_sec":[round(p-o,2) for p,o in zip(predicted,observed)]}


def measure(entry, position):
    epoch, sec, saturn, elevation = load_ephemeris(entry)
    centuries = (epoch-J2000).total_seconds()/(36525*86400)
    pole = direction(40.589-0.036*centuries, 83.537-0.004*centuries)
    ra, dec, erra, edec, corr = position
    radius, limb = surfaces(saturn, direction(ra, dec), pole)
    near = local_minimum(sec, radius)
    t_limb = local_minimum(sec, limb)
    levels = {"a_outer": A_OUTER, "a_inner": A_INNER, "encke_outer": ENCKE_OUTER,
              "encke_inner": ENCKE_INNER}
    times = {label: crossings(sec, radius, level) for label, level in levels.items()}
    result = {"year": entry["year"], "site": entry["station"], "provider": entry["provider"],
              "min_ring_km": round(near.fun, 2), "min_ring_utc": format_utc(epoch, near.x),
              "min_limb_ratio": round(t_limb.fun, 6),
              "one_bar_limb": [format_utc(epoch, t) for t in crossings(sec, limb, 1)]}
    for label, ts in times.items():
        result[label] = [format_utc(epoch, t) for t in ts]
    if entry["year"] == "2030":
        # Test the size of the foreground monopole bending, excluded in the
        # straight-ray solution. The thin-lens image is radially displaced on
        # the sky; at origins within one ring radius it is locally straight
        # to sub-km accuracy. This excludes refraction, multipoles and waves.
        source = direction(ra, dec)
        impact = (saturn @ source)[:, None]*source - saturn
        impact2 = np.sum(impact**2, axis=1)
        einstein2 = 4 * 37931207.8 * (saturn @ source) / 299792.458**2
        scale = (1 + np.sqrt(1 + 4*einstein2/impact2))/2
        bent, _ = surfaces(saturn - (scale-1)[:, None]*impact, source, pole)
        result["monopole_min_ring_km"] = round(local_minimum(sec, bent).fun, 2)
        for name, edge in levels.items():
            result["monopole_"+name] = [format_utc(epoch,q) for q in crossings(sec,bent,edge)]
    if elevation is not None and times["a_outer"]:
        result["a_edge_elevation_deg"] = [round(float(np.interp(t,sec,elevation)),2)
                                          for t in times["a_outer"]]
    # Sample the boundary of a radius-three ellipse in catalogue coordinate
    # error units. RFC gives sigma-alpha without cos(delta); offset RA itself.
    if times["a_outer"] and entry["year"] == "2030":
        chol = np.linalg.cholesky([[erra**2,corr*erra*edec],[corr*erra*edec,edec**2]])
        differences = []
        for angle in np.arange(72)*2*np.pi/72:
            dra, ddec = (3/3600000) * (chol @ np.array([np.cos(angle),np.sin(angle)]))
            rr, _ = surfaces(saturn, direction(ra+dra,dec+ddec),pole)
            new = crossings(sec, rr, A_OUTER)
            if len(new) != len(times["a_outer"]):
                raise ValueError("Astrometric sensitivity changed a chord classification")
            differences.extend(np.abs(np.array(new)-times["a_outer"]))
        result["max_catalog_3sigma_edge_shift_sec"] = round(float(max(differences)),2)
    return result


def main():
    catalog = json.loads((EVIDENCE/"catalog.json").read_text())
    positions = {}
    for line in catalog["rfc_2026c_cat.txt"]["lines"]:
        r = line.split()
        positions[r[1]] = (15*sexagesimal(":".join(r[3:6])), sexagesimal(":".join(r[6:9])),
                           *map(float,r[9:12]))
    manifest = json.loads((EVIDENCE/"manifest.json").read_text())
    output = [recover_1989(row) if row["year"]=="1989" else measure(row, positions["J0413+1903" if row["year"]=="2030" else "J0225+1134"])
              for row in manifest["ephemerides"]]
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
