"""Reproduce the ASASSN-23hd field registration and pre-eruption photometry."""
import csv
import hashlib
import json
from xml.etree import ElementTree
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"
TARGET = "4662989888116855808"


def read_image(name):
    """Decode the single two-dimensional signed-int16 HDU in these OGLE cutouts."""
    payload = (INPUTS / name).read_bytes()
    header = {}
    for start in range(0, len(payload), 80):
        card = payload[start:start + 80].decode("ascii")
        if card.startswith("END "):
            break
        if card[8:10] == "= ":
            header[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'")
    assert header["BITPIX"] == "16" and header["NAXIS"] == "2"
    width, height = int(header["NAXIS1"]), int(header["NAXIS2"])
    offset = ((start + 80 + 2879) // 2880) * 2880
    raw = np.frombuffer(payload, ">i2", width * height, offset).reshape(height, width)
    return raw * float(header.get("BSCALE", 1)) + float(header.get("BZERO", 0)), header["DATE-OBS"]


def star_positions(data):
    # The reference cutout includes a masked strip set to 65535.
    clean = data.copy()
    clean[clean >= 32000] = np.median(clean[clean < 32000])
    smoothed = gaussian_filter(clean, 1.55)
    peaks = (smoothed == maximum_filter(smoothed, size=9)) & (smoothed > np.median(smoothed) + 20)
    locations = np.argwhere(peaks)
    locations = locations[(locations.min(axis=1) > 8) &
                          ((np.array(data.shape) - locations).min(axis=1) > 8)]
    stars = []
    for y, x in locations:
        yy, xx = np.mgrid[y-5:y+6, x-5:x+6]
        radius2 = (xx - x)**2 + (yy - y)**2
        cut = clean[yy, xx]
        weights = np.maximum(cut - np.median(cut[radius2 >= 25]), 0) * np.exp(-radius2 / 18)
        if weights.sum() > 0:
            stars.append([np.average(xx, weights=weights), np.average(yy, weights=weights)])
    return np.array(stars)


def read_catalog():
    with (INPUTS / "gaia_field.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in ["ra", "dec", "phot_g_mean_mag", "bp_rp", "pmra", "pmdec"]:
            row[key] = float(row[key]) if row[key] else np.nan
    return rows


def coordinates(rows, target, epoch):
    # A local tangent-plane approximation over a 165-arcsec cutout. PMRA includes cos(dec).
    east = (np.array([r["ra"] for r in rows]) - target["ra"]) * np.cos(np.deg2rad(target["dec"])) * 3600
    north = (np.array([r["dec"] for r in rows]) - target["dec"]) * 3600
    pm = np.nan_to_num([[r["pmra"], r["pmdec"]] for r in rows])
    return np.c_[east, north] + pm * (epoch - 2016) / 1000


def plate(data, date, rows, target):
    stars = star_positions(data)
    source_xy = cKDTree(stars)
    epoch = 2000 + (datetime.fromisoformat(date) - datetime(2000, 1, 1, 12)).days / 365.25
    sky = coordinates(rows, target, epoch)
    design = np.c_[sky, np.ones(len(sky))]
    eligible = ((np.array([r["phot_g_mean_mag"] for r in rows]) < 20)
                & (np.linalg.norm(sky, axis=1) > 3))
    # OGLE north-up, east-left, ~0.26 arcsec/pixel; fit an unconstrained affine solution.
    predicted = sky * np.array([-1, 1]) / 0.26 + 320
    for tolerance in [7, 4, 2, 1, 0.8]:
        distances, nearest = source_xy.query(predicted)
        keep = (distances < tolerance) & eligible
        assert keep.sum() >= 15 and len(set(nearest[keep])) == keep.sum()
        measured = stars[nearest[keep]]
        params = np.linalg.lstsq(design[keep], measured, rcond=None)[0]
        predicted = design @ params
    residual = predicted[keep] - measured
    scale = np.sqrt(abs(np.linalg.det(np.linalg.inv(params[:2]))))
    return params, {
        "reference_stars": int(keep.sum()),
        "radial_rms_arcsec": float(scale * np.sqrt(np.mean(np.sum(residual**2, axis=1)))),
        "scale_arcsec_per_pixel": float(scale),
        "epoch": epoch,
    }


def centroid(image, center, half_size):
    # Elliptical Gaussian plus tilted background, with positive source amplitude.
    x0, y0 = np.rint(center).astype(int)
    yy, xx = np.mgrid[y0-half_size:y0+half_size+1, x0-half_size:x0+half_size+1]
    cut = image[yy, xx]
    def residual(par):
        bg, amp, x, y, sx, sy, rho, bx, by = par
        u = (xx - x) / sx
        v = (yy - y) / sy
        gauss = np.exp(-0.5 * (u*u + v*v - 2*rho*u*v) / (1 - rho*rho))
        return (bg + amp*gauss + bx*(xx-x0) + by*(yy-y0) - cut).ravel()
    init = [np.median(cut), cut.max() - np.median(cut), *center, 1.8, 1.8, 0, 0, 0]
    bounds = ([-10000, 0, x0-3, y0-3, 0.8, 0.8, -.7, -20, -20],
              [10000, 10000, x0+3, y0+3, 5, 5, .7, 20, 20])
    result = least_squares(residual, init, bounds=bounds, loss="soft_l1", f_scale=10)
    assert result.success
    return result.x[2:4]


def register(rows, target):
    reference, ref_date = read_image("reference.fits")
    science, obs_date = read_image("science.fits")
    difference, diff_date = read_image("difference.fits")
    assert obs_date == diff_date and obs_date.startswith("2014-10-06")
    results = {}
    for kind, plate_image, plate_date, psf_image in [
        ("difference", reference, ref_date, difference),
        ("science", science, obs_date, science),
        ("quiescent_reference", reference, ref_date, reference),
    ]:
        params, solution = plate(plate_image, plate_date, rows, target)
        centroid_epoch = 2014.763 if kind == "difference" else solution["epoch"]
        target_sky = coordinates([target], target, centroid_epoch)[0]
        projected = np.r_[target_sky, 1.] @ params
        offsets = {}
        for n in [5, 7, 9]:
            measured = centroid(psf_image, projected, n)
            east, north = (measured - projected) @ np.linalg.inv(params[:2])
            offsets[n] = dict(pixel_xy=measured.tolist(), east_arcsec=east, north_arcsec=north,
                              radial_arcsec=np.hypot(east, north))
        results[kind] = dict(**solution, gaia_pixel_xy=projected.tolist(), fits=offsets)
    return results


def light_curve():
    rows = np.loadtxt(INPUTS / "ogle_photometry.dat")
    detections = rows[rows[:, 2] > 0]
    groups = np.split(detections, np.flatnonzero(np.diff(detections[:, 0]) > 45) + 1)
    episodes = []
    for g in groups:
        if len(g) < 2:
            continue
        jd, mag, err = g[np.argmin(g[:, 1]), :3]
        date = datetime(2000, 1, 1, 12, tzinfo=timezone.utc) + timedelta(days=jd-2451545.0)
        episodes.append(dict(detections=len(g), hjd_first=g[0, 0], hjd_last=g[-1, 0],
                             brightest_hjd=jd, brightest_date=date.date().isoformat(),
                             brightest_I=mag, error_mag=err))
    return {"rows":len(rows), "positive_error_rows":len(detections),
            "flagged_rows":int((rows[:, 2] < 0).sum()), "isolated_detections":sum(len(g)==1 for g in groups),
            "episodes":episodes}


def gaia_light_curve():
    table = ElementTree.parse(INPUTS / "gaia_epochs.xml")
    columns = [field.attrib["name"] for field in table.findall(".//{*}FIELD")]
    cells = [item.text.split() for item in table.findall(".//{*}TD")]
    values = dict(zip(columns, cells))
    assert values["source_id"] == [TARGET]
    records = []
    for i in range(int(values["n_transits"][0])):
        if values["variability_flag_g_reject"][i].lower() == "true":
            continue
        time = 2455197.5 + float(values["g_transit_time"][i])
        date = datetime(2000, 1, 1, 12) + timedelta(days=time - 2451545.)
        entry = {"bjd_tcb":time, "date":date.date().isoformat(),
                 "G":float(values["g_transit_mag"][i]),
                 "error_G":float(2.5 / np.log(10) / float(values["g_transit_flux_over_error"][i]))}
        for band in ["bp", "rp"]:
            if values[f"variability_flag_{band}_reject"][i].lower() == "false":
                entry[band.upper()] = float(values[f"{band}_mag"][i])
        records.append(entry)
    # Count distinct visits, retaining paired field-of-view transits as intranight checks.
    groups = []
    for row in records:
        if not groups or row["bjd_tcb"] - groups[-1][-1]["bjd_tcb"] > 1:
            groups.append([])
        groups[-1].append(row)
    visits = [{"date":g[0]["date"], "n_transits":len(g),
               "median_G":float(np.median([r["G"] for r in g]))} for g in groups]
    return {"transits":records, "visits":visits}


def main():
    manifest = json.loads((ROOT / "provenance.json").read_text())
    for filename, source in manifest["inputs"].items():
        actual = hashlib.sha256((INPUTS / filename).read_bytes()).hexdigest()
        if actual != source["sha256"]:
            raise ValueError(f"Input checksum changed for {filename}")
    rows = read_catalog()
    target = next(row for row in rows if row["source_id"] == TARGET)
    other = [row for row in rows if row is not target]
    distances = np.linalg.norm(coordinates(other, target, 2016), axis=1)
    closest = other[np.argmin(distances)]
    report = {"gaia_dr3":TARGET, "gaia_ra":target["ra"], "gaia_dec":target["dec"],
              "images":register(rows, target), "photometry":light_curve(), "gaia_photometry":gaia_light_curve(),
              "local_field":{"other_sources_within_60_arcsec":int((distances <= 60).sum()),
                             "nearest_other_source":closest["source_id"],
                             "nearest_other_arcsec":float(distances.min())}}
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
