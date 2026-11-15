"""Back-project measured angular motions against Green's October 2024 SNR centers."""
import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import tarfile
from urllib.request import urlretrieve

import numpy as np

HERE = Path(__file__).resolve().parent
ATNF = "https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz"
GREEN = "https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html"
MAS_PER_RAD = 180 / np.pi * 3_600_000
SECONDS_YEAR = 365.25 * 86400


def sexagesimal(value):
    parts = value.split(":")
    return (-1 if value.startswith("-") else 1) * sum(abs(float(x)) / 60**i for i, x in enumerate(parts))


def basis(lon, lat):
    a, d = np.deg2rad([lon, lat])
    return np.array([(np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)),
                     (-np.sin(a), np.cos(a), 0),
                     (-np.cos(a)*np.sin(d), -np.sin(a)*np.sin(d), np.cos(d))])


def parse_pulsars(data):
    pulsars = []
    for entry in data.split("@"):
        fields = {p[0]: p[1:] for line in entry.splitlines() if len(p := line.split()) > 1}
        if "PSRJ" not in fields:
            continue
        def value(k):
            return float(fields[k][0])
        def error(k):
            return float(Decimal(fields[k][1]) * Decimal(10)**Decimal(fields[k][0]).as_tuple().exponent)
        try:
            ratio = value("P0")/value("P1") if "P0" in fields and "P1" in fields else -value("F0")/value("F1")
        except (KeyError, ZeroDivisionError):
            continue
        age = ratio / (2 * SECONDS_YEAR)
        if not 0 < age <= 3e6:
            continue
        if all(k in fields for k in ["RAJ","DECJ","PMRA","PMDEC"]):
            axes = basis(15*sexagesimal(fields["RAJ"][0]), sexagesimal(fields["DECJ"][0]))
            keys = ["PMRA","PMDEC"]
        elif all(k in fields for k in ["ELONG","ELAT","PMELONG","PMELAT"]):
            axes = basis(value("ELONG"), value("ELAT"))
            e = np.deg2rad(23.4392911)
            axes = axes @ np.array([[1,0,0],[0,np.cos(e),np.sin(e)],[0,-np.sin(e),np.cos(e)]])
            keys = ["PMELONG","PMELAT"]
        else:
            continue
        motion = np.array([value(k) for k in keys])
        if np.linalg.norm(motion) == 0 or not all(len(fields[k]) > 2 and fields[k][1].isdigit() for k in keys):
            continue
        errors = np.array([error(k) for k in keys])
        if (errors <= 0).any():
            continue
        pulsars.append((fields["PSRJ"][0], age, axes, motion, errors))
    return pulsars


def parse_shells(data):
    shells = []
    for line in data.splitlines():
        if 'HREF="snrs.G' not in line:
            continue
        p = re.sub("<[^>]*>", "", line).split()
        if not p[8].startswith(("S","C")):
            continue
        direction = basis(15*sexagesimal(":".join(p[2:5])), sexagesimal(":".join(p[5:7])))[0]
        diameters = [float(x) for x in p[7].replace("?", "").split("x")]
        radius = np.sqrt(diameters[0]*diameters[-1]) / 2
        shells.append(("G"+p[0]+p[1], direction, radius))
    return shells


def write_csv(name, fields, rows):
    with (HERE / name).open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(fields)
        w.writerows(rows)


def run(db, green):
    rng = np.random.default_rng(20261115)
    pulsars, shells = parse_pulsars(db), parse_shells(green)
    directions = np.array([s[1] for s in shells])
    radii = np.array([s[2] for s in shells])
    trials = 10000
    null_counts = np.zeros(trials, dtype=int)
    pairs = []
    uncertainty_rows = []
    for name, age, axes, motion, errors in pulsars:
        projection = directions @ axes.T
        idx = np.where(projection[:,0] > np.cos(np.deg2rad(8)))[0]
        if not len(idx):
            continue
        offset = -MAS_PER_RAD * projection[idx,1:] / projection[idx,:1]
        speed = np.linalg.norm(motion)
        sn = speed / np.linalg.norm(errors)
        r = radii[idx]*60000
        distance = np.linalg.norm(offset, axis=1)/r
        rel = offset @ motion / speed**2
        impact = np.linalg.norm(offset-rel[:,None]*motion, axis=1)/r
        lower, upper = max(0,0.02*age), min(190000,6*age)
        for j in np.where((rel>lower) & (rel<upper) & (impact<1.3))[0]:
            shell = shells[idx[j]][0]
            selected = bool(sn>=5 and distance[j]>1 and impact[j]<0.35)
            pairs.append([name, shell, round(age), round(rel[j]), round(sn,3),
                          round(distance[j],4), round(impact[j],4), selected])
            if selected:
                draw = rng.normal(motion, errors, size=(30000,2))
                t = (draw @ offset[j]) / np.sum(draw**2,axis=1)
                d = np.linalg.norm(offset[j]-t[:,None]*draw,axis=1)/r[j]
                q = np.quantile(t,[.16,.5,.84])
                uncertainty_rows.append([name,shell,*(round(v) for v in q),round(float(np.median(d)),4),
                                         round(float(np.mean((t>lower)&(t<upper)&(d<.35))),4)])
        # Whole-screen trials keep the sky, each motion amplitude and age fixed.
        if sn >= 5:
            phi = rng.uniform(0,2*np.pi,trials)
            unit = np.stack([np.cos(phi),np.sin(phi)],axis=1)
            for j in np.where(distance>1)[0]:
                along = unit @ offset[j]
                t = along/speed
                perp = np.abs(unit[:,0]*offset[j,1] - unit[:,1]*offset[j,0])/r[j]
                null_counts += (t>lower)&(t<upper)&(perp<.35)
    pairs.sort(key=lambda row:row[6])
    write_csv("pairs.csv",["pulsar","remnant","spin_age_yr","closest_lookback_yr","motion_snr",
                           "present_radii","impact_radii","exterior_shortlist"],pairs)
    write_csv("uncertainty.csv",["pulsar","remnant","t16_yr","t50_yr","t84_yr","median_impact_radii",
                                 "fraction_passing_time_and_impact"],uncertainty_rows)
    summary = dict(pulsars=len(pulsars),shell_or_composite_remnants=len(shells),loose_pairs=len(pairs),
                   outside_shortlist=len(uncertainty_rows),null_trials=trials,
                   shuffled_mean=float(null_counts.mean()),
                   shuffled_at_least_observed=int((null_counts>=len(uncertainty_rows)).sum()),
                   atnf_db_sha256=hashlib.sha256(db.encode()).hexdigest(),
                   green_html_sha256=hashlib.sha256(green.encode()).hexdigest())
    (HERE/"screen.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    cache = parser.parse_args().cache
    cache.mkdir(exist_ok=True,parents=True)
    for name, url in [("atnf.tgz",ATNF),("snr.html",GREEN)]:
        if not (cache/name).exists():
            urlretrieve(url,cache/name)
    with tarfile.open(cache/"atnf.tgz") as archive:
        db = archive.extractfile("psrcat_tar/psrcat.db").read().decode("utf-8")
    run(db,(cache/"snr.html").read_bytes().decode("utf-8"))
