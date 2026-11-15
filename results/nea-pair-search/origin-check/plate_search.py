"""Compare the predicted November 1988 positions to STScI DSS and GAVO plate centers."""
import argparse
import csv
from datetime import datetime, timedelta
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from encounters import C, orbit_trajectory, planet_ephemeris

DSS = "https://outerspace.stsci.edu/download/attachments/176435492/SurveyPlateListLatest.csv?version=1&modificationDate=1718908770237&api=v2"
QUERY = "SELECT * FROM wfpdb.main WHERE epoch BETWEEN 1988.845 AND 1988.884"
WFP = "https://dc.g-vo.org/tap/sync?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=QUERY))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("kernel", type=Path)
    cli.add_argument("gravity_fit", type=Path)
    cli.add_argument("acceleration_fits", type=Path)
    cli.add_argument("--snapshots", type=Path, help="use local dss.csv and platepos.csv")
    args = cli.parse_args()
    epoch = float(json.loads(args.gravity_fit.read_text())["epoch_jd_tdb"])
    state = json.loads(args.acceleration_fits.read_text())["all/A1"]["state"]
    planets = planet_ephemeris(args.kernel, epoch, .5)
    old_arc, _ = orbit_trajectory(state[:6], state[6], state[7], planets, epoch)
    obliquity = np.deg2rad(23.439291111)
    cos_e, sin_e = np.cos(obliquity), np.sin(obliquity)

    def prediction(time):
        # Time differences of order one minute between plate UT and TDB, plus
        # ignored topocentric parallax, are small compared to degree-scale misses.
        t = (time-datetime(2000,1,1,12)).total_seconds()/86400+2451545-epoch
        earth = planets(t)[2]
        rough = old_arc[2](t)[:3]
        r = old_arc[2](t-np.linalg.norm(rough-earth)/C)[:3]
        relative = r-earth
        delta, radius = np.linalg.norm(relative), np.linalg.norm(r)
        phase = np.arccos(np.clip(r@relative/radius/delta,-1,1))
        tan = np.tan(phase/2)
        phi = .85*np.exp(-3.33*tan**.63)+.15*np.exp(-1.87*tan**1.22)
        mag = 26+5*np.log10(radius*delta)-2.5*np.log10(phi)
        x,y,z = relative/delta
        unit = np.array([x,cos_e*y-sin_e*z,sin_e*y+cos_e*z])
        return unit, mag

    report = {}
    for survey, url, file in [("DSS", DSS, "dss.csv"), ("WFPDB", WFP, "platepos.csv")]:
        raw = (args.snapshots/file).read_bytes() if args.snapshots else urlopen(url,timeout=90).read()
        catalog = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        rows, selected = [], 0
        for r in catalog:
            if survey == "DSS":
                if not r["UTobservation"].startswith("1988-11"): continue
                selected += 1
                center = (r["RightAscension"],r["Declination"])
                date = datetime.fromisoformat(r["UTobservation"])
                ident, exptime = r["plateLabel"], None
            else:
                selected += 1
                center = (r["raj2000"],r["dej2000"])
                date = datetime(1988,1,1)+timedelta(days=(float(r["epoch"])-1988)*366)
                ident = r["wfpdbid"]
                exptime = float(r["exptime"]) if r["exptime"] else None
                date += timedelta(seconds=(exptime or 0)/2)
            if not all(center): continue
            a,b = np.deg2rad([float(value) for value in center])
            field = [np.cos(a)*np.cos(b), np.sin(a)*np.cos(b), np.sin(b)]
            u,mag = prediction(date)
            angle = np.rad2deg(np.arccos(np.clip(u@field,-1,1)))
            rows.append(dict(plate=ident, coordinate_evaluation_utc=date.isoformat(timespec="seconds"),
                             exposure_seconds=exptime, field_center_deg=list(map(float,center)),
                             asteroid_ra_dec_deg=[round(np.rad2deg(np.arctan2(u[1],u[0]))%360,4),
                                                  round(np.rad2deg(np.arcsin(u[2])),4)],
                             center_distance_deg=round(angle,4), illustrative_V=round(mag,2)))
        report[survey] = {"source":url,"snapshot_sha256":hashlib.sha256(raw).hexdigest(),
                          "rows_in_time_window":selected,"rows_with_centers":len(rows),
                          "nearest_centers":sorted(rows,key=lambda r:r["center_distance_deg"])[:3]}
    first = report["WFPDB"]["nearest_centers"][0]
    # Actual plate size is from its WFPDB row. Focal length is documented by
    # the primary Asiago study linked in the receipt. Diagonal uses ideal rectilinear projection.
    report["nearest_plate_geometry"] = {
        "id":first["plate"], "side_cm":12, "focal_length_cm":100,
        "half_diagonal_deg":round(np.rad2deg(np.arctan(np.sqrt(2)*6/100)),4),
        "source":"https://openaccess.inaf.it/bitstream/20.500.12386/27630/1/aa32131-17.pdf",
        "scope":"No verified third-apparition observation. Catalog centers alone cannot establish a detection. Unknown centers, plate depths and unindexed archives remain untested.",
    }
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
