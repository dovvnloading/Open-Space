"""Replay a bounded search for Gaia BH3 stellar alignments. Standard library only."""
import csv
import hashlib
import io
import json
import math
import random
from pathlib import Path

BH_ID = "4318465066420528000"
RA, DEC = 294.8278502301, 14.9309190869
MU = (-28.317, -155.221)
PLX, MASS = 1.6933, 32.7
MAS_PER_DEG = 3600000
COS_DEC = math.cos(math.radians(DEC))
AXES = ("ra", "dec", "pmra", "pmdec")


def offsets(row, epoch=2016.0):
    x = (float(row["ra"]) - RA) * COS_DEC * MAS_PER_DEG
    y = (float(row["dec"]) - DEC) * MAS_PER_DEG
    return x - MU[0]*(epoch-2017.5), y - MU[1]*(epoch-2017.5)


def closest(x, y, vx, vy, epoch):
    dt = -(x*vx + y*vy)/(vx*vx + vy*vy)
    return epoch + dt, math.hypot(x+vx*dt, y+vy*dt)


def predict(row, dx=0, dy=0, dpmx=0, dpmy=0, epoch=2016.0,
            lens_mu=MU, lens_pos=(0,0)):
    x, y = offsets(row, epoch)
    x += dx - lens_pos[0] - (lens_mu[0]-MU[0])*(epoch-2017.5)
    y += dy - lens_pos[1] - (lens_mu[1]-MU[1])*(epoch-2017.5)
    return closest(x, y, float(row["pmra"]) + dpmx - lens_mu[0],
                   float(row["pmdec"]) + dpmy - lens_mu[1], epoch)


def deflection(b, mass, lens_plx, distance):
    theta2 = 8.144 * mass * (lens_plx-1000/distance)
    return theta2*b/(b*b+2*theta2)


def cholesky(row):
    """Lower factor of Gaia's marginal position / proper-motion correlation matrix."""
    n = len(AXES)
    lower = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            v = 1 if i == j else float(row[f"{AXES[j]}_{AXES[i]}_corr"])
            v -= sum(lower[i][k]*lower[j][k] for k in range(j))
            lower[i][j] = math.sqrt(v) if i == j else v/lower[j][j]
    return lower


def interval(values):
    values.sort()
    return [round(values[int(q*(len(values)-1))], 4) for q in (.025,.5,.975)]


def analyze(field_text, data):
    field = list(csv.DictReader(io.StringIO(field_text)))
    eligible = [r for r in field if r["source_id"] != BH_ID and r["pmra"] and
                r["pmdec"] and r["ruwe"] and float(r["ruwe"]) < 1.4]
    in_window = [(predict(r), r) for r in eligible]
    in_window = [(p,r) for p,r in in_window if 1900 <= p[0] <= 2125]
    close = sorted([(p,r) for p,r in in_window if p[1]<1000], key=lambda item:item[0])
    lens_g = float(next(r["phot_g_mean_mag"] for r in field if r["source_id"] == BH_ID))
    output = {"field_rows":len(field), "eligible":len(eligible),
              "closest_approach_in_window":len(in_window), "within_one_arcsec":len(close),
              "field_sha256":hashlib.sha256(field_text.encode()).hexdigest(),"events":[]}
    rng = random.Random(1179)
    for (year, b), row in close:
        sid = row["source_id"]
        extra = next(r for r in data["extra"] if r["source_id"] == sid)
        distance = data["distances"][sid]
        rgeo = float(distance["rgeo"])
        factor = cholesky(extra)
        epochs, impacts, shifts = [], [], []
        for _ in range(15000):
            normals = [rng.gauss(0,1) for axis in AXES]
            delta = [float(extra[axis+"_error"])*sum(scale*g for scale,g in zip(coeff,normals))
                     for axis, coeff in zip(AXES,factor)]
            mu = (rng.gauss(MU[0],.067),rng.gauss(MU[1],.111))
            pos = (rng.gauss(0,.054),rng.gauss(0,.074))
            t, beta = predict(row,*delta,lens_mu=mu,lens_pos=pos)
            epochs.append(t)
            impacts.append(beta)
            shifts.append(deflection(beta,rng.gauss(MASS,.82),rng.gauss(PLX,.0164),rgeo))
        def ps_gap(p):
            dt = 2000+(float(p["epochMean"])-51544.5)/365.25-2016
            dx = (float(p["raMean"])-float(row["ra"]))*COS_DEC*MAS_PER_DEG
            dy = (float(p["decMean"])-float(row["dec"]))*MAS_PER_DEG
            return math.hypot(dx-dt*float(row["pmra"]),dy-dt*float(row["pmdec"]))
        ps = min(data["ps1"],key=ps_gap)
        dr2_id = "4318465066396437632" if sid == "4318465066396437504" else sid
        dr2 = next(r for r in data["dr2"] if r["source_id"]==dr2_id)
        y2,b2 = predict(dr2,epoch=2015.5)
        raw_px = float(row["parallax"])
        delta_g = float(row["phot_g_mean_mag"]) - lens_g
        shift = deflection(b,MASS,PLX,rgeo)
        # Falsifiable relative position prediction at an epoch far from lensing peaks.
        x, y = offsets(row)
        dt = 2035-2016
        pred = [x + dt*(float(row["pmra"])-MU[0]), y+dt*(float(row["pmdec"])-MU[1])]
        output["events"].append({
            "source":sid,"year":round(year,4),"separation_mas":round(b,3),
            "bh_shift_mas":round(shift,4), "year_simulation":interval(epochs),
            "separation_mas_simulation":interval(impacts), "shift_mas_simulation_fixed_distance":interval(shifts),
            "shift_with_raw_parallax_mas":round(deflection(b,MASS,PLX,1000/raw_px),4),
            "shift_at_distance_q16_q84_mas":[round(deflection(b,MASS,PLX,float(distance[k])),4)
                                             for k in ("b_rgeo","B_rgeo")],
            "dr2_year_and_impact_mas":[round(y2,3),round(b2,2)],
            "ps1_id":ps["objID"], "ps1_detections":int(ps["nDetections"]),
            "ps1_position_residual_mas":round(ps_gap(ps),2),
            "g_contrast_mag":round(delta_g,3),
            "2035_barycentric_east_north_mas":[round(p,2) for p in pred]})
    return output


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    data = json.loads((here/"followup.json").read_text())
    result = analyze((here/"field.csv").read_text(),data)
    expected = json.loads((here/"measurements.json").read_text())
    if result != expected:
        raise SystemExit("Calculated values differ from the deposited results")
    print(json.dumps(result,indent=2))
