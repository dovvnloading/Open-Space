import csv, concurrent.futures
from datetime import datetime, timedelta
from scipy.optimize import minimize_scalar
from scipy.interpolate import CubicHermiteSpline
import numpy as np
from net import CACHE, craft, vec
from screen import AU, date


def calc(c, s, x, vel, details=False):
    # Solve in emission time so the radial shell reaches the spacecraft radius at receipt.
    v = vel * 86400 / AU
    pos = s(x)
    rs = np.linalg.norm(pos, axis=-1)
    em = x
    for _ in range(10):
        co = c(em)
        rc = np.linalg.norm(co, axis=-1)
        em = x - (rs - rc) / v
    co = c(em)
    rc = np.linalg.norm(co, axis=-1)
    diff = np.linalg.norm(pos - co * (rs / rc)[..., None], axis=-1)
    if details:
        return diff, rs - rc, em, rc
    return np.where((rs - rc > 0) & (em >= c.x[0]) & (em <= c.x[-1]), diff, np.inf)


def loadarc(z):
    name = z["craft"]
    cmd, a, b = craft[name]
    sdata = vec(name, cmd, a, b, "1 d")
    s = CubicHermiteSpline(sdata[:, 0], sdata[:, 1:4], sdata[:, 4:7], extrapolate=False)
    t = datetime.fromisoformat(z["date"])
    dt = 33 if name == "ulysses" else 16
    # Account for drift in approximate periodic comet elements.
    w = 42
    start = (t - timedelta(days=w + dt)).date().isoformat()
    stop = (t + timedelta(days=w)).date().isoformat()
    cdata = vec(
        f"eph_{z['spkid']}_{start}_{stop}", f"DES={z['pdes']};CAP;", start, stop, "6 h"
    )
    c = CubicHermiteSpline(cdata[:, 0], cdata[:, 1:4], cdata[:, 4:7], extrapolate=True)
    mid = float(z["jd"])
    lo = max(mid - w, s.x[0], c.x[0] + dt)
    hi = min(mid + w - 1, s.x[-1], c.x[-1])
    return c, s, lo, hi


def run(z):
    c, s, lo, hi = loadarc(z)
    if lo >= hi:
        return []
    day = np.arange(lo, hi, 0.08)
    out = []
    for v in (250, 300, 350, 400, 450, 500, 600, 750, 900):
        ds = calc(c, s, day, v)
        best = np.argmin(ds)
        if not np.isfinite(ds[best]):
            continue
        # scipy bounded precision depends on size of abscissa; optimize relative to initial time.
        origin = day[best]
        val = minimize_scalar(
            lambda y: float(calc(c, s, origin + y, v)),
            bounds=(max(lo - origin, -0.12), min(hi - origin, 0.12)),
            method="bounded",
            options={"xatol": 1e-8},
        )
        x = origin + val.x
        dist, lag, em, rc = calc(c, s, x, v, True)
        if dist < 0.13 and lag > 0 and em >= c.x[0]:
            out.append(
                dict(
                    craft=z["craft"],
                    comet=z["comet"],
                    spkid=z["spkid"],
                    guess=z["date"],
                    speed=v,
                    date=date(x),
                    miss_km=round(dist * AU),
                    radial_au=round(lag, 5),
                    rh=round(rc, 5),
                    jd=round(x, 8),
                    emit_jd=round(em, 8),
                )
            )
    return out


if __name__ == "__main__":
    rows = list(csv.DictReader(open(CACHE / "screen.csv")))
    # This is a triage subset of the coarse search; archival orbit validity requires separate checks.
    work = [
        z
        for z in rows
        if not z["comet"].startswith(("D/", "5D/"))
        and (float(z["miss_au"]) < 0.03 or z["craft"] in ("wind", "solo", "stereo-a"))
        and (z["craft"] != "earth" or z["date"] > "2026-11-15")
    ]
    out = []
    receipts = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        jobs = {pool.submit(run, z): z for z in work}
        for future in concurrent.futures.as_completed(jobs):
            z = jobs[future]
            try:
                result = future.result()
                out.extend(result)
                status = "retained" if result else "no retained model"
            except (OSError, ValueError) as err:
                status = "unavailable: " + str(err).replace("\n", " ")[:260]
            receipts.append(
                dict(craft=z["craft"], comet=z["comet"], guess=z["date"], status=status)
            )
    out.sort(key=lambda z: (z["miss_km"], z["craft"], z["jd"], z["speed"]))
    for file, rows in [
        ("refined.csv", out),
        (
            "refinement-status.csv",
            sorted(receipts, key=lambda x: (x["craft"], x["guess"], x["comet"])),
        ),
    ]:
        if not rows:
            continue
        with open(CACHE / file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    print(f"{len(receipts)} requests, {len(out)} retained models.")
