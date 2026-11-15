import csv, json, math
import numpy as np
from net import CACHE, craft, vec
from datetime import datetime, timedelta

AU = 149597870.7
K = 0.01720209895


def date(j):
    return (datetime(2000, 1, 1, 12) + timedelta(days=float(j) - 2451545)).isoformat(
        timespec="minutes"
    )


def states(q, e, inc, node, arg):
    v = np.arange(-math.pi, math.pi, 0.009)[1:]
    c = np.cos(v)
    keep = 1 + e * c > q * (1 + e) / 2.8
    v = v[keep]
    r = q * (1 + e) / (1 + e * np.cos(v))
    if abs(1 - e) < 1e-6:
        d = np.tan(v / 2)
        t = (2 * q**3) ** 0.5 / K * (d + d**3 / 3)
    elif e < 1:
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(v / 2))
        t = (q / (1 - e)) ** 1.5 / K * (E - e * np.sin(E))
    else:
        H = 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(v / 2))
        t = (q / (e - 1)) ** 1.5 / K * (e * np.sinh(H) - H)
    ci, co = np.cos([inc, node])
    si, so = np.sin([inc, node])
    u = v + arg
    xyz = np.column_stack(
        [
            co * np.cos(u) - so * np.sin(u) * ci,
            so * np.cos(u) + co * np.sin(u) * ci,
            np.sin(u) * si,
        ]
    )
    return t, r, xyz


def main():
    arcs = {}
    for name, (cmd, a, b) in craft.items():
        d = vec(name, cmd, a, b, "1 d")
        arcs[name] = d
    j = json.load(open(CACHE / "sbdb.json"))
    out = []
    n = 0
    for row in j["data"]:
        id, pdes, fullname, epoch, q, e, i, om, w, tp, m1, k1, period = row
        if not all([q, e, i, om, w, tp]) or not 0.08 < float(q) < 2.2:
            continue
        if fullname.strip().endswith("(SOHO)") and not m1:
            continue
        q, e, tp = map(float, [q, e, tp])
        per = 2 * np.pi * (q / (1 - e)) ** 1.5 / K if e < 1 else 1e15
        t, r, unit = states(q, e, *np.deg2rad([float(i), float(om), float(w)]))
        tmin, tmax = 2448150, 2462900
        k0 = math.ceil((tmin - tp - t[-1]) / per)
        k1 = math.floor((tmax - tp - t[0]) / per)
        for k in range(k0, k1 + 1):
            peri = tp + k * per
            launch = peri + t
            n += 1
            for name, d in arcs.items():
                if launch[-1] < d[0, 0] - 40 or launch[0] > d[-1, 0]:
                    continue
                V = 400 * 86400 / AU
                dt0 = launch - d[0, 0]

                def pos(idx):
                    idx = np.clip(idx, 0, len(d) - 1.000001)
                    ix = idx.astype(int)
                    f = (idx - ix)[:, None]
                    return d[ix, 1:4] * (1 - f) + d[ix + 1, 1:4] * f

                s = pos(dt0)
                rad = np.linalg.norm(s, axis=1)
                for _ in range(2):
                    lag = (rad - r) / V
                    s = pos(dt0 + lag)
                    rad = np.linalg.norm(s, axis=1)
                miss = np.linalg.norm(s - unit * rad[:, None], axis=1)
                miss[
                    (lag < 0.02)
                    | (lag > 30)
                    | (dt0 + lag < 0)
                    | (dt0 + lag > len(d) - 1)
                ] = np.inf
                ix = np.argmin(miss)
                if miss[ix] < 0.048:
                    out.append(
                        dict(
                            craft=name,
                            spkid=id,
                            comet=fullname.strip(),
                            pdes=pdes,
                            m1=m1,
                            peri=round(peri, 4),
                            jd=round(launch[ix] + lag[ix], 5),
                            date=date(launch[ix] + lag[ix]),
                            miss_au=round(miss[ix], 6),
                            dist_au=round(rad[ix] - r[ix], 4),
                            rh=round(r[ix], 4),
                        )
                    )
    out.sort(key=lambda z: z["miss_au"])
    with open(CACHE / "screen.csv", "w") as f:
        w = csv.DictWriter(f, fieldnames=out[0])
        w.writeheader()
        w.writerows(out)
    print("arcs, apparitions, leads", len(arcs), n, len(out))
    for z in out[:25]:
        print(z["craft"], z["comet"], z["date"], z["miss_au"], z["dist_au"], z["m1"])


if __name__ == "__main__":
    main()
