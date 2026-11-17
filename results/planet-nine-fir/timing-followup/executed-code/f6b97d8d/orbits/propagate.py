"""Conditional Phan et al. (2025) orbit grid. Units: au, days, J2000 equatorial."""
import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import erfa
import numpy as np

# DE440 Sun plus planetary systems, expressed in au^3/day^2.
GM_KM = [132712440041.279419, 22031.868551, 324858.592000, 398600.435507,
         4902.800118, 42828.375816, 126712764.100000, 37940584.841800,
         5794556.400000, 6836527.100580, 975.500000]
MU = sum(GM_KM) * 86400**2 / 149597870.7**3
C = 173.144632674240       # light speed in au/day
OBLIQUITY = np.radians(23.43928)
ARCSEC = 180 / np.pi * 3600


def epoch(s):
    t = datetime.fromisoformat(s.replace("Z", "+00:00"))
    if t.tzinfo is None:
        t = t.replace(tzinfo=timezone.utc)
    t = t.astimezone(timezone.utc)
    utc = erfa.dtf2d("UTC", t.year, t.month, t.day, t.hour, t.minute,
                     t.second + t.microsecond/1e6)
    tt1, tt2 = erfa.taitt(*erfa.utctai(*utc))
    return (tt1 - 2451545.0) + tt2


def observers(t):
    """Geocentre heliocentric and SSB positions (au, BCRS). TT approximates TDB."""
    epochs, inv = np.unique(np.atleast_1d(t), return_inverse=True)
    h, b = erfa.epv00(2451545., epochs)
    shape = np.shape(t) + (3,)
    return h["p"][inv].reshape(shape), b["p"][inv].reshape(shape)


def earth(t):
    return observers(t)[1]


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def norm(a):
    return np.linalg.norm(a, axis=-1)


def ray(observer, n, r):
    """Position on positive observed ray at the specified heliocentric distance."""
    dot = np.sum(observer*n,axis=-1)
    dist = -dot + np.sqrt(dot**2 + r**2 - norm(observer)**2)
    return observer + dist[..., None] * n, dist


def kepler(rvec, vel, dt):
    """Universal Kepler equation with short-arc Stumpff series, stable at E=0."""
    r = norm(rvec)
    rv = np.sum(rvec * vel, axis=-1) / np.sqrt(MU)
    alpha = 2 / r - np.sum(vel**2, axis=-1) / MU
    h = np.sqrt(MU) * dt
    chi = h / r
    for _ in range(4):
        z = alpha * chi**2
        if np.max(np.abs(z)) > 0.05:
            raise ValueError("Propagation beyond validated short-arc domain")
        cs = 1/2 + z * (-1/24 + z * (1/720 + z * (-1/40320 + z/3628800)))
        ss = 1/6 + z * (-1/120 + z * (1/5040 + z * (-1/362880 + z/39916800)))
        f = rv * chi**2 * cs + (1 - alpha*r) * chi**3 * ss + r*chi - h
        df = rv * chi * (1-z*ss) + (1-alpha*r) * chi**2 * cs + r
        chi -= f / df
    fg_f = 1 - chi**2 * cs / r
    fg_g = dt - chi**3 * ss / np.sqrt(MU)
    return fg_f[..., None]*rvec + fg_g[..., None]*vel


def family(inputs, n_time=9, n_radius=21, n_rate=41, coadd_endpoint_fraction=0., parameters=None, angular_offsets=None):
    a, b = [epoch(x) for x in inputs["iras_coadd_interval_utc"]]
    tb = epoch(inputs["akari_midpoint_utc"])
    # Bracket endpoint range changes slightly beyond the local escape speed.
    if parameters is None:
        t1, r2, fraction = [q.ravel() for q in np.meshgrid(
            np.linspace(a, b, n_time), np.linspace(*inputs["r_akari_au"], n_radius),
            np.linspace(-1.05, 1.05, n_rate), indexing="ij")]
    else:
        # Unit-cube draws parameterize the same support without coincident grid lines.
        u_time,u_radius,u_rate=np.asarray(parameters).T
        t1=a+(b-a)*u_time
        low,high=inputs["r_akari_au"]
        r2=low+(high-low)*u_radius
        fraction=2.1*u_rate-1.05
    n1, n2 = [unit(*inputs[k]) for k in ("iras_radec_deg", "akari_radec_deg")]
    if angular_offsets is not None:
        directions=[]
        for normal,key,offsets in zip((n1,n2),("iras_radec_deg","akari_radec_deg"),angular_offsets):
            ra=np.radians(inputs[key][0])
            east=np.array([-np.sin(ra),np.cos(ra),0.])
            north=np.cross(normal,east)
            vec=normal+offsets[:,0,None]*east/ARCSEC+offsets[:,1,None]*north/ARCSEC
            directions.append(vec/norm(vec)[:,None])
        n1,n2=directions
    r1 = r2 - fraction * np.sqrt(2*MU/r2) * (tb-t1)
    h1, b1 = observers(t1)
    if coadd_endpoint_fraction:
        # At each fixed mean scan date, compare a single epoch with a coadd
        # of the first and last epochs, weighted to that same mean date.
        x = ((t1-a)/(b-a))[:, None]
        endpoint = (1-x)*earth(a) + x*earth(b)
        b_adjust = coadd_endpoint_fraction*(endpoint-b1)
        h1 += b_adjust
        b1 += b_adjust
    h2, b2 = observers(tb)
    pos1, dist1 = ray(h1, n1, r1)
    pos2, dist2 = ray(h2, n2, r2)
    pos1 += b1-h1
    pos2 += b2-h2
    dt = t1 - dist1/C - (tb - dist2/C)
    v = (pos1 - pos2) / dt[:, None]
    for _ in range(4):
        v += (pos1 - kepler(pos2, v, dt)) / dt[:, None]
    residual = norm(kepler(pos2, v, dt) - pos1) / dist1 * ARCSEC
    energy = np.sum(v**2, axis=-1)/2 - MU/norm(pos2)
    keep = energy < 0
    angmom = np.cross(pos2[keep], v[keep])
    ecl_z = np.array([0., -np.sin(OBLIQUITY), np.cos(OBLIQUITY)])
    inc = np.degrees(np.arccos(np.clip(angmom @ ecl_z / norm(angmom), -1, 1)))
    info = dict(
        tested=len(r2), bound=int(keep.sum()),
        inclination_ecliptic_deg=[float(inc.min()), float(inc.max())],
        max_endpoint_residual_arcsec=float(residual[keep].max()),
        boundary_rate_reached=bool(np.any(np.abs(fraction[keep]) >= 1.05 - 1e-12)),
        grid=[n_time, n_radius, n_rate] if parameters is None else None, coadd_endpoint_fraction=coadd_endpoint_fraction)
    return pos2[keep], v[keep], tb - dist2[keep]/C, info


def direction(states, t, observer=None):
    r, v, t0, _ = states
    obs = earth(t) if observer is None else observer
    dt = t - t0 - norm(r - obs)/C
    for _ in range(2):
        q = kepler(r, v, dt) - obs
        dt = t - t0 - norm(q)/C
    q = kepler(r, v, dt) - obs
    return q / norm(q)[:, None]


def sky(states, t, observer=None):
    q = direction(states, t, observer)
    ra = np.degrees(np.arctan2(q[:, 1], q[:, 0])) % 360
    dec = np.degrees(np.arcsin(q[:, 2]))
    return ra, dec


def default_epochs():
    return sorted({f"{y}-{m:02d}-01" for y in range(2006, 2027) for m in range(1, 13)}
                  | {"2026-11-17", "2026-12-31"})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", type=Path, default=Path(__file__).with_name("orbit_input.json"))
    p.add_argument("--epochs", type=Path, help="One ISO UTC exposure midpoint per line")
    p.add_argument("--out", type=Path, default=Path.cwd())
    p.add_argument("--grid", type=int, nargs=3, default=(9, 21, 41), metavar=("TIMES","RADII","RATES"))
    p.add_argument("--coadd-endpoint-fraction", type=float, choices=[0., .5, 1.], default=0.)
    p.add_argument("--tracks", type=Path, help="Optional gzip CSV of every track at the requested epochs")
    args = p.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    states = family(json.loads(args.inputs.read_text()), *args.grid, args.coadd_endpoint_fraction)
    dates = args.epochs.read_text().splitlines() if args.epochs else default_epochs()
    dates = [d.strip() for d in dates if d.strip()]
    if args.tracks:
        import gzip
        tracks = gzip.open(args.tracks, "wt", newline="")
        tw = csv.writer(tracks)
        tw.writerow(["utc", "orbit_index", "ra_deg", "dec_deg"])
    with (args.out / "ephemeris.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["utc", "ra_min_deg", "ra_max_deg", "dec_min_deg", "dec_max_deg",
                    "speed_min_as_per_hour", "speed_max_as_per_hour"])
        for d in dates:
            t = epoch(d)
            q = direction(states, t)
            ra = np.degrees(np.arctan2(q[:, 1], q[:, 0])) % 360
            dec = np.degrees(np.arcsin(q[:, 2]))
            later = direction(states, t + 1/24)
            speed = 2 * np.arcsin(norm(later-q)/2) * ARCSEC
            w.writerow([d, *[f"{x:.6f}" for x in (ra.min(), ra.max(), dec.min(), dec.max(), speed.min(), speed.max())]])
            if args.tracks:
                tw.writerows((d, j, f"{a:.8f}", f"{b:.8f}") for j, (a, b) in enumerate(zip(ra, dec)))
    if args.tracks:
        tracks.close()
    (args.out / "orbit_summary.json").write_text(json.dumps(states[-1], indent=2) + "\n")


if __name__ == "__main__":
    main()
