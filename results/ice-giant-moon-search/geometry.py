"""Spacecraft Horizons queries and observer-projected ice-giant orbit templates."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

AU = 149597870.7
ARCSEC = 180 * 3600 / np.pi
C = 299792.458
CENTER = {"JWST": "500@-170", "HST": "500@-48"}
PLANET = {"Uranus": 799, "Neptune": 899}
MOONS = {
    "Uranus": [*range(701, 728), 75051, 75052],
    "Neptune": [*range(801, 815), 85051, 85052],
}
# GM [km3/s2], reference radius [km], J2, J4. J2,J4 are scaled to R below.
FIELD = {"Uranus": (5793951.3, 25559., 3.5107e-3, -34.17e-6),
         "Neptune": (6835100., 25225., 3.40843e-3, -33.40e-6)}
# Conservative radial envelopes including the eccentric epsilon ring and full arcs.
RINGS = {"Uranus": [(26800, 51650), (65400, 69200), (89200, 106200)],
         "Neptune": [(41000, 43000), (53150, 57300), (61850, 62050), (62850, 63020)]}


def query(body, telescope, jd, cache):
    """Return ICRF RA, Dec [deg], angular diameter [arcsec], r and delta [AU].

    Retain the complete API response. Split epoch lists to keep GET requests small.
    One request at a time, as required by SSD API fair-use terms.
    """
    unique, reverse = np.unique(np.asarray(jd).ravel().round(9), return_inverse=True)
    chunks = []
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    for start in range(0, len(unique), 60):
        times = unique[start:start+60]
        # Horizons uses a quoted multi-epoch list separated by commas.
        params = dict(format="json", COMMAND=f"'{body}'", CENTER=f"'{CENTER[telescope]}'",
                      MAKE_EPHEM="'YES'", OBJ_DATA="'NO'", EPHEM_TYPE="'OBSERVER'",
                      TLIST="'" + ",".join(f"{j:.9f}" for j in times) + "'",
                      QUANTITIES="'1,13,19,20'", CAL_FORMAT="'JD'", ANG_FORMAT="'DEG'",
                      EXTRA_PREC="'YES'", CSV_FORMAT="'YES'", TIME_TYPE="'UT'")
        url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
        key = hashlib.sha256(url.encode()).hexdigest()
        target = cache / f"{body}-{key}.json"
        if target.exists():
            raw = target.read_bytes()
        else:
            with urlopen(url, timeout=100) as response:
                raw = response.read()
        payload = json.loads(raw)
        result = payload.get("result", "")
        if "$$SOE" not in result or "$$EOE" not in result:
            raise ValueError(f"Horizons {body} {telescope}: " + payload.get("error", result)[-600:])
        header, rows = result.split("$$SOE")
        if f"({body})" not in header.split("Center body name:")[0]:
            raise ValueError(f"unexpected Horizons target for {body}")
        data = list(csv.reader(io.StringIO(rows.split("$$EOE")[0].strip())))
        # Blank daylight/visibility indicators occupy positions 1 and 2.
        def number(v):
            return float(v) if v.strip() not in {"n.a.", ""} else np.nan
        values = np.array([[number(row[k]) for k in (0, 3, 4, 5, 6, 8)] for row in data])
        if values.shape != (len(times), 6) or np.max(np.abs(values[:, 0] - times)) > 1e-7:
            raise ValueError(f"incomplete or reordered Horizons response for {body}")
        if not np.isfinite(values[:, [1, 2, 4, 5]]).all():
            raise ValueError(f"missing Horizons coordinates or distances for {body}")
        if not target.exists():
            temp = target.with_suffix(".tmp")
            temp.write_bytes(raw)
            temp.replace(target)
        chunks.append(values[:, 1:])
    return np.concatenate(chunks)[reverse].reshape(np.shape(jd) + (5,))


def direction(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)), axis=-1)


def tangent_basis(ra, dec):
    los = direction(ra, dec)
    east = direction(np.asarray(ra) + 90, np.zeros_like(ra))
    return los, east, np.cross(los, east)


def pole(planet, jd):
    if planet == "Uranus":
        # Retrograde spin: regular moons' angular momentum is opposite IAU north.
        return direction(77.311, 15.175)
    n = np.deg2rad(357.85 + 52.316 * (jd - 2451545.) / 36525.)
    return direction(299.36 + 0.70*np.sin(n), 43.46 - 0.51*np.cos(n))


def basis(planet, jd):
    h = pole(planet, jd)
    u = np.cross([0., 0., 1.], h)
    u /= np.linalg.norm(u)
    return u, np.cross(h, u), h


def sky_xy(ra, dec, center_ra, center_dec, scale, center_xy):
    """Exact forward gnomonic projection. East left, north up (array y increases north)."""
    los, east, north = tangent_basis(center_ra, center_dec)
    vec = direction(ra, dec)
    denom = (vec * los).sum(axis=-1)
    x = center_xy[0] - (vec*east).sum(axis=-1) / denom * ARCSEC/scale
    y = center_xy[1] + (vec*north).sum(axis=-1) / denom * ARCSEC/scale
    return np.stack([np.where(denom > 0, x, np.nan), np.where(denom > 0, y, np.nan)], axis=-1)


def ring_mask(planet, ra, dec, delta, jd, scale, center, shape, pad):
    """Rasterize projected annuli in the image and dilate in sky pixels."""
    from scipy.ndimage import distance_transform_edt
    los, east, north = tangent_basis(ra, dec)
    h = pole(planet, jd)
    yy, xx = np.indices(shape)
    km_px = delta*AU*scale/ARCSEC
    e, n = -(xx-center[0])*km_px, (yy-center[1])*km_px
    # For nearly edge-on rings, sampled projected boundaries still survive rasterization.
    mask = np.zeros(shape, dtype=bool)
    h_los = h @ los
    if abs(h_los) > 1e-8:
        depth = -(e*(h @ east) + n*(h @ north))/h_los
        radius = np.sqrt(e*e + n*n + depth*depth)
        for lo, hi in RINGS[planet]:
            mask |= (radius >= lo) & (radius <= hi)
    u, v, _ = basis(planet, jd)
    for radius in np.ravel(RINGS[planet]):
        angle = np.linspace(0, 2*np.pi, max(256, int(4*np.pi*radius/km_px)+1))
        pos = radius*(np.cos(angle)[:, None]*u + np.sin(angle)[:, None]*v)
        x = np.rint(center[0] - pos@east/km_px).astype(int)
        y = np.rint(center[1] + pos@north/km_px).astype(int)
        inside = (x>=0) & (x<shape[1]) & (y>=0) & (y<shape[0])
        mask[y[inside], x[inside]] = True
    return (distance_transform_edt(~mask) <= pad) if mask.any() else mask


def frequency(planet, radius):
    gm, ref, j2, j4 = FIELD[planet]
    q = (ref/radius)**2
    return np.sqrt(gm/radius**3 * (1 + 1.5*j2*q - 1.875*j4*q*q))


def integrate(planet, a, e, inc, node, peri, mean, times, frame):
    """J2+J4 integration of explicit initial osculating Kepler elements."""
    from scipy.integrate import solve_ivp
    gm, ref, j2, j4 = FIELD[planet]
    u, v, h = frame
    n = np.cos(node)*u + np.sin(node)*v
    m = np.cross(h, n)*np.cos(inc) + h*np.sin(inc)
    p = n*np.cos(peri) + m*np.sin(peri)
    q = -n*np.sin(peri) + m*np.cos(peri)
    E = mean
    for _ in range(20):
        dE = (E-e*np.sin(E)-mean)/(1-e*np.cos(E))
        E -= dE
        if abs(dE) < 1e-14:
            break
    pos = a*((np.cos(E)-e)*p + np.sqrt(1-e*e)*np.sin(E)*q)
    vel = np.sqrt(gm/a)/(1-e*np.cos(E))*(-np.sin(E)*p + np.sqrt(1-e*e)*np.cos(E)*q)

    def deriv(t, state):
        r = state[:3]
        d = np.linalg.norm(r)
        s = r@h/d
        acc = -gm*r/d**3
        for k, j, polynomial, derivative in (
            (2, j2, (3*s*s-1)/2, 3*s),
            (4, j4, (35*s**4-30*s*s+3)/8, (140*s**3-60*s)/8),
        ):
            acc += gm/d**2*j*(ref/d)**k * (((k+1)*polynomial+s*derivative)*r/d-derivative*h)
        return np.r_[state[3:], acc]

    coords = np.empty(times.shape+(3,))
    coords[times==0] = pos
    for sign in (-1,1):
        mask = times*sign>0
        if not np.any(mask):
            continue
        unique, reverse = np.unique(times[mask]*sign, return_inverse=True)
        solution = solve_ivp(deriv, (0, sign*unique[-1]), np.r_[pos, vel],
                             method="DOP853", rtol=1e-11, atol=1e-10, t_eval=unique*sign)
        if not solution.success:
            raise RuntimeError(solution.message)
        coords[mask] = solution.y[:3].T[reverse]
    return coords


def generate(cube, output, elements=None, grid=None, model=None, epoch=None):
    from astropy.time import Time
    with np.load(cube, allow_pickle=False) as f:
        planet = str(f["planet"])
        ra, dec, delta = [f[k] for k in ("path_ra", "path_dec", "path_delta")]
        jd = Time(f["path_mjd"], format="mjd", scale="utc").tdb.jd - delta*AU/C/86400.
        ids, center, scale = f["exposure_id"], f["center_xy"], float(f["scale_arcsec"])
    epoch = float(np.median(jd) if epoch is None else epoch)
    times = (jd - epoch)*86400.
    los, east, north = tangent_basis(ra, dec)
    frame = basis(planet, epoch)
    if elements:
        rows = np.genfromtxt(elements, delimiter=",", names=True, dtype=float, ndmin=1)
        trials = np.column_stack([rows[k] for k in ("a_km", "e", "i_deg", "node_deg", "peri_deg", "mean_deg")])
    else:
        lo, hi, step, phase_step = grid
        if lo <= 0 or hi < lo or step <= 0 or not 0 < phase_step <= 360:
            raise ValueError("invalid radial or angular grid")
        a, phase = np.meshgrid(np.arange(lo, hi+step*1e-9, step), np.arange(0,360,phase_step))
        trials = np.zeros((a.size, 6))
        trials[:,0], trials[:,-1] = a.ravel(), phase.ravel()
    model = model or ("osculating" if elements else "circular")
    x, y = [np.empty((len(trials),)+jd.shape, dtype=np.float32) for _ in range(2)]
    for index, (a,e,inc,node,peri,mean) in enumerate(trials):
        if not np.isfinite([a,e,inc,node,peri,mean]).all() or a*(1-e) <= (25559 if planet=="Uranus" else 24764) or not 0 <= e < 1:
            raise ValueError("invalid or planet-intersecting trial")
        angles = np.deg2rad([inc,node,peri,mean])
        if model == "circular":
            if any([e,inc,node,peri]):
                raise ValueError("circular grid has e=i=node=peri=0")
            phase = angles[-1] + frequency(planet,a)*times
            pos = a*(np.cos(phase)[...,None]*frame[0] + np.sin(phase)[...,None]*frame[1])
        else:
            pos = integrate(planet, a,e,*angles,times,frame)
        distance = delta*AU+(pos*los).sum(axis=-1)
        x[index] = center[0] - (pos*east).sum(axis=-1)/distance*ARCSEC/scale
        y[index] = center[1] + (pos*north).sum(axis=-1)/distance*ARCSEC/scale
    mid = jd.shape[1]//2
    np.savez_compressed(output, x=x[:,:,mid], y=y[:,:,mid], x_path=x, y_path=y,
                        orbit_id=np.arange(len(trials)).astype(str), exposure_id=ids,
                        elements=trials, a_km=trials[:,0], phase_deg=trials[:,-1],
                        epoch_jd_tdb=epoch, planet=planet, model=model)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cube")
    p.add_argument("output")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--grid", type=float, nargs=4, metavar=("RMIN","RMAX","DR_KM","DPHASE_DEG"))
    g.add_argument("--elements", help="CSV: a_km,e,i_deg,node_deg,peri_deg,mean_deg")
    p.add_argument("--epoch-jd-tdb", type=float)
    args = p.parse_args()
    generate(args.cube, args.output, args.elements, args.grid, epoch=args.epoch_jd_tdb)
