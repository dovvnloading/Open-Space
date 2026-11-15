"""Rebuild orbital timing, width, source-position and historical checks for TIC 35481236."""
import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy.optimize import brentq, least_squares, minimize_scalar

ROOT = Path(__file__).resolve().parent
PERIOD = 1.8434336
PRIMARY = 1493.19101
SECONDARY = 1492.38524 + PERIOD
SEED = 35481236


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def cached_product(row, cache):
    path = cache / row["cache_name"]
    if not path.exists():
        with urllib.request.urlopen(row["url"], timeout=150) as stream:
            payload = stream.read()
        if hashlib.sha256(payload).hexdigest() != row["sha256"]:
            raise ValueError(f"Archive bytes changed for {row['url']}")
        path.write_bytes(payload)
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise ValueError(f"Checksum mismatch for {path}")
    return path


def shape(x, center, width, exponent, depth, level, slope):
    return level + slope*x - depth*np.exp(-0.5*abs((x-center)/width)**exponent)


def bisector(x, y, params):
    center, _, _, depth, level, slope = params
    z = (level+slope*x-y)/depth
    midpoints = []
    for fraction in (0.30, 0.45, 0.60, 0.75):
        crossings = []
        for side in (-1, 1):
            use = (side*(x-center)>0) & (abs(z-fraction)<0.24)
            if use.sum() < 3:
                break
            gradient, intercept = np.polyfit(x[use], z[use], 1)
            if side*gradient >= 0:
                break
            crossings.append((fraction-intercept)/gradient)
        if len(crossings) == 2:
            midpoints.append(np.mean(crossings))
    return np.median(midpoints) if len(midpoints)>=3 else np.nan


def time_eclipses(t, flux, source, sector):
    events = []
    for kind, reference in (("primary", PRIMARY), ("secondary", SECONDARY)):
        first = int(np.floor((t.min()-reference)/PERIOD))
        last = int(np.ceil((t.max()-reference)/PERIOD))
        for cycle in range(first, last+1):
            epoch = reference + cycle*PERIOD
            use = abs(t-epoch)<0.17
            x, y = t[use]-epoch, flux[use]
            if len(y)<10 or min(sum(x<-.08), sum(x>.08))<2:
                continue
            cadence = np.median(np.diff(x))
            level = np.percentile(y,84)
            depth = level-np.percentile(y,4)
            if np.max(np.diff(x))>max(1.8*cadence,.012) or depth<.10:
                continue
            params = least_squares(
                lambda p: shape(x,*p)-y, [0,.04,2,depth,level,0],
                bounds=([-.07,.02,1,.12,.7,-.7],[.06,.09,5,.8,1.3,.7]),
                loss="soft_l1", f_scale=.0025
            ).x
            offset, width, exponent, amplitude, _, _ = params
            events.append(dict(
                product=source, sector=sector, kind=kind, cycle=cycle,
                time_btjd=epoch+offset, bisector_btjd=epoch+bisector(x,y,params),
                fwhm_min=2880*width*(2*np.log(2))**(1/exponent),
                depth=amplitude, cadences=len(y)
            ))
    return events


def pair_events(events):
    primary = {r["cycle"]: r for r in events if r["kind"]=="primary"}
    results = []
    for r in events:
        n = r["cycle"]
        if r["kind"]!="secondary" or n not in primary or n+1 not in primary:
            continue
        phase = {}
        for name, key in (("phase","time_btjd"), ("bisector_phase","bisector_btjd")):
            t0, t1 = primary[n][key], primary[n+1][key]
            phase[name] = (r[key]-t0)/(t1-t0)
        results.append(dict(product=r["product"], sector=r["sector"], cycle=n,
                            secondary_btjd=r["time_btjd"], **phase))
    return results


def localize(path, sector):
    with fits.open(path) as hdus:
        tab = hdus[1].data
        mask = (tab["QUALITY"]==0) & np.isfinite(tab["TIME"])
        t, flux = tab["TIME"][mask], tab["FLUX"][mask].astype(float)
        wcs = WCS(hdus[2].header, naxis=2)
        tx, ty = wcs.world_to_pixel_values(hdus[0].header["RA_OBJ"],hdus[0].header["DEC_OBJ"])
    y, x = np.indices(flux.shape[1:])
    aperture = np.hypot(x-tx,y-ty)<3
    out = []
    for kind, reference in (("primary",PRIMARY),("secondary",SECONDARY)):
        if kind == "primary":
            shift = .004 if sector==108 else .0036
        else:
            shift = -.0195 if sector==108 else -.0156
        images = []
        for n in range(int(np.floor((t.min()-reference)/PERIOD)),int(np.ceil((t.max()-reference)/PERIOD))+1):
            dt = t-(reference+n*PERIOD+shift)
            on = abs(dt)<.022
            left, right = (dt>-.17)&(dt<-.12), (dt>.12)&(dt<.17)
            if min(on.sum(),left.sum(),right.sum())<10:
                continue
            images.append((np.nanmedian(flux[left],axis=0)+np.nanmedian(flux[right],axis=0))/2-
                          np.nanmedian(flux[on],axis=0))
        difference = np.median(images,axis=0)
        def centroid(grid):
            return np.sum(grid[aperture]*difference[aperture])/np.sum(difference[aperture])
        cx, cy = centroid(x), centroid(y)
        out.append(dict(sector=sector,kind=kind,eclipses=len(images),target_x=float(tx),
                        target_y=float(ty),difference_x=cx,difference_y=cy,
                        separation_pix=np.hypot(cx-tx,cy-ty)))
    return out


def observed_phase(ecc, omega_deg, inclination=84.2):
    # True anomaly f=u-omega. Minimize sky-projected squared separation in
    # each conjunction's neighborhood, then convert f to elapsed mean anomaly.
    omega = np.deg2rad(omega_deg)
    sin_i = np.sin(np.deg2rad(inclination))
    def distance2(u):
        return ((1-ecc*ecc)/(1+ecc*np.cos(u-omega)))**2*(1-(sin_i*np.sin(u))**2)
    def mean_anomaly(u):
        f = u-omega
        E = 2*np.arctan2(np.sqrt(1-ecc)*np.sin(f/2),np.sqrt(1+ecc)*np.cos(f/2))
        return E-ecc*np.sin(E)
    u1,u2 = [minimize_scalar(distance2,bounds=(c-.4,c+.4),method="bounded").x
             for c in (np.pi/2,3*np.pi/2)]
    return ((mean_anomaly(u2)-mean_anomaly(u1))%(2*np.pi))/(2*np.pi)


def historical_check(path):
    measurements = []
    for line in path.read_text().splitlines():
        if line.startswith("#dataset="):
            field = int(line.split("=")[1].split(";")[0])
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split()
        if field==3 and cols[-2] in ("A","B") and float(cols[1])<20:
            measurements.append((float(cols[0])-7000, float(cols[1])))
    t, mag = np.array(measurements).T
    flux = 10**(-.4*(mag-np.median(mag)))
    # The two types are far apart. Fit each folded eclipse with its own local
    # baseline, width, amplitude and shape. Time scales shift both centers.
    def phase_for(times, values):
        centers = []
        for reference in [PRIMARY,SECONDARY]:
            x = (times-reference+PERIOD/2)%PERIOD-PERIOD/2
            use = abs(x)<.21
            x, y = x[use], values[use]
            level = np.percentile(y,80)
            depth = max(.05,level-np.percentile(y,6))
            z = least_squares(lambda p:shape(x,*p)-y,[0,.043,2,depth,level,0],
                              bounds=([-.10,.013,1,.08,.80,-1],[.12,.11,6,.65,1.25,1]),
                              loss="soft_l1", f_scale=.018).x
            centers.append(z[0])
        return (SECONDARY+centers[1]-PRIMARY-centers[0])/PERIOD
    rng = np.random.default_rng(SEED)
    out = []
    for tag, mask in (("all", np.ones(len(t),dtype=bool)), ("early",t<np.median(t)),("late",t>=np.median(t))):
        times, values = t[mask],flux[mask]
        # Resample whole observing nights to keep same-night point correlations.
        nights, inv = np.unique(np.floor(times+.5),return_inverse=True)
        groups = [np.flatnonzero(inv==i) for i in range(len(nights))]
        draws = []
        for _ in range(200):
            indices = np.concatenate([groups[i] for i in rng.integers(len(groups),size=len(groups))])
            draws.append(phase_for(times[indices],values[indices]))
        low, high = np.percentile(draws,[2.5,97.5])
        out.append(dict(sample=tag,points=len(times),nights=len(nights),median_hjd=2457000+np.median(times),
                        phase=phase_for(times,values),bootstrap_low=low,bootstrap_high=high))
    return out



def endpoint_repeatability(pairs):
    """Retain alternating cycles so no two triplets share a primary eclipse."""
    product = "SPOC PDCSAP_FLUX"
    def sample(sector, parity, key):
        return [row[key] for row in pairs if row["product"]==product and
                row["sector"]==sector and row["cycle"]%2==parity]
    out = []
    for early in (6,7):
        for method in ("phase","bisector_phase"):
            for p1 in (0,1):
                for p2 in (0,1):
                    a,b = sample(early,p1,method),sample(108,p2,method)
                    out.append(dict(early_sector=early, method=method,
                                    early_parity=p1,late_parity=p2,early_count=len(a),late_count=len(b),
                                    contraction_minutes=(np.median(a)-np.median(b))*PERIOD*1440))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("archive-cache"))
    parser.add_argument("--output", type=Path, default=Path("replay"))
    args = parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    args.output.mkdir(parents=True,exist_ok=True)
    events, pairs, pixels, sector_rows = [], [], [], []
    with (ROOT/"sources.csv").open() as stream:
        sources = list(csv.DictReader(stream))
    for row in sources:
        path = cached_product(row,args.cache)
        kind = row["product"]
        sector = int(row["sector"] or 0)
        if kind=="pixels":
            pixels.extend(localize(path,sector))
        elif kind=="ASAS-3":
            write_csv(args.output/"ground.csv",historical_check(path))
        elif kind in {"SPOC","QLP"}:
            tab = fits.getdata(path)
            keys = ["PDCSAP_FLUX","SAP_FLUX"] if kind=="SPOC" else [
                "DET_FLUX" if "DET_FLUX" in tab.names else "KSPSAP_FLUX"]
            for key in keys:
                t, f = tab["TIME"], tab[key].astype(float)
                good = (tab["QUALITY"]==0) & np.isfinite(t+f)
                ev = time_eclipses(t[good], f[good]/np.median(f[good]), kind+" "+key, sector)
                pr = pair_events(ev)
                events.extend(ev)
                pairs.extend(pr)
                phases = [r["phase"] for r in pr]
                bisections = [r["bisector_phase"] for r in pr if np.isfinite(r["bisector_phase"])]
                summary = dict(product=kind+" "+key, sector=sector,events=len(ev),triplets=len(pr),
                               median_btjd=np.median([r["secondary_btjd"] for r in pr]),
                               phase=np.median(phases),phase_sd=np.std(phases,ddof=1),
                               bisector_phase=np.median(bisections) if bisections else None)
                for ek in ("primary","secondary"):
                    summary[ek+"_fwhm_min"] = np.median([r["fwhm_min"] for r in ev if r["kind"]==ek])
                sector_rows.append(summary)
    for fname, rows in (("eclipses.csv",events),("triplets.csv",pairs),
                        ("sectors.csv",sector_rows),("localization.csv",pixels)):
        write_csv(args.output/fname,rows)
    write_csv(args.output/"endpoint-checks.csv",endpoint_repeatability(pairs))
    pdcsap = {r["sector"]:r for r in sector_rows if r["product"]=="SPOC PDCSAP_FLUX"}
    a,b = pdcsap[6],pdcsap[108]
    dt = b["median_btjd"]-a["median_btjd"]
    result = dict(delta_minutes=(b["phase"]-a["phase"])*PERIOD*1440,baseline_days=dt)
    for e in [.12,.13,.14]:
        angles = [brentq(lambda w: observed_phase(e,w)-f,0,90) for f in (a["phase"],b["phase"])]
        result[str(e)] = dict(omega_early_deg=angles[0],omega_late_deg=angles[1],
                              rate_deg_per_year=(angles[1]-angles[0])*365.25/dt)
    gaia = Table.read(ROOT/"gaia-spectroscopy.xml",format="votable",use_names_over_ids=True)[0]
    if int(gaia["source_id"])!=2946311659129414144 or gaia["nss_solution_type"]!="SB2":
        raise ValueError("Unexpected Gaia spectroscopic source or solution")
    ge = float(gaia["eccentricity"])
    angles = [brentq(lambda w: observed_phase(ge,w)-f,0,90) for f in (a["phase"],b["phase"])]
    result["gaia_spectroscopic_geometry"] = dict(
        eccentricity=ge, eccentricity_error=float(gaia["eccentricity_error"]),
        spectroscopic_omega_deg=float(gaia["arg_periastron"]),
        spectroscopic_omega_error_deg=float(gaia["arg_periastron_error"]),
        tess_omega_early_deg=angles[0], tess_omega_late_deg=angles[1],
        rate_deg_per_year=(angles[1]-angles[0])*365.25/dt)
    mass = 3.10*1.98847e30
    sma = (6.67430e-11*mass*(PERIOD*86400/(2*np.pi))**2)**(1/3)
    result["gr_deg_per_year"] = np.rad2deg(6*np.pi*6.67430e-11*mass/(sma*(1-.13**2)*299792458**2))*365.25/PERIOD
    with (args.output/"orbit.json").open("w") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")


if __name__=="__main__":
    main()
