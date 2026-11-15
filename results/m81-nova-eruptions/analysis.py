"""Recalculate the M81 catalogue search and Gaia registrations. See README.md."""
import argparse
import csv
import hashlib
import itertools
import re
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from lxml import html
from PIL import Image
from scipy import optimize, spatial

HERE = Path(__file__).resolve().parent
CACHE = HERE / "inputs"
RA = 15 * (9 + 55 / 60 + 54.39 / 3600)
DEC = 69 + 50.2 / 3600
CENTER = (15 * (9 + 55 / 60 + 33.173 / 3600), 69 + 3 / 60 + 55.06 / 3600)
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode({
    "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
    "QUERY": f"SELECT source_id,ra,dec,pmra,pmdec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.05))",
})
URLS = {
    "mpe.html": "https://www.mpe.mpg.de/~m31novae/opt/m81/M81_table.html",
    "subaru.jpg": "https://www.rochesterastronomy.org/sn2000/nm8101n1.jpg",
    "ondrejov.png": "https://www.asu.cas.cz/~asteroid/PNV_J09555439+6900502.png",
    "rozhen.png": "https://drive.google.com/uc?export=download&id=1M7QGltC1WvtNHV82sjtbV3KSj6-S2DNj",
    "gaia.csv": GAIA,
    "hubble.csv": "https://catalogs.mast.stsci.edu/api/v0.1/hsc/v3/summary.csv?" +
        urllib.parse.urlencode({"ra": RA, "dec": DEC, "radius": 0.001, "pagesize": 1000}),
}
# Each plate's rough Gaia-matching seed is (arcsec/pixel, rotation radians,
# x, y of the nominal sky position). It only locates reference-star windows;
# the fit uses measured star centres, with six independent affine coefficients.
PLATES = [
    ("Subaru", "subaru.jpg", 2000.98, (0.20159, -0.00054, 383.3, 366.46), (385, 365)),
    ("Ondrejov", "ondrejov.png", 2020.56, (0.52467, 0.019, 309.8, 490.04), (310, 490)),
    ("Rozhen R", "rozhen.png", 2020.56, (0.44633, -0.04071, 287.2, 455.69), (287, 454)),
    ("Rozhen Ha", "rozhen.png", 2020.56, (0.45066, -0.03234, 275.77, 454.56), (276, 454)),
]


def tsv(name):
    with (HERE / name).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(name, header, rows):
    with (HERE / name).open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fetch():
    CACHE.mkdir(exist_ok=True)
    for name, url in URLS.items():
        loc = CACHE / name
        if loc.exists():
            continue
        data = urllib.request.urlopen(url, timeout=90).read()
        if name in ("gaia.csv", "hubble.csv"):
            # TAP does not guarantee a stable row order. Canonicalize before hashing.
            lines = data.decode().splitlines()
            data = ("\n".join([lines[0], *sorted(lines[1:])]) + "\n").encode()
        loc.write_bytes(data)


def check_hashes():
    for row in tsv("source-hashes.tsv"):
        data = (CACHE / row["file"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Input differs from the study snapshot: {row['file']}")


def parse_catalogue():
    rows = []
    for tr in html.parse(str(CACHE / "mpe.html")).xpath("//tr[td]"):
        cells = [" ".join(c.text_content().split()) for c in tr.xpath("td")]
        if "no M81 nova" in cells[0]:
            continue
        fields = re.findall(r"[+-]?\d+(?:\.\d*)?", cells[1])
        if len(fields) != 6:
            print(f"Excluding malformed coordinate for {cells[0]}: {cells[1]}")
            continue
        h, m, s, d, dm, ds = map(float, fields)
        rows.append((cells[0], cells[4], f"{15*(h+m/60+s/3600):.9f}",
                     f"{d+dm/60+ds/3600:.9f}", cells[2]))
    write_tsv("catalogue.tsv", ["event", "jd", "ra_deg", "dec_deg", "error_as_listed_arcsec"], rows)


def prob_stats():
    entries = tsv("catalogue.tsv")
    coords = np.array([[float(r["ra_deg"]), float(r["dec_deg"])] for r in entries])
    # Project all points on the tangent plane about the centre of M81.
    ra, dec = np.radians(coords.T)
    ra0, dec0 = np.radians(CENTER)
    proj = np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0) + np.sin(dec0)*np.sin(dec)
    xy = np.rad2deg(np.column_stack((np.cos(dec)*np.sin(ra-ra0),
        np.cos(dec0)*np.sin(dec)-np.sin(dec0)*np.cos(dec)*np.cos(ra-ra0)))/proj[:,None])*3600
    jd = np.array([float(r["jd"]) for r in entries])
    eligible = [(i, j) for i, j in itertools.combinations(range(len(entries)), 2) if abs(jd[i]-jd[j]) > 90]
    distance = lambda i, j: float(np.linalg.norm(xy[i]-xy[j]))
    ranked = sorted(eligible, key=lambda pair: distance(*pair))
    i, j = ranked[0]
    d = distance(i, j)
    print(f"{len(entries)} entries, {len(eligible)} pairs >90 days apart")
    for a, b in ranked[:6]:
        print(entries[a]["event"], entries[b]["event"], round(distance(a,b),3))
    print(f"Discovery-date interval {abs(jd[i]-jd[j]):.2f} d")
    # ATel 13440 independently gives a sexagesimal position for the one
    # malformed discovery-catalogue row. Include it in a sensitivity run.
    repair_ra, repair_dec = 15*(9+55/60+48.50/3600), 69+7/60+33.2/3600
    rra, rdec = np.radians([repair_ra, repair_dec])
    scale = np.cos(dec0)*np.cos(rdec)*np.cos(rra-ra0) + np.sin(dec0)*np.sin(rdec)
    repair_pos = np.rad2deg(np.array([np.cos(rdec)*np.sin(rra-ra0),
        np.cos(dec0)*np.sin(rdec)-np.sin(dec0)*np.cos(rdec)*np.cos(rra-ra0)])/scale)*3600
    extra = np.linalg.norm(xy[np.abs(jd-2458879.51)>90]-repair_pos,axis=1)
    print(f"With ATel 13440 repair: {len(entries)+1} entries, {len(eligible)+len(extra)} pairs; "
          f"closest added pair {min(extra):.3f} arcsec; original pair rank {1+sum(extra<d)}")

    old, new = (i, j) if jd[i] < jd[j] else (j, i)
    previous = xy[jd < jd[new]-90]
    distances = np.linalg.norm(previous-xy[new],axis=1)
    for r in [30, 60, 120]:
        others = int(np.sum((distances <= r) & (distances > d+1e-6)))
        print(f"Local: {others} other earlier entries within {r} arcsec; expected {others*d*d/(r*r-d*d):.6f} within catalogue offset")
    # Use all possible pairs for a look-elsewhere estimate. For each entry in
    # turn, displace it locally; count only eligible, earlier counterparts.
    # Fixed radii, uniform alternative directions, preserving actual sky density.
    ang = np.arange(720)*2*np.pi/720
    tree = spatial.cKDTree(xy)
    for radius in [5, 10, 20]:
        vectors = radius * np.column_stack((np.cos(ang),np.sin(ang)))
        counts = []
        for v in vectors:
            neighbors = tree.query_ball_point(xy+v, d)
            counts.append(sum(sum(jd[p] < jd[k]-90 for p in js) for k,js in enumerate(neighbors)))
        print(f"Shift {radius} arcsec: mean pairs {np.mean(counts):.3f}, fraction >=1 {np.mean(np.array(counts)>0):.3f}")

    for k in [5,10,20]:
        # Each later event is assigned a local scale from earlier eruptions.
        # Omit its nearest neighbor from the k background tracers to avoid
        # shrinking that scale because of the candidate itself.
        earlier = [xy[jd<t-90] for t in jd]
        kj = []
        for z,pos in zip(earlier,xy):
            if len(z) < k+1:
                kj.append(np.nan)
            else:
                kj.append(np.sort(np.linalg.norm(z-pos,axis=1))[k])
        kj = np.asarray(kj)
        scores = [(distance(a,b)/kj[a if jd[a]>jd[b] else b],a,b) for a,b in eligible
                  if np.isfinite(kj[a if jd[a]>jd[b] else b])]
        scores.sort()
        ratio = d/kj[new]
        print(f"Density-scaled screen (k={k}): rank {1+sum(e[0]<ratio-1e-8 for e in scores)}, score {ratio:.5f}")
        for radius in [5,10,20]:
            # Recompute nearest-neighbor distances for every eligible older
            # sample under every displacement; same search and statistic.
            shifts = radius*np.column_stack((np.cos(ang),np.sin(ang)))
            mins = np.full(len(shifts),np.inf)
            for pos,z in zip(xy,earlier):
                if len(z)>=k+1:
                    dist=spatial.cKDTree(z).query(pos+shifts,k=[1,k+1])[0]
                    mins=np.minimum(mins,dist[:,0]/dist[:,1])
            print(f"  {radius} arcsec shift: fraction at least this extreme {np.mean(mins<=ratio):.3f}")


def star_fit(im, x, y, radius=6):
    x, y = round(x), round(y)
    yy, xx = np.mgrid[-radius:radius+1, -radius:radius+1]
    patch = im[y-radius:y+radius+1, x-radius:x+radius+1]
    mask = (xx**2+yy**2 <= radius**2)
    z, xx, yy = patch[mask], xx[mask], yy[mask]
    low, high = np.percentile(patch,[20,99])
    def residual(p):
        cx, cy, amp, sig, bg, dx, dy = p
        return bg+dx*xx+dy*yy+amp*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*sig**2))-z
    fit = optimize.least_squares(residual, [0,0,max(high-low,5),2,low,0,0],
        bounds=([-4,-4,2,0.8,-200,-12,-12],[4,4,1000,5,300,12,12]),
        loss="soft_l1", f_scale=3).x
    if max(abs(fit[:2])) > 3.9 or fit[3] >= 4.95:
        return None
    return fit[:2]+[x,y], fit[2], fit[3]


def aperture(im, target, radius):
    x,y = target
    yy,xx = np.mgrid[-radius:radius+1,-radius:radius+1]
    z = im[y-radius:y+radius+1,x-radius:x+radius+1]
    mask = xx**2+yy**2 <= radius**2
    weights = np.maximum(z-np.median(z[~mask]),0)*mask
    return np.array([x+np.sum(xx*weights)/weights.sum(),y+np.sum(yy*weights)/weights.sum()])


def fit_affine(stars):
    return np.linalg.lstsq(np.c_[stars[:,:2],np.ones(len(stars))], stars[:,2:4],rcond=None)[0]


def controls(stars,target):
    m = fit_affine(stars)
    resid = np.c_[stars[:,:2],np.ones(len(stars))]@m-stars[:,2:4]
    sky = np.r_[target,1]@m
    held = []
    for i,row in enumerate(stars):
        other = np.delete(stars,i,axis=0)
        held.append(np.linalg.norm(np.r_[row[:2],1]@fit_affine(other)-row[2:4]))
    return sky, [round(float(v),4) for v in np.std(resid,axis=0)], round(float(np.median(held)),4)


def plate_measurements(name, fname, year, seed, target, gaia):
    im = np.asarray(Image.open(CACHE / fname).convert("L"),dtype=float)
    if name == "Rozhen R":
        im = im[:,:555]
    if name == "Rozhen Ha":
        im = im[:,565:]
    s,t,x0,y0 = seed
    refs = []
    for row in gaia:
        west = -(float(row["ra"])-RA)*np.cos(np.radians(DEC))*3600
        south = -(float(row["dec"])-DEC)*3600
        west -= float(row["pmra"] or 0)*(year-2016)/1000
        south -= float(row["pmdec"] or 0)*(year-2016)/1000
        x = x0+(west*np.cos(t)-south*np.sin(t))/s
        y = y0+(west*np.sin(t)+south*np.cos(t))/s
        if min(x,y)<18 or x>im.shape[1]-18 or y>im.shape[0]-18:
            continue
        if (x0-15<x<x0+130 and abs(y-y0)<6) or (abs(x-x0)<5 and y0-5<y<y0+85):
            continue
        if name=="Subaru" and min(abs(x-632),abs(x-730))<15:
            continue
        if name.startswith("Rozhen") and (y<105 or y>520 or x>im.shape[1]-50):
            continue
        f = star_fit(im,x,y)
        if f is not None and f[1]>17 and 1.15<f[2]<4:
            refs.append([*f[0],west,south,int(row["source_id"])])
    stars = np.array([r[:4] for r in refs])
    ids = np.array([str(r[4]) for r in refs])
    for limit in [1.2,0.75,0.55]:
        residual = np.c_[stars[:,:2],np.ones(len(stars))]@fit_affine(stars)-stars[:,2:4]
        keep = np.linalg.norm(residual,axis=1)<limit
        stars,ids = stars[keep],ids[keep]
    tgt = aperture(im,target,6) if name.startswith("Rozhen") else star_fit(im,*target)[0]
    return stars,ids,tgt


def full():
    fetch()
    check_hashes()
    parse_catalogue()
    with (CACHE/"gaia.csv").open() as f:
        gaia = sorted(csv.DictReader(f),key=lambda row:float(row["phot_g_mean_mag"]))
    with (CACHE/"hubble.csv").open() as f:
        hubble = list(csv.DictReader(f))
    write_tsv("hubble.tsv", ["match_id","ra_deg","dec_deg","f814w","n_visits"],
              [(r["MatchID"],r["MatchRA"],r["MatchDec"],r["A_F814W"],r["NumVisits"]) for r in hubble])
    rows,objects = [],[]
    for name,fn,year,seed,t in PLATES:
        stars,ids,tgt = plate_measurements(name,fn,year,seed,t,gaia)
        rows.extend((name, sid, *(f"{v:.8f}" for v in row)) for row,sid in zip(stars,ids))
        objects.append((name,*(f"{v:.8f}" for v in tgt)))
    write_tsv("references.tsv",["plate","gaia_dr3_source_id","x","y","west_as","south_as"],rows)
    write_tsv("targets.tsv",["plate","x","y"],objects)


def measure():
    refs = tsv("references.tsv")
    skies = []
    for target in tsv("targets.tsv"):
        name = target["plate"]
        stars = np.array([[float(r[k]) for k in ("x","y","west_as","south_as")]
                          for r in refs if r["plate"]==name])
        sky,rms,held = controls(stars,[float(target["x"]),float(target["y"])])
        skies.append(sky)
        print(f"{name}: {len(stars)} reference stars, west/south {sky.round(3)}, "
              f"scatter {rms}, median leave-one-out radial residual {held:.3f} arcsec")
    hsc = tsv("hubble.tsv")
    hxy = np.array([[-(float(r["ra_deg"])-RA)*np.cos(np.radians(DEC))*3600,
                     -(float(r["dec_deg"])-DEC)*3600] for r in hsc])
    for p, xy in zip(PLATES, skies):
        print(f"{p[0]}: {np.sum(np.linalg.norm(hxy-xy,axis=1)<0.5)} HSC MatchIDs within 0.5 arcsec")
    for p,xy in zip(PLATES[1:],skies[1:]):
        print(f"Subaru to {p[0]}: {np.linalg.norm(xy-skies[0]):.3f} arcsec")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-images",action="store_true",help="download fixed inputs and repeat source measurements")
    if parser.parse_args().from_images:
        full()
    prob_stats()
    measure()
