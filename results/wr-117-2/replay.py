"""Reproduce the WR 117-2 timing search and independent image measurements."""
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
RA, DEC = 277.77770833, -9.80038889


def rows(archive, name):
    return Table.read(archive.read(name).decode().splitlines(), format="ascii.ipac")


def visit_table(table, strict=False):
    clean = []
    for r in table:
        cc = str(r["cc_flags"])[:2]
        marked = cc != "00" if strict else any(c.isupper() for c in cc)
        if marked or str(r["moon_masked"]).zfill(2) != "00" or r["qual_frame"] < 5 or r["dist"] > 2:
            continue
        values = [r[k] for k in ("mjd","w1mpro","w2mpro","w1sigmpro","w2sigmpro","ra","dec")]
        if any(np.ma.is_masked(v) for v in values) or max(values[3:5]) > 0.15:
            continue
        clean.append(values)
    if not clean:
        return np.empty((0, 6))
    arr = np.array(sorted(clean))
    runs = np.split(arr, np.where(np.diff(arr[:, 0]) > 20)[0]+1)
    return np.array([list(np.median(run[:, [0,1,2,5,6]], axis=0)) + [len(run)] for run in runs if len(run)>=4])


def sine_design(t, period):
    angle = (np.asarray(t)-57500)*2*np.pi/period
    return np.column_stack([np.ones(len(angle)), np.sin(angle), np.cos(angle)])


def predictions(visits, period, times):
    coeff = np.linalg.lstsq(sine_design(visits[:,0], period), visits[:,1:3], rcond=None)[0]
    return sine_design(times, period) @ coeff


def look_elsewhere(visits):
    """Shuffle whole two-band visits together and refit the full frequency grid."""
    rng = np.random.default_rng(1172)
    y = visits[:,1:3] - np.mean(visits[:,1:3],axis=0)
    y /= np.sqrt(np.sum(y*y,axis=0))
    frequencies = np.linspace(1/2000,1/50,12000)
    angles = 2*np.pi*frequencies[:,None]*(visits[:,0]-57500)
    a,b = np.sin(angles),np.cos(angles)
    a -= a.mean(axis=1,keepdims=True)
    a /= np.linalg.norm(a,axis=1,keepdims=True)
    b -= b.mean(axis=1,keepdims=True)
    b -= a*np.sum(a*b,axis=1,keepdims=True)
    b /= np.linalg.norm(b,axis=1,keepdims=True)
    projection = np.stack([a,b]).reshape(-1,len(y))
    def peaks(samples):
        dots = projection @ samples.reshape(len(y),-1)
        powers = np.sum(dots.reshape(2,len(frequencies),2,-1)**2,axis=(0,2))/2
        return powers.max(axis=0)
    observed = peaks(y[:,:,None])[0]
    exceedances=0
    for _ in range(32):
        sims=np.stack([y[rng.permutation(len(y))] for _ in range(128)],axis=2)
        exceedances += int(np.sum(peaks(sims)>=observed))
    return {"shuffles":4096,"exceedances":exceedances,"seed":1172}


def fit_aliases(visits):
    times = visits[:,0]
    centered = visits[:,1:3] - np.mean(visits[:,1:3], axis=0)
    norm = np.sum(centered**2, axis=0)
    periods = 1/np.linspace(1/2000,1/50,12000)
    def score(p):
        difference = visits[:,1:3] - predictions(visits,p,times)
        return np.mean(1-np.sum(difference**2,axis=0)/norm)
    scores = np.array([score(p) for p in periods])
    local = np.flatnonzero((scores[1:-1]>scores[:-2]) & (scores[1:-1]>scores[2:]))+1
    best = sorted(local, key=lambda i:-scores[i])[:5]
    return [(float(periods[i]),float(scores[i])) for i in best]


def measure_image(content, obs_id):
    with fits.open(io.BytesIO(content), memmap=False) as hdu:
        hdr = hdu["IMAGE"].header
        x,y = WCS(hdr,hdu).all_world2pix([[RA,DEC]],0)[0]
        image = hdu["IMAGE"].data
        if not (0<=round(x)<image.shape[1] and 0<=round(y)<image.shape[0]):
            return None
        # W is a detector-coordinate lookup with no celestial SIP distortion.
        spectral = WCS(hdr,hdu,key="W")
        wave = spectral.wcs_pix2world([[x,y]],0)[0,0]
        yy,xx = np.indices(image.shape)
        radius = np.hypot(xx-x,yy-y)
        sky = np.nanmedian(image[(radius>4) & (radius<7)])
        flag = int(hdu["FLAGS"].data[round(y),round(x)])
        fluxes = []
        for threshold in (2,2.5,3):
            aperture = radius<threshold
            val = np.sum((image[aperture]-sky).astype('float64'))
            fluxes.append(val*1e6*(6.15/206264.806)**2)
        return [obs_id, float(hdr["MJD-AVG"]),wave,flag,*fluxes]


def matched_fluxes(measurements):
    visits = [[r for r in measurements if low<r[1]<high and 4.5<=r[2]<=4.95
               and r[3] & ~(1<<21)==0 and np.all(np.isfinite(r[4:]))]
              for low,high in [(60900,61000),(61100,61200),(61250,61350)]]
    matches = []
    for c in sorted(visits[2],key=lambda r:r[2]):
        a,b = [min(v,key=lambda r:abs(c[2]-r[2])) for v in visits[:2]]
        if max(abs(c[2]-a[2]),abs(c[2]-b[2]))<=0.025:
            matches.append([a,b,c])
    # A visit can contribute multiple exposures at nearby wavelengths.
    # Give each 0.05 micron bin equal weight.
    bins = {}
    for three in matches:
        key = int(np.floor((three[2][2]-4.5)/0.05))
        bins.setdefault(key,[]).append(np.array([r[4:] for r in three]))
    flux = np.stack([np.median(bins[k],axis=0) for k in sorted(bins)])
    ratios = flux/flux[:,0:1,:]
    return matches, flux, np.median(ratios,axis=0)


def tablewrite(path, header, data):
    with path.open("w",newline="") as out:
        writer = csv.writer(out,lineterminator="\n")
        writer.writerow(header)
        writer.writerows(data)


def main():
    with zipfile.ZipFile(ROOT/"catalog-tables.zip") as archive:
        table = rows(archive,"neowise.tbl")
        v = visit_table(table)
        strict = visit_table(table,True)
        ctrls = {}
        for i in range(1,9):
            c = visit_table(rows(archive,f"control-{i}.tbl"))
            if len(c):
                ctrls[str(i)] = [len(c), *np.ptp(c[:,1:3],axis=0)]
        early = []
        for mission in ("4band","3band"):
            tab = rows(archive,f"wr117-allsky_{mission}_p1bs_psd.tbl")
            mask=(np.asarray(tab["dist"])<=2) & np.array([not any(ch.isupper() for ch in str(c)[:2]) for c in tab["cc_flags"]])
            early.append([float(np.median(tab[k][mask])) for k in ("mjd","w1mpro","w2mpro")])
    aliases = fit_aliases(v)
    shuffle = look_elsewhere(v)
    early=np.array(early)
    summary = {"neowise_n":len(table),"retained_n":int(v[:,5].sum()),"visits":len(v),
               "ranges":np.ptp(v[:,1:3],axis=0).tolist(),"strict_visits":len(strict),
               "strict_ranges":np.ptp(strict[:,1:3],axis=0).tolist(),"controls":ctrls,
               "early_2010":early.tolist(),"alias_fits":[],"whole_visit_shuffle":shuffle}
    tablewrite(ROOT/"visits.csv",["mjd","w1_mag","w2_mag","ra_deg","dec_deg","exposures"],v)
    for p,s in aliases:
        summary["alias_fits"].append({"period_d":p,"mean_band_r2":s,
            "early_predicted":predictions(v,p,early[:,0]).tolist()})
    source_rows = list(csv.DictReader((ROOT/"spherex-manifest.csv").open()))
    measured = []
    with zipfile.ZipFile(ROOT/"spherex-cutouts.zip") as archive:
        for r in source_rows:
            if not r["cutout_sha256"]:
                continue
            data=archive.read(r["obs_id"]+".fits")
            if hashlib.sha256(data).hexdigest()!=r["cutout_sha256"]:
                raise ValueError(f'Checksum mismatch for {r["obs_id"]}')
            row=measure_image(data,r["obs_id"])
            if row is not None: measured.append(row)
    tablewrite(ROOT/"spherex-photometry.csv",["obs_id","mjd","wavelength_um","central_flags","jy_r2","jy_r2_5","jy_r3"],measured)
    matches,flux,ratios=matched_fluxes(measured)
    tablewrite(ROOT/"matched-images.csv",["first_id","second_id","third_id","first_um","second_um","third_um"],
               [[r[0] for r in m]+[r[2] for r in m] for m in matches])
    summary["spherex"]={"queried":len(source_rows),"included_cutouts":sum(bool(r["cutout_sha256"]) for r in source_rows),
        "on_cutout":len(measured),"matched_triplets":len(matches),"wavelength_bins":len(flux),
        "median_flux_jy_r2_5":np.median(flux[:,:,1],axis=0).tolist(),
        "median_ratios_r2_r2_5_r3":ratios.tolist(),"per_bin_ratio_r2_5":(flux[:,:,1]/flux[:,0:1,1]).tolist()}
    (ROOT/"result.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":
    import warnings
    from astropy import log
    log.setLevel("ERROR")
    warnings.filterwarnings("ignore",message=".*invalid value.*")
    main()
