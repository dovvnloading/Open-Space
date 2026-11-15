"""Reproduce the RGZ × DESI screen, position controls, and two field measurements."""
import csv
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs"


def separation(ra, dec, ra2, dec2):
    """Great-circle separation in arcseconds, using decimal degrees."""
    hav = (np.sin(np.radians(dec-dec2)/2)**2 +
           np.cos(np.radians(dec))*np.cos(np.radians(dec2))*np.sin(np.radians(ra-ra2)/2)**2)
    return 7200*np.degrees(np.arcsin(np.sqrt(hav)))


def scale(z):
    distance = (299792.458/70) * quad(lambda x: (0.3*(1+x)**3+0.7)**-0.5, 0, z)[0]
    return distance/(1+z)*math.pi/648000


def catalog(path):
    lines = [line for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join([lines[0], *lines[3:]])), sep="\t", skipinitialspace=True)


def crossmatch(hosts, spectra, dec_offset=0):
    matches = []
    # DESI sometimes observes the same target across more than one program.
    spectra = spectra.query("z > 0 and zwarn == 0 and deltachi2 > 25").sort_values("deltachi2", ascending=False)
    for source in spectra.drop_duplicates("targetid").itertuples():
        dist = separation(hosts.ra, hosts.dec+dec_offset, source.mean_fiber_ra, source.mean_fiber_dec)
        h = hosts.loc[dist.idxmin()]
        if dist.min() < 3:
            matches.append(dict(CatID=int(h.CatID), rgz=h.rgz, targetid=int(source.targetid),
                                separation_arcsec=float(dist.min()), z=source.z,
                                size_mpc=scale(source.z)*h.LAE))
    return pd.DataFrame(matches)


def image_cutout(path):
    """Read the single 280x280 float image in the frozen Legacy Surveys TAN cutout."""
    raw = path.read_bytes()
    keys = {}
    for i in range(0, len(raw), 80):
        card = raw[i:i+80].decode("ascii")
        if card[:8].strip() == "END":
            start = (i+80+2879)//2880*2880
            break
        if card[8:10] == "= ":
            keys[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
    else:
        raise ValueError("Missing FITS END card")
    assert keys['BITPIX'] == '-32' and keys['NAXIS'] == '2'
    n, m = int(keys['NAXIS1']), int(keys['NAXIS2'])
    return np.frombuffer(raw, '>f4', count=n*m, offset=start).reshape(m,n).astype(float), keys


def field(row, spectrum):
    tag = row.rgz.split("_")[1][:7]
    first = catalog(INPUT / (tag+"first.dat"))
    ra, dec = spectrum.mean_fiber_ra, spectrum.mean_fiber_dec
    dist = separation(ra, dec, first._RAJ2000, first._DEJ2000)
    # The side of the nucleus containing each outer hotspot was selected on FIRST.
    pairs = {"J000314": ["J000310.0+111503","J000319.3+111612"],
             "J085850": ["J085845.0+262258","J085855.2+262338"]}
    outer = first.set_index("FIRST").loc[pairs[tag]]
    west, east = outer.iloc[0], outer.iloc[1]
    baseline = separation(west._RAJ2000, west._DEJ2000, east._RAJ2000, east._DEJ2000)
    core = first.loc[dist.idxmin()]
    points = np.column_stack(((outer._RAJ2000-ra)*np.cos(np.radians(dec)), outer._DEJ2000-dec))
    angle = np.degrees(np.arccos(np.clip(np.dot(*points)/np.prod(np.linalg.norm(points,axis=1)), -1, 1)))
    image, h = image_cutout(INPUT/(tag+"_vlass.fits"))
    xs, ys = np.indices(image.shape)[::-1]
    med = np.median(image)
    noise = 1.4826*np.median(np.abs(image-med))
    peaks = []
    for _, r in outer.iterrows():
        x = (r._RAJ2000-float(h['CRVAL1']))*np.cos(np.radians(float(h['CRVAL2'])))/float(h['CD1_1'])+float(h['CRPIX1'])-1
        y = (r._DEJ2000-float(h['CRVAL2']))/float(h['CD2_2'])+float(h['CRPIX2'])-1
        local = np.where((xs-x)**2+(ys-y)**2 <= 16, image, -np.inf)
        yp,xp = np.unravel_index(np.argmax(local),local.shape)
        peaks.append((int(xp),int(yp),round(float(local[yp,xp]/noise),1)))
    vl_sep = math.dist(peaks[0][:2],peaks[1][:2])
    opt = pd.read_csv(INPUT/(tag+"opt.dat"))
    opt["sep"] = separation(ra,dec,opt.ra,opt.dec)
    return dict(rgz=row.rgz, CatID=int(row.CatID), targetid=int(spectrum.targetid),
                z=spectrum.z, zerr=spectrum.zerr, redrock_delta_chi2=spectrum.deltachi2,
                core=core.FIRST, core_peak_mJy_per_beam=core.Fpeak,
                core_fiber_sep_arcsec=float(dist.min()),
                optical_DEV_fiber_sep_arcsec=float(opt.loc[opt.type=='DEV','sep'].min()),
                neighbor_PSF_fiber_sep_arcsec=float(opt.loc[opt.type=='PSF','sep'].min()),
                RGZ_LAE_arcsec=row.LAE, RGZ_box_Mpc=scale(spectrum.z)*row.LAE,
                FIRST_peak_names=pairs[tag], FIRST_peak_sep_arcsec=baseline,
                FIRST_peak_Mpc=baseline*scale(spectrum.z), opening_angle_deg=angle,
                VLASS_peaks_xy_snr=peaks, VLASS_sep_arcsec=vl_sep, VLASS_peak_Mpc=vl_sep*scale(spectrum.z))


def main():
    hosts = pd.read_csv(INPUT/"selection.csv")
    spectra = pd.read_csv(INPUT/"desi.csv")
    matches = crossmatch(hosts,spectra)
    shifted = crossmatch(hosts,pd.read_csv(INPUT/"shifted_desi.csv"),30/3600)
    near = pd.read_csv(INPUT/"catalog_neighbors.csv")
    nearest = near.groupby("CatID").separation_arcsec.min()
    matches["nearest_giant_arcsec"] = matches.CatID.map(nearest)
    matches = matches.sort_values("size_mpc",ascending=False)
    matches.to_csv(HERE/"screen.csv",index=False)
    giants = matches.query("size_mpc > 0.7")
    print(f"{len(hosts)} hosts, {len(matches)} matches, {len(shifted)} shifted matches, "
          f"{len(giants)} over 0.7 Mpc, {sum(giants.nearest_giant_arcsec<=180)} with catalog hits")
    results=[]
    for cat_id in (8724,142770):
        row = hosts[hosts.CatID==cat_id].iloc[0]
        target = matches[matches.CatID==cat_id].targetid.iloc[0]
        spectrum = spectra[spectra.targetid==target].iloc[0]
        results.append(field(row,spectrum))
    (HERE/"measurements.json").write_text(json.dumps(results,indent=2)+"\n")


if __name__ == "__main__":
    main()
