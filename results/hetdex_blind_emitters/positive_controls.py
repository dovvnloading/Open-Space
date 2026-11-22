"""Recover distinct published He II systems; QA never determines the thresholds."""
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

import spectral as spec
from acquisition import plain
from desi_lines import desi_lines
from vetting import Cache, DESI, DESI_BACKUP, PRIOR, separation


def desi_positive():
    """First catalog-order, secure, moderate-z [O II] galaxy as an external fit check."""
    cache = Cache(Path(os.environ.get("RUNNER_TEMP","/tmp"))/"hetdex")
    name = DESI.rsplit("/",1)[1]
    try:
        path = cache.fetch("catalog",name,[DESI,DESI_BACKUP])
        with fits.open(path,memmap=True,character_as_bytes=True) as hdu:
            table = hdu["INFO"].data
            for idx in range(len(table)):
                row = table[idx]
                z = float(row["vi_z"])
                typ = str(plain(row["source_type"])).strip().lower()
                if not (0.05 <= z <= 0.45 and float(row["vi_quality"]) >= 3 and typ=="oii"):
                    continue
                center = spec.vac_to_air(3728.4835*(1+z))
                measurements = desi_lines(hdu,idx,center,"unknown")
                return {"status":"tested","vac_row":idx,"vi_z":z,"vi_quality":float(row["vi_quality"]),
                        "targetid":int(row["targetid"]),"source_type":typ,"air_wave_predicted":center,
                        "line_measurements":measurements["lines"],"input":cache.inputs[name]}
        return {"status":"no_secure_oii_row","input":cache.inputs[name]}
    except (OSError,ValueError,KeyError,IndexError) as error:
        return {"status":"unavailable","error":str(error),"input":cache.inputs.get(name)}


def positive_controls(cat,out):
    """Search complete spectra without prior information; match after detection."""
    from scan import channels, json_default
    ra,dec = np.asarray(cat.info["ra"],float), np.asarray(cat.info["dec"],float)
    identifiers = np.asarray(cat.info["detectid"])
    results = []
    dx = float(np.median(np.diff(cat.wave)))
    for ident,r,d,z in PRIOR:
        box = (np.abs((ra-r+180)%360-180)*np.cos(np.deg2rad(d)) < 3.01/3600) & (np.abs(dec-d) < 3.01/3600)
        nearby = np.flatnonzero(box)
        offsets = separation(ra[nearby],dec[nearby],r,d)
        matched = np.union1d(nearby[offsets<=3.0],np.flatnonzero(identifiers==ident))
        expected = spec.vac_to_air(spec.UV["heii"]*(1+z))
        entries = []
        for i in matched:
            meta = cat.metadata(int(i))
            y,e = cat.read(i,i+1)
            peaks = spec.search(cat.wave,y[0],e[0],density_unit=1e-17)
            # ±500 km/s plus rounding of a 3-decimal literature redshift.
            near = [p for p in peaks if abs(p["wave"]-expected) <= expected*500/spec.C+0.85]
            peak = min(near,key=lambda p:abs(p["wave"]-expected),default=None)
            zcat = meta.get("z_hetdex")
            zcat = float(zcat) if zcat is not None and np.isfinite(zcat) and zcat>=0 else None
            selected = channels(cat.wave,y[0],e[0],peak,zcat) if peak is not None else []
            f, se, ok = spec._inputs(cat.wave, y, e, density_unit=1e-17)
            _, _, scales = spec._screen(f,se,ok,dx)
            f,se = f[0],se[0]*scales[0]
            # Prespecified width/center independent of which blind peaks survived.
            fixed = {k:spec._force(cat.wave,f,se,dx,rest*(1+z),2.4) for k,rest in spec.UV.items()}
            # Offset/width exploration is diagnostic. Report the number of trials;
            # grid S/N is not an additional fixed-center detection significance.
            trials = []
            for velocity in (-450,-300,-150,0,150,300,450):
                for width in (2.4,4.,6.5):
                    m = spec._force(cat.wave,f,se,dx,spec.UV["heii"]*(1+z)*(1+velocity/spec.C),width)
                    if m["status"]=="measured":
                        trials.append(dict(m,offset_km_s=velocity,sigma_A=width))
            reasons=[]
            if peak is None:
                reasons.append("no nearby blind peak at fitted S/N >= 7")
            else:
                if peak["sky_near"]: reasons.append("sky proximity")
                if peak["single_pixel"]: reasons.append("single-pixel flag")
                if peak["z_heii"] is None: reasons.append("outside He II redshift interval")
                anchored = (zcat is not None and peak["z_heii"] is not None and
                            abs(zcat-peak["z_heii"]) < .008)
                if not anchored and (peak["lya_pre_snr"] or 0) < 4.0:
                    reasons.append("no redshift anchor and Ly-alpha precursor below 4")
            entries.append({"row":int(i),"detectid":int(identifiers[i]),"source_id":int(cat.src_ids[i]),
                            "shotid":int(cat.shot_ids[i]),"separation_arcsec":float(separation(ra[i],dec[i],r,d)),
                            "exact_detectid":bool(identifiers[i]==ident),"source_type":meta.get("source_type"),
                            "catalog_z":zcat,"blind_peak":peak,"blind_channels":selected,
                            "heii_screen_passed":"heii" in selected,"screen_reasons":reasons,
                            "at_literature_z_sigma2p4":fixed,"heii_grid_valid_trials":len(trials),
                            "heii_grid_best":max(trials,key=lambda x:x["snr"],default=None)})
        results.append({"published_detectid":ident,"ra":r,"dec":d,"published_z":z,
                        "heii_expected_air":expected,"matching_catalog_rows":entries,
                        "matched":bool(entries),"detected_near_heii":any(e["blind_peak"] is not None for e in entries),
                        "recovered_heii_screen":any(e["heii_screen_passed"] for e in entries)})
    report = {"release":cat.release,"publication":"10.3847/1538-4357/aea6d9, Table 1",
              "targets":len(results),"position_or_id_matches":sum(t["matched"] for t in results),
              "target_blind_heii_peaks":sum(t["detected_near_heii"] for t in results),
              "target_heii_screen_recoveries":sum(t["recovered_heii_screen"] for t in results),
              "controls":results}
    # The external VAC is identical for both releases. Run the independent check once.
    if cat.release=="sc2":
        report["desi_known_oii"]=desi_positive()
    dest = Path(out)/"positive-controls.json"
    dest.write_text(json.dumps(report,indent=2,default=json_default,allow_nan=False)+"\n")
    return report
