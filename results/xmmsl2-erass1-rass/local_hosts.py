"""Test bright galaxy coincidences against the surrounding Legacy Surveys field."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
import json
import math
from pathlib import Path

from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
import requests

DL = "https://datalab.noirlab.edu/query/query"
R_FIELD, R_INNER = 210., 45.
CONTROL_RADII = (55.,85.,115.,145.,175.)
CONTROL_ANGLES = np.arange(0,360,15.)
GALAXY = {"EXP","DEV","SER","REX","COMP"}


def rate_candidates(center, pos, brightness, is_gal, field_brightness, gate, sigma):
    # Pcc follows the local Poisson probability of at least one galaxy this
    # bright within r. Include astrometric uncertainty in the area, even for
    # exact centroid coincidences (Bloom et al. 2002, AJ 123, 1111).
    sep=center.separation(pos).arcsec
    possible=np.flatnonzero(is_gal & (sep<=gate))
    if not len(possible):
        return []
    area=np.pi*(R_FIELD**2-R_INNER**2)
    out=[]
    for i in possible:
        radius=max(float(sep[i]),1.515*sigma)
        n=int((field_brightness>=brightness[i]).sum())
        # Add one object when estimating the field density. This avoids zero
        # probabilities for a source brighter than every comparison galaxy.
        density=(n+1)/area
        p=-math.expm1(-np.pi*radius**2*density)
        out.append((p,int(i),float(sep[i]),n))
    return sorted(out)


def measure(source):
    sid=source["detection_uid"]
    x,y=source["ra"],source["dec"]
    coord=SkyCoord(x*u.deg,y*u.deg)
    gate=min(35,max(5,2.45*source["poserr_arcsec"]))
    sql=("SELECT ls_id,type,ra,dec,flux_z,flux_r FROM ls_dr10.tractor "
         f"WHERE q3c_radial_query(ra,dec,{x:.9f},{y:.9f},{R_FIELD/3600:.9f}) AND flux_z>0")
    rec={"detection_uid":sid,"ra":x,"dec":y,"gate_arcsec":gate,"query":sql}
    rows=[]
    try:
        r=requests.get(DL,params={"sql":sql,"ofmt":"csv","async":"False"},
                       headers={"X-DL-AuthToken":"anonymous.0.0.anon_access"},timeout=(12,65))
        r.raise_for_status()
        rec.update(sha256=hashlib.sha256(r.content).hexdigest(),bytes=len(r.content),url=r.url)
        reader=csv.DictReader(io.StringIO(r.text))
        if not {"ls_id","ra","dec","type","flux_z"}<=set(reader.fieldnames or ()):
            raise ValueError("invalid Data Lab field response")
        frame=pd.DataFrame(list(reader))
        if not len(frame):
            rec["status"]="empty field or outside imaging footprint"
            return rec,rows
        frame=frame.drop_duplicates("ls_id").reset_index(drop=True)
        if len(frame)>=100000:
            raise ValueError("field export may be truncated")
        brightness=pd.to_numeric(frame.flux_z,errors="coerce").to_numpy()
        pos=SkyCoord(frame.ra.to_numpy(dtype=float)*u.deg,frame.dec.to_numpy(dtype=float)*u.deg)
        is_gal=frame.type.str.strip().isin(GALAXY).to_numpy() & np.isfinite(brightness)
        distance=coord.separation(pos).arcsec
        annulus=(distance>R_INNER)&(distance<=R_FIELD)
        ref=brightness[annulus & is_gal]
        # Empty sectors near imaging edges cannot calibrate a sparse field.
        angles=coord.position_angle(pos).wrap_at(360*u.deg).degree%360
        counts=np.histogram(angles[annulus],bins=np.arange(0,361,45))[0]
        rec.update(field_rows=len(frame),annulus_galaxies=len(ref),sector_counts=counts.tolist())
        if np.any(counts<3):
            rec["status"]="insufficient field or imaging coverage"
            return rec,rows
        sigma=source["poserr_arcsec"]
        candidates=rate_candidates(coord,pos,brightness,is_gal,ref,gate,sigma)
        offsets=coord.directional_offset_by(np.tile(CONTROL_ANGLES,len(CONTROL_RADII))*u.deg,
                                      np.repeat(CONTROL_RADII,len(CONTROL_ANGLES))*u.arcsec)
        blank=[]
        for c in offsets:
            alts=rate_candidates(c,pos,brightness,is_gal,ref,gate,sigma)
            blank.append(alts[0][0] if alts else 1.)
        rec.update(status="ok",control_positions=len(blank),control_scores=blank)
        for rank,(p,i,d,n) in enumerate(candidates,start=1):
            k=int((np.array(blank)<=p).sum())
            rows.append({"detection_uid":sid,"rank":rank,"ls_id":frame.iloc[i].ls_id,
                         "type":frame.iloc[i].type,"host_ra":frame.iloc[i].ra,"host_dec":frame.iloc[i].dec,
                         "sep_arcsec":d,"sigma_arcsec":sigma,"flux_z_nanomaggies":brightness[i],
                         "field_galaxies_as_bright":n,"pcc":p,"shifted_at_least_as_good":k,
                         "shifted_total":len(blank),"empirical_fraction":(k+1)/(len(blank)+1),
                         "number_of_galaxy_alternatives":len(candidates)})
    except Exception as e:
        rec.update(status="query or catalogue error",error=str(e))
    return rec,rows


def main(out):
    inputs=pd.read_csv(out/"input.csv",low_memory=False)
    audit=[]
    alt=[]
    # Repeat slew epochs at identical refined positions reuse the field.
    coords=inputs.groupby(["ra","dec","poserr_arcsec"],sort=False)
    reps=[g.iloc[0].to_dict() for _,g in coords]
    groups={g.iloc[0].detection_uid:g.detection_uid.tolist() for _,g in coords}
    with ThreadPoolExecutor(max_workers=3) as pool:
        for rec,rows in pool.map(measure,reps):
            rec["applies_to"]=groups[rec["detection_uid"]]
            audit.append(rec)
            for sid in rec["applies_to"]:
                alt.extend([dict(r,detection_uid=sid) for r in rows])
    dest=out/"published_hosts"
    dest.mkdir(exist_ok=True)
    pd.DataFrame(alt,columns=["detection_uid","rank","ls_id","type","host_ra","host_dec","sep_arcsec","sigma_arcsec","flux_z_nanomaggies","field_galaxies_as_bright","pcc","shifted_at_least_as_good","shifted_total","empirical_fraction","number_of_galaxy_alternatives"]).to_csv(dest/"local_galaxies.csv",index=False)
    (dest/"local_fields.jsonl").write_text("".join(json.dumps(r)+"\n" for r in audit))
    print("Local optical fields:",len(audit),"positions,",len(alt),"per-detection alternatives",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    main(p.parse_args().out)
