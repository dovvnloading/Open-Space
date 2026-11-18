"""Join by detector run/event, supersede seven V4 positions and test four controls."""
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from contours import contains, prepare

ROOT=Path(__file__).resolve().parent


def vectors(ra,dec):
    return np.stack([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)],axis=-1)


def polygon_area(xy):
    v=vectors(*xy.T)
    center=np.sum(v[:-1],axis=0)
    center/=np.linalg.norm(center)
    left,right=v[:-1],v[1:]
    determinant=np.cross(left,right)@center
    denom=1+left@center+np.sum(left*right,axis=1)+right@center
    return float(abs(np.sum(2*np.arctan2(determinant,denom)))*(180/np.pi)**2)


def polygon_contains(xy,ra,dec):
    p=vectors(np.radians(ra),np.radians(dec))
    v=vectors(*xy.T)
    # Restrict to the hemisphere covering this small, contiguous contour;
    # otherwise the antipodal point can also have a winding number of magnitude one.
    if np.mean(v[:-1],axis=0)@p<=0:
        return False
    tangent=v-(v@p)[:,None]*p
    tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    a,b=tangent[:-1],tangent[1:]
    return bool(abs(np.arctan2(np.cross(a,b)@p,np.sum(a*b,axis=1)).sum())>np.pi)


def main():
    maps=[json.loads(row) for row in (ROOT/"icecat1_likelihood90.jsonl").read_text().splitlines()]
    revised={r["alert_id"]:np.array(r["vertices_rad"]) for r in map(json.loads,(ROOT/"icecat2_published90.jsonl").read_text().splitlines())}
    source=list(csv.DictReader((ROOT/"icecat2_v1.tsv").open(),delimiter="\t"))
    new_by_key={r["RUNID"]+"_"+r["EVENTID"]:r for r in source}
    output=[]
    for m in maps:
        if m["cosmic_ray_veto"]:
            continue
        k=m["run_event"]
        r=new_by_key.get(k)
        if r is None:
            output.append(dict(m,stream=m["stream"].split("-")[-1].upper(),localization="IceCat-1 V4 likelihood pixels"))
            continue
        assert m['alert_id']==r["NAME"] and not m["cosmic_ray_veto"]
        verts=revised[r["NAME"]]
        assert np.all(np.isfinite(verts)) and np.array_equal(verts[0],verts[-1])
        output.append(dict(run_event=k,alert_id=r["NAME"],utc=m["utc"],stream=r["I3TYPE"].split("-")[-1].upper(),
                           ra_deg=float(r["RA [deg]"]),dec_deg=float(r["DEC [deg]"]),area_deg2=polygon_area(verts),
                           published_area_deg2=float(r["CONTOUR AREA (90%) [sqdeg]"]),polygons_deg=[np.degrees(verts).tolist()],
                           localization="IceCat-2 preliminary V1 publisher 90% contour"))
    assert len(output)==340 and len({r["run_event"] for r in output})==340 and sum("polygons_deg" in r for r in output)==7
    data="".join(json.dumps(m,sort_keys=True,separators=(",",":"))+"\n" for m in output).encode()
    (ROOT/"archival_maps.jsonl.gz").write_bytes(gzip.compress(data,compresslevel=9,mtime=0))
    controls=[]
    originals={m["alert_id"]:m for m in maps}
    for fn in ["historical_optical_pairs.csv","historical_txs.csv"]:
        for r in csv.DictReader((ROOT.parent/"controls"/fn).open()):
            ra,dec=float(r["ra"]),float(r["dec"])
            m=originals[r["alert"]];p=prepare(m["uniq90"])
            controls.append(dict(source=r["source"],alert_id=r["alert"],source_ra_deg=ra,source_dec_deg=dec,
                                 icecat1_v4_likelihood90=bool(contains(p,ra,dec)),icecat2_prelim_v1_contour90=polygon_contains(revised[r["alert"]],ra,dec)))
    with (ROOT/"native_control_membership.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(controls[0]),lineterminator="\n");w.writeheader();w.writerows(controls)


if __name__=="__main__":
    main()
