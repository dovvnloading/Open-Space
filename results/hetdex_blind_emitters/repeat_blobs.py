"""Repeat-geometry search for Ly-alpha bridges in frozen HETDEX two-knot proposals."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from spectral import BAD_AIR, SKY_AIR
from vet_cubes import cube_test
from vetting import Ancillary, Cache, atom, chord, enrich, external_check, number, separation, tap_batch, validate_local, hsc_photometry


SKY_TOLERANCE = 1.25  # arcsec, applied independently to both observed knots
WAVE_TOLERANCE = 3.5  # air Angstrom, applied independently to both observed lines
KNOT_KEYS = """release row_index source_id detectid shotid ra dec ra_det dec_det ifuslot fwhm wave snr flux flux_err
knot_detectid secondary_detectid secondary_ra secondary_dec secondary_wave secondary_snr knot_separation_arcsec
z_lya source_type agn_flag z_hetdex z_hetdex_src z_hetdex_conf flag_best sigma_A spectral_EW0_A spectral_EW0_lower_A""".split()


def dumps(obj):
    return json.dumps(atom(obj),separators=(",",":"),allow_nan=False)


def coords(p):
    return ((p["ra"],p["dec"]),(p["secondary_ra"],p["secondary_dec"]))


def wavelength(p):
    return (p["wave"],p["secondary_wave"])


def source_rows(path,counts):
    rows = []
    with gzip.open(path,"rt") as f:
        for line,text in enumerate(f,1):
            r = json.loads(text)
            p,v = r["peak"],r["vet"]
            if p["channel"] != "blob" or p.get("secondary_detectid") is None:
                continue
            counts["original_two_knot_proposals"] += 1
            if not all(np.isfinite(number(x)) for pair in coords(p) for x in pair):
                counts["invalid_endpoint_coordinates"] += 1
                continue
            # Preserve the prior decision even for excluded/published rows so a
            # previously identified object cannot enter through another visit.
            prior = {k:v[k] for k in ("disposition","reasons","LAN","cube","passed","low_z_families") if k in v}
            rec = {k:p[k] for k in KNOT_KEYS if k in p}
            rec["channel"] = "blob"
            rows.append((line,rec,prior))
    return rows


def geometric_groups(rows,counts):
    if not rows:
        return []
    w = np.array([wavelength(p) for _,p,_ in rows],dtype=float)
    shot = np.array([p["shotid"] for _,p,_ in rows])
    pos = np.array([coords(p) for _,p,_ in rows],dtype=float)
    xyz = chord(pos[:,:,0],pos[:,:,1])
    mid = xyz.sum(axis=1)
    mid /= np.linalg.norm(mid,axis=1)[:,None]
    tol = 2*np.sin(np.deg2rad(SKY_TOLERANCE/3600)/2)
    tree=cKDTree(mid)
    edges=[]
    parent={}
    def root(x):
        while parent.setdefault(x,x)!=x:
            parent[x]=parent[parent[x]]
            x=parent[x]
        return x
    # A bright single pointing can generate many nearby pair midpoints. Discard
    # same-shot neighbors before building any temporary coordinate-pair arrays.
    for start in range(0,len(rows),512):
        neighbors=tree.query_ball_point(mid[start:start+512],tol*1.000001)
        for local,near in enumerate(neighbors):
            a=start+local
            j=np.asarray(near,dtype=np.int64)
            j=j[(j>a)&(shot[j]!=shot[a])]
            if not len(j):
                continue
            counts["cross_shot_midpoint_coincidences"]+=len(j)
            dist0=np.linalg.norm(xyz[a]-xyz[j],axis=2)
            dist1=np.linalg.norm(xyz[a]-xyz[j,::-1],axis=2)
            delta0=w[a]-w[j]
            delta1=w[a]-w[j,::-1]
            direct=dist0.max(axis=1)<=tol
            reverse=dist1.max(axis=1)<=tol
            for delta in (-50.,50.):
                key="displaced_minus50A_geometric_edges" if delta<0 else "displaced_plus50A_geometric_edges"
                matches=(direct & (abs(delta0-delta).max(axis=1)<=WAVE_TOLERANCE)) | (reverse & (abs(delta1-delta).max(axis=1)<=WAVE_TOLERANCE))
                counts[key]+=int(matches.sum())
            same=direct & (abs(delta0).max(axis=1)<=WAVE_TOLERANCE)
            swap=reverse & (abs(delta1).max(axis=1)<=WAVE_TOLERANCE)
            for k in np.flatnonzero(same|swap):
                b=int(j[k])
                flipped=bool(swap[k] and (not same[k] or dist1[k].max()<dist0[k].max()))
                ds=dist1[k] if flipped else dist0[k]
                dw=delta1[k] if flipped else delta0[k]
                edges.append((a,b,flipped,np.rad2deg(2*np.arcsin(ds/2))*3600,dw))
                ra,rb=root(a),root(b)
                parent[max(ra,rb)]=min(ra,rb)
    grouped = defaultdict(list)
    for e in edges:
        grouped[root(e[0])].append(e)
    counts.update(geometry_consistent_edges=len(edges),geometry_components=len(grouped),
                  geometry_matched_rows=len(parent))
    return list(grouped.values())


def lan_matches(p,aux):
    if "LAN" not in aux.tables:
        return []
    _,tab,idx,tree = aux.tables["LAN"]
    xyz = chord(*np.array(coords(p)).T)
    mid = xyz.sum(axis=0)
    mid /= np.linalg.norm(mid)
    probes = np.vstack((xyz,mid))
    rows = set()
    for vec in probes:
        rows.update(tree.query_ball_point(vec,2*np.sin(np.deg2rad(25./3600)/2)))
    result=[]
    for j in rows:
        r=tab[int(idx[j])]
        z=number(r["z_hetdex"])
        if abs(z-p["z_lya"])>.015:
            continue
        direction=chord(float(r["ra"]),float(r["dec"]))
        separation_arcsec=np.rad2deg(2*np.arcsin(np.linalg.norm(probes-direction,axis=1)/2))*3600.
        # Explicit angular units for the isophotal area. Uncertain units give only the positional test.
        area=number(r["area_r_iso_circ"],0.)
        unit=str(tab.columns["area_r_iso_circ"].unit or "").lower()
        radius=math.sqrt(area/math.pi) if area>0 and ("arcsec" in unit or "asec" in unit) else 0.
        if separation_arcsec.min()<=3. or (number(r["flag_resolved"])==1 and separation_arcsec.min()<=radius):
            result.append(dict(detectid=int(r["detectid"]),z=z,flag_resolved=number(r["flag_resolved"],None),
                               distances_arcsec=separation_arcsec.tolist(),isophotal_radius_arcsec=radius))
    return result


def observation_issue(p,v):
    if v["disposition"] in ("excluded","published","candidate"):
        return "first_pass_"+v["disposition"]
    for w in wavelength(p):
        if any(abs(w-s)<=6 for s in SKY_AIR) or any(lo<=w<=hi for lo,hi in BAD_AIR):
            return "catalog_sky_or_calibration_window"
    if min(number(p.get("snr")),number(p.get("secondary_snr")))<7:
        return "knot_below_seven_sigma"
    width=number(p.get("fwhm"),0.)
    # A bridge must have a nonzero length after masking 1.3 seeing FWHM around each knot.
    if .5<width<5 and p["knot_separation_arcsec"]<=2.6*width:
        return "no_bridge_between_seeing_masks"
    if str(p.get("source_type","")).strip().lower()=="agn" or number(p.get("agn_flag")) in (0,1):
        return "hetdex_agn_or_broad_line_source"
    return None


def safe_cube(p,prior,cache):
    c=prior.get("cube") or {}
    if c.get("status")=="measured" and c.get("bridge") is not None:
        return dict(c,measurement_source="first-pass frozen cube")
    try:
        return dict(cube_test(p,cache.fetch),measurement_source="repeat-geometry cube")
    except Exception as e:
        return dict(status="unavailable",error=str(e))


def measurable(c):
    return all(number((c.get(k) or {}).get("sn"))>-math.inf for k in
               ("central_aperture","bridge","background_wavelength_test"))


def clean_bridge(c):
    return (c.get("central_above_5sigma") and c.get("bridge_above_5sigma") and
            not c.get("background_suspect") and number((c.get("background_wavelength_test") or {}).get("sn"))>-math.inf)


def knot_catalog(cat, ident):
    k=np.searchsorted(cat.sorted_det_ids,ident)
    if k==len(cat.sorted_det_ids) or int(cat.sorted_det_ids[k])!=ident:
        return {}
    d=cat.detinfo[int(cat.det_order[k])]
    pick={key:atom(d[cat.det_columns[key]]) for key in
         ("detectid","shotid","source_id","ra_det","dec_det","wave","sn","flux","flux_err","fiber_id","multiframe",
          "x_raw","y_raw","flag_best","flag_badamp","flag_pixmask","flag_badfib","flag_baddet") if key in cat.det_columns}
    sid=pick.get("source_id")
    if sid is None:
        return pick
    a,b=np.searchsorted(cat.sorted_src_ids,sid,side="left"),np.searchsorted(cat.sorted_src_ids,sid,side="right")
    rr=cat.src_order[a:b]
    rr=rr[cat.shot_ids[rr]==pick["shotid"]]
    if not len(rr):
        return pick
    meta=cat.metadata(int(rr[0]))
    pick["associated_row"]=int(rr[0])
    for key in ("source_type","z_hetdex","z_hetdex_src","z_hetdex_conf","agn_flag","fwhm"):
        pick[key]=meta.get(key)
    if meta.get("ra_det") is not None and meta.get("dec_det") is not None:
        pick["extraction_offset_arcsec"]=float(separation(pick["ra_det"],pick["dec_det"],meta["ra_det"],meta["dec_det"]))
    return pick


def adjudicate_edge(cat,aux,cache,rows,a,b,swapped,cubes):
    i,p,_=rows[a]
    _,q,_=rows[b]
    # Primary direction must be represented by a true source-centred extraction.
    pk,sk=[knot_catalog(cat,int(p[k])) for k in ("knot_detectid","secondary_detectid")]
    qk,tk=[knot_catalog(cat,int(q[k])) for k in ("knot_detectid","secondary_detectid")]
    if swapped:
        qk,tk=tk,qk
    knots=[pk,sk,qk,tk]
    absences=[]
    hard=[]
    for knot in knots:
        if not knot or knot.get("associated_row") is None:
            absences.append("knot classification missing")
        if str(knot.get("source_type")).strip().lower()=="agn" or number(knot.get("agn_flag")) in (0,1):
            hard.append("a repeated knot belongs to an HETDEX AGN or broad-line candidate")
        label=str(knot.get("source_type","")).strip().lower()
        ref=str(knot.get("z_hetdex_src","")).lower()
        if label in ("oii","lzg","star") and number(knot.get("z_hetdex_conf"))>=.8 and any(
                t in ref for t in ("sdss","desi","diagnos","specz","2em")):
            if number(knot.get("z_hetdex"))>=0 and abs(number(knot["z_hetdex"])-p["z_lya"])>.02:
                hard.append("a knot has a secure foreground classification")
        if number(knot.get("flag_best"),1)==0:
            absences.append("a knot detection failed its HETDEX quality flag")
    for x,y in ((pk,qk),(sk,tk)):
        different=(x.get("multiframe") and y.get("multiframe") and x["multiframe"]!=y["multiframe"]) or (
            number(x.get("y_raw"))>-math.inf and number(y.get("y_raw"))>-math.inf and abs(x["y_raw"]-y["y_raw"])>=4)
        if not different:
            absences.append("persistent detector/fibre response not excluded at both knots")
    if number(pk.get("extraction_offset_arcsec"),math.inf)>1.25:
        absences.append("1D spectral extraction is displaced from primary knot")
    p={**p,"catalog":cat.metadata(int(p["row_index"]))}
    p.update(enrich(p,cat))
    prev = getattr(cat,"confirmed_repeat_rows",None)
    idx=pk.get("associated_row")
    alt=qk.get("associated_row")
    cat.confirmed_repeat_rows={p["row_index"]:(alt,) if idx==p["row_index"] and alt is not None and
                               number(qk.get("extraction_offset_arcsec"),math.inf)<=1.25 else ()}
    try:
        v=validate_local(p,cat,aux)
    finally:
        if prev is None:
            del cat.confirmed_repeat_rows
        else:
            cat.confirmed_repeat_rows=prev
    if v["disposition"]=="unresolved":
        if "LAN" in aux.errors:
            v["LAN_catalog_unavailable"]=True
        qa = dict(q)
        if swapped:
            qa.update(ra=q["secondary_ra"],dec=q["secondary_dec"],ra_det=q["secondary_ra"],dec_det=q["secondary_dec"],
                      secondary_ra=q["ra"],secondary_dec=q["dec"],wave=q["secondary_wave"])
        # Supply the actual independently matched astrometric geometry, not the original sky positions.
        v["repeats"]=[dict(shotid=q["shotid"],row_index=alt,snr=q["secondary_snr"] if swapped else q["snr"],
                 flux=qk.get("flux")*1e-17 if qk.get("flux") is not None else None,
                 repeat_cube_meta=qa,measurement="independently catalogued knot in matched geometry")]
        if "independent HETDEX shot" not in v["passed"]:
            v["passed"].append("independent HETDEX shot")
        p["rescue_audit"]={"cube":cubes[a]}
        from focus import batched_matches, trace_spectra
        found=[dict(peak=k,first={},same_visit_two_line_repeats=[]) for k in (p,swapped_view(p))]
        try:
            bypoint,cil=batched_matches(found)
            sil=trace_spectra(found,bypoint,cache.root/"repeat-pair-spectra")
            p["rescue_audit"]["SPARCL"]=found[0]["external"]
            v["secondary_spectroscopy"]=found[1]["external"]
            v["knot_spectroscopy_receipts"]=dict(catalog=cil,spectra=sil)
            for m in found[1]["external"].get("matches",[]):
                if number(m.get("redshift_warning"))!=0:
                    continue
                if str(m.get("spectype","")).upper() in ("QSO","STAR"):
                    hard.append("independent spectroscopy classifies the second knot as QSO or star")
                zz=number(m.get("redshift"))
                if zz>=0 and abs(zz-p["z_lya"])>.02:
                    absences.append("independent spectral redshift at secondary knot disagrees")
        except Exception as e:
            v["secondary_spectroscopy"]=dict(status="unavailable",error=str(e))
        try:
            v["secondary_HSC"]=hsc_photometry(swapped_view(p))
        except Exception as e:
            v["secondary_HSC"]=dict(status="unavailable",error=str(e))
        points=[p,{"ra":p["secondary_ra"],"dec":p["secondary_dec"]}]
        with ThreadPoolExecutor(max_workers=3) as pool:
            replies=dict(zip(("NED","SIMBAD","Legacy"),pool.map(lambda s:tap_batch(s,points),("NED","SIMBAD","Legacy"))))
        for name,grp in replies.items():
            v[name]=grp[0]
        v=external_check(p,v,cache)
        v["secondary_archive_cones"]={name:grp[1] for name,grp in replies.items()}
        for name,grp in replies.items():
            ans=grp[1]
            if ans.get("status")!="queried":
                absences.append(name+" query unavailable at secondary knot")
            for m in ans.get("matches",[]):
                if number(m.get("sep_arcsec"),math.inf)>1.25:
                    continue
                ty=str(m.get("pretype",m.get("otype",""))).strip().upper()
                if ty in ("QSO","AGN","SY1","SY2","SEYFERT","BLL","BLAZAR","BLA"):
                    hard.append("catalogued active nucleus at secondary knot")
                zz=number(m.get("z",m.get("rvz_redshift")))
                if zz>=0 and abs(zz-p["z_lya"])>.02 and (m.get("rvz_bibcode") or m.get("refcode")):
                    absences.append("secondary counterpart redshift disagrees; association needs review")
    if hard:
        v["disposition"]="excluded"
    elif absences and v["disposition"]=="candidate":
        v["disposition"]="unresolved"
    v.setdefault("reasons",[]).extend(sorted(set(hard+absences)))
    v.update(confirmed_knot_catalogs=knots,matched_second_shot_original_line=rows[b][0])
    p.pop("catalog",None)
    p.pop("rescue_audit",None)
    return {"first_pass_line":i,"peak":p,"vet":v}


def swapped_view(p):
    return {**p,"ra":p["secondary_ra"],"dec":p["secondary_dec"],"ra_det":p["secondary_ra"],
            "dec_det":p["secondary_dec"],"wave":p["secondary_wave"],"secondary_ra":p["ra"],"secondary_dec":p["dec"]}


def find_repeated_blobs(indir,out,release,cache_root):
    indir,out,cache_root = map(Path,(indir,out,cache_root))
    out.mkdir(exist_ok=True,parents=True)
    path=indir/"vetted.jsonl.gz"
    with path.open("rb") as f:
        digest=hashlib.file_digest(f,"sha256").hexdigest()
    counts=Counter()
    rows=source_rows(path,counts)
    groups=geometric_groups(rows,counts)
    cache=Cache(cache_root)
    failures={}
    regular_fetch=cache.fetch
    def cached_fetch(kind,name,urls):
        if name in failures:
            raise OSError(failures[name])
        try:
            return regular_fetch(kind,name,urls)
        except OSError as e:
            failures[name]=str(e)
            raise
    cache.fetch=cached_fetch
    stages=Counter()
    with ExitStack() as files:
        ledger=files.enter_context(gzip.open(out/"systems.jsonl.gz","wt"))
        overlay=files.enter_context(gzip.open(out/"vetted.jsonl.gz","wt"))
        aux=None
        cat=None
        restoration_error=None
        for gi,edges in enumerate(groups):
            ids=sorted(set(x for a,b,*_ in edges for x in (a,b)))
            obj={"group":gi,"release":release,
                 "geometry_edges":[dict(first_pass_lines=[rows[a][0],rows[b][0]],swapped=bool(sw),
                     endpoint_separation_arcsec=ds,endpoint_wavelength_deltas_A=dw) for a,b,sw,ds,dw in edges],
                 "observations":[]}
            notes,problems={},{}
            for j in ids:
                line,p,v=rows[j]
                issues=observation_issue(p,v)
                notes[j]=dict(first_pass_line=line,peak=p,original_vetting=v,local_gate=issues)
                obj["observations"].append(notes[j])
                if issues:
                    problems[j]=issues
            candidate_edges=[e for e in edges if e[0] not in problems and e[1] not in problems]
            if candidate_edges and any(rows[j][2].get("disposition") in ("published","candidate") for j in ids):
                candidate_edges=[]
                stage="previously catalogued or first-pass candidate in geometry group"
            else:
                stage="no clean independent geometry after spectral and resolution cuts"
            if candidate_edges:
                if aux is None:
                    aux=Ancillary(cache)
                    files.callback(aux.close)
                for j in ids:
                    notes[j]["LAN_spatial_matches"]=lan_matches(rows[j][1],aux)
                if any(any(m.get("flag_resolved")==1 for m in n["LAN_spatial_matches"]) or
                       number((n["original_vetting"].get("LAN") or {}).get("flag_resolved"))==1
                       for n in notes.values()):
                    candidate_edges=[]
                    stage="previously resolved HETDEX LAN at knot, bridge midpoint or published isophote"
            if candidate_edges:
                counts["systems_requiring_two_shot_cube_check"]+=1
                scans=sorted(set(x for a,b,*_ in candidate_edges for x in (a,b)),
                             key=lambda j:(rows[j][1]["shotid"],str(rows[j][1].get("ifuslot")),j))
                cubes={}
                coverage=Counter()
                for j in scans:
                    _,p,v=rows[j]
                    c=safe_cube(p,v,cache)
                    notes[j]["cube"]=c
                    cubes[j]=c
                    if (c.get("bridge") or {}).get("sn") is not None:
                        coverage["bridge_measured"]+=1
                    if clean_bridge(c):
                        coverage["clean_primary_and_bridge"]+=1
                obj["cube_observations"]=dict(coverage)
                counts.update(coverage)
                available=[e for e in candidate_edges if all(measurable(cubes[j]) for j in e[:2])]
                counts["systems_with_measurable_bridge_in_distinct_shots"]+=bool(available)
                best=[e for e in candidate_edges if all(clean_bridge(cubes[j]) for j in e[:2])]
                stage=("measured bridge below joint central, sky and two-shot threshold" if len(available)==len(candidate_edges)
                       else "incomplete independent cube coverage or sky measurements")
                if best:
                    counts["systems_with_primary_and_bridge_confirmed_twice"]+=1
                    # Confirm the fainter knot in each shot and the bridge centered on
                    # its own independently fitted line, using the same empirical control windows.
                    needed=set(x for a,b,*_ in best for x in (a,b))
                    for j in sorted(needed):
                        rev=safe_cube(swapped_view(rows[j][1]),{},cache)
                        notes[j]["secondary_centered_cube"]=rev
                    fully_measured=all(measurable(notes[j]["secondary_centered_cube"]) for j in needed)
                    best=[e for e in best if all(clean_bridge(notes[j]["secondary_centered_cube"]) for j in e[:2])]
                    stage=("measured secondary-centered trial below two-shot threshold" if fully_measured else
                           "incomplete secondary-centered cube coverage or background")
                if best:
                    counts["systems_with_both_knots_and_bridge_confirmed_twice"]+=1
                    # Select one actual matched edge per connected system; endpoints of
                    # every edge tested directly, never inferred through a graph path.
                    best.sort(key=lambda e:(-min((cubes[j]["bridge"] or {})["sn"] for j in e[:2]),e[0],e[1]))
                    a,b,sw,*_=best[0]
                    if cat is None and restoration_error is None:
                        from acquisition import Spectra, prepare
                        try:
                            cat=files.enter_context(Spectra(prepare(release,cache_root,out),release))
                        except Exception as e:
                            restoration_error=str(e)
                    if cat is None:
                        stage="unresolved; release spectral restoration failed"
                        obj["restoration_error"]=restoration_error
                    else:
                        try:
                            result=adjudicate_edge(cat,aux,cache,rows,a,b,sw,cubes)
                            result.update(repeat_blob_group=gi,geometry_lines=[rows[x][0] for x in ids],
                                first_pass_disposition=rows[a][2]["disposition"],
                                cube_evidence_ledger="systems.jsonl.gz",geometry_edge_lines=[rows[a][0],rows[b][0]])
                            overlay.write(dumps(result)+"\n")
                            overlay.flush()
                            stage="source_vetting_"+result["vet"]["disposition"]
                        except Exception as e:
                            stage="unresolved; source vetting failed"
                            obj["vetting_error"]=str(e)
            obj["status"]=stage
            stages[stage]+=1
            ledger.write(dumps(obj)+"\n")
            if (gi+1)%50==0:
                print(release,"repeat geometry systems",gi+1,"of",len(groups),dict(stages),flush=True)
        cat_errors=dict(aux.errors) if aux is not None else {}
    summary={"release":release,"input":"../vetted.jsonl.gz","input_sha256":digest,"counts":dict(counts),
        "system_statuses":dict(stages),"ancillary_catalog_errors":cat_errors,"failed_remote_inputs":failures,
        "scope":"Cohorts are connected components of directly matching catalog knot pairs in different shots of the same release. Each accepted pair has both endpoints within 1.25 arcsec and both centroids within 3.5 air Angstrom, permitting reversed brightness order. A prior resolved LAN association excludes that system from novelty follow-up. Every distinct retained geometry is tested; no rank cap. Incomplete cube or archive coverage remains unresolved.",
        "geometric_control":"Hold both sky coordinates fixed and shift both knot wavelengths in the second shot by -50 or +50 A. These counts are positional/line coincidences before source or cube vetting.",
        "spatial_tests":"Each independent shot needs a positive central aperture and two-Moffat-subtracted bridge at five empirical sigma, clean outer background, and the same requirements centered at the secondary knot. A nominal PSF residual is conditional on the observed seeing and adopted PSF. The experiment does not select single-center diffuse halos or objects missing from catalog knots.",
        "products":["systems.jsonl.gz","vetted.jsonl.gz","inputs.json"]}
    (out/"summary.json").write_text(json.dumps(atom(summary),indent=2,allow_nan=False)+"\n")
    (out/"inputs.json").write_text(dumps(cache.inputs)+"\n")
    return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--release",choices=("sc1","sc2"),required=True)
    parser.add_argument("--cache",type=Path,required=True)
    args=parser.parse_args()
    find_repeated_blobs(args.input,args.output,args.release,args.cache)
