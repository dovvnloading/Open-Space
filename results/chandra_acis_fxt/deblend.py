"""Compare two fixed, nonoverlapping Chandra apertures at original trial and archived CSC positions."""
import argparse
import gzip
import hashlib
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.optimize import brentq
from scipy.stats import fisher_exact

import campaign
import diagnostics
import followup
import ingest_primary as ingest
import pointed_visits
from archive import dump
from migrate import INVENTORY
from report import bounds
from search import _orientation, _sky_wcs
from target_history import CSC
from vet import separation
from visits import _inside_segment

PRIOR = "45011722864-1"
PHASES = ("before", "trigger", "after")
# Original Git blobs of the complete CSC request/response receipts, including URLs and entity SHA-256s.
CSCS = {
    22929: ("target-00", "2CXO J114044.9-263004", "55a7e569670ea81b471a097e02715023aaabc04e",
            "cd16a56dc6db330bd582d8dfb1a43201130bb48f"),
    23779: ("target-01", "2CXO J125704.0+275119", "2e8e93680c79c269ec899e33846195b70db1b9f1",
            "ea1a6ccae7faeb2c52a0b325dce46c0de972ba26"),
}
need = ingest.require
sha = ingest.sha


def csc_rows(root, oid):
    target, name, master_blob, assoc_blob = CSCS[oid]
    area = root/"pointed_visits"/PRIOR/"history"/target
    tables, proofs = {}, {}
    for kind, blob in (("master", master_blob), ("detections", assoc_blob)):
        path = area/f"csc_{kind}.json"
        raw = path.read_bytes()
        need(hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest() == blob, "Changed CSC request receipt")
        rec = campaign.read(path)
        response = rec["response"]
        gz = area/f"csc_{kind}.xml.gz"
        need(rec["status"] == "ok" and rec["endpoint"] == CSC and rec["parameters"]["version"] == "rel2.2" and
             response["raw_entity_gzip"] == gz.name and sha(gz) == response["gzip_sha256"], "Changed CSC archive")
        xml = gzip.decompress(gz.read_bytes())
        need(len(xml) == response["bytes"] and hashlib.sha256(xml).hexdigest() == response["sha256"], "Changed CSC bytes")
        doc = ET.fromstring(xml)
        need(all(v.get("value","").upper()=="OK" for v in doc.findall(".//{*}INFO")
                 if v.get("name","").upper()=="QUERY_STATUS"), "Incomplete CSC response")
        table, = doc.findall(".//{*}TABLE")
        cols = [(e.get("name") or e.get("ID")).lower().split(".")[-1] for e in table.findall("{*}FIELD")]
        rows = [dict(zip(cols, ((e.text or "").strip() for e in row), strict=True)) for row in table.findall(".//{*}TR")]
        need(rows == rec["rows"] and len(rows) == rec["row_count"], "CSC row replay differs")
        tables[kind] = [r for r in rows if r["name"] == name]
        proofs[kind] = {"receipt": path.relative_to(root).as_posix(), "git_blob": blob, "receipt_sha256": sha(path),
                        "archived_body": gz.relative_to(root).as_posix(), **response}
    master, = tables["master"]
    same = [r for r in tables["detections"] if int(r["obsid"]) == oid]
    need(len(same)>0, "Prior same-observation CSC association missing")
    return master, same, proofs


def nearest_ellipse(east, north, major, minor, pa):
    """Distance in the CSC tangent plane from a point to the filled error ellipse."""
    need(0 < minor <= major and np.isfinite([east,north,major,minor,pa]).all(), "Invalid CSC error ellipse")
    phi = np.radians(pa)
    coord = np.array([east*np.sin(phi)+north*np.cos(phi), east*np.cos(phi)-north*np.sin(phi)])
    axes = np.array([major, minor])
    if np.sum((coord/axes)**2) <= 1:
        return 0.
    def residual(lam):
        return np.sum((axes*coord/(lam+axes**2))**2)-1
    limit = major*float(np.linalg.norm(coord))
    lam = brentq(residual, 0., limit+major**2)
    return float(np.linalg.norm(coord*lam/(lam+axes**2)))


def association(master, sky, radius):
    csc = SkyCoord(float(master["ra"])*u.deg, float(master["dec"])*u.deg)
    center = SkyCoord(float(sky[0])*u.deg, float(sky[1])*u.deg)
    east, north = csc.spherical_offsets_to(center)
    d = nearest_ellipse(east.arcsec, north.arcsec, float(master["err_ellipse_r0"]), float(master["err_ellipse_r1"]),
                        float(master["err_ellipse_ang"]))
    sep = float(csc.separation(center).arcsec)
    return {"center_icrs_deg": list(map(float, sky)), "radius_arcsec": float(radius), "separation_arcsec": sep,
            "center_in_csc_ellipse": d==0, "csc_center_in_circle": sep<=radius,
            "distance_to_filled_csc_ellipse_arcsec": d, "csc_ellipse_intersects_circle": d<=radius}


def geometry(wcs, event, master):
    trial = np.array([event["aperture_x"], event["aperture_y"]], dtype=float)
    sky = wcs.all_pix2world([trial], 1)[0]
    position = [float(master["ra"]), float(master["dec"])]
    csc = wcs.all_world2pix([position], 1)[0]
    need(np.isfinite(np.r_[trial,csc,sky]).all() and separation(*position,*wcs.all_pix2world([csc],1)[0])<1e-4,
         "Invalid or non-invertible input WCS")
    d = float(np.linalg.norm(csc-trial))
    original = event["aperture_radius_px"]
    radius = min(original, .45*d)
    need(0<radius<d/2 and original>0, "Centers do not admit disjoint apertures")
    mid = (trial+csc)/2
    angles = np.arange(72)*(2*np.pi/72)
    boundary = wcs.all_pix2world(trial+original*np.column_stack([np.cos(angles),np.sin(angles)]), 1)
    scale = np.array([separation(*sky,*p) for p in boundary])/original
    lo, hi = map(float, [scale.min(), scale.max()])
    extraction = association(master, sky, lo*original)
    extraction.update(radius_range_arcsec=[lo*original, hi*original],
                      csc_center_in_circle=bool(d<=original),
                      ellipse_intersection_over_scale_range=[extraction["csc_ellipse_intersects_circle"],
                                                             extraction["distance_to_filled_csc_ellipse_arcsec"]<=hi*original])
    return {"wcs_origin": 1, "event_wcs_header": wcs.to_header().tostring(sep="\n",endcard=False,padding=False),
            "trial_xy": trial.tolist(), "csc_xy": csc.tolist(), "midpoint_xy": mid.tolist(),
            "separation_px": d, "core_radius_px": radius, "original_radius_px": original,
            "annulus_px": [d/2+2*radius,d/2+4*radius], "arcsec_per_sky_px_range": [lo,hi],
            "extraction_circle": extraction,
            "localization_circle": association(master,[event["ra_deg"],event["dec_deg"]],event["position_error_arcsec"])}


def fraction(k, n):
    return {"k": int(k), "n": int(n), "fraction": k/n if n else None, "cp95": list(bounds(k,n))}


def describe(values):
    if not len(values):
        return None
    return {"min_px": float(min(values)), "median_px": float(np.median(values)), "max_px": float(max(values))}


def measure(data, geom, event, lo, hi):
    by_phase, tags = {}, []
    x, y = data["x"], data["y"]
    def distance(key):
        cx, cy = geom[key]
        return np.hypot(x-cx,y-cy)
    dt, dc, dm = (distance(key) for key in ("trial_xy","csc_xy","midpoint_xy"))
    r = geom["core_radius_px"]
    rin, rout = geom["annulus_px"]
    src, csc = dt<=r, dc<=r
    ann = (dm>=rin)&(dm<=rout)
    original = dt<=geom["original_radius_px"]
    need(not np.any((src&csc)|(src&ann)|(csc&ann)), "Apertures or local background overlap")
    scale = r*r/(rout*rout-rin*rin)
    matrix = _orientation(data)
    limits = ([int(data[k].min()) for k in ("cx","cy")], [int(data[k].max()) for k in ("cx","cy")])
    duty = float(event["live_fraction"])
    need(0<duty<=1 and [lo,hi]==event["diagnostics"]["gti_met_s"], "Invalid CCD time support")
    for phase, left, right in zip(PHASES, (lo,event["start_met_s"],event["stop_met_s"]),
                                  (event["start_met_s"],event["stop_met_s"],hi), strict=True):
        time = (data["t"]>=left)&(data["t"]<right)
        a,b,n,o = (int(np.sum(time & mask)) for mask in (src,csc,ann,original))
        need(o==event[{"before":"before_counts","trigger":"on_counts","after":"after_counts"}[phase]],
             "Cannot reproduce the original photometry")
        # The independent-visit support guard uses a 4-radius bounding disk; choose it to contain our entire annulus.
        cover = _inside_segment(data[time], *geom["midpoint_xy"], rout/4, left, right, duty, matrix, limits)
        accepted = np.zeros(len(data), dtype=bool)
        for t0,t1 in cover["accepted_bins"]:
            accepted |= (data["t"]>=t0)&(data["t"]<t1)
        by_phase[phase] = {"met_s":[left,right],"gt_live_s":(right-left)*duty, "trial_core":a, "csc_core":b,
            "shared_annulus":n, "nominal_background_per_core":scale*n, "original_aperture":o,
            "original_photons_in_trial_core":int(np.sum(time&original&src)),
            "original_photons_in_csc_core":int(np.sum(time&original&csc)),
            "original_photons_outside_both_cores":int(np.sum(time&original&~src&~csc)),
            "original_photons_nearer_csc":int(np.sum(time&original&(dc<dt))),
            "trial_fraction_among_both_cores":fraction(a,a+b),
            "original_distances_to_trial":describe(dt[time&original]), "original_distances_to_csc":describe(dc[time&original]),
            "common_geometry":{"supported_live_s":sum((b-a)*duty for a,b in cover["accepted_bins"]),
                "accepted_bins":cover["accepted_bins"],"rejected_bins":cover.get("rejected_bins"),
                "reason":cover.get("reason"),"neighborhood_quadrant_counts":cover.get("neighborhood_quadrant_counts"),
                "trial_core":int(np.sum(accepted&src)),"csc_core":int(np.sum(accepted&csc)),
                "shared_annulus":int(np.sum(accepted&ann))}}
        for i in np.flatnonzero(time&(src|csc|original|ann)):
            row = data[i]
            tags.append({"events_row_zero_based":int(row["row"]),"phase":phase,"time_met_s":float(row["t"]),
                "x":float(x[i]),"y":float(y[i]),"chipx":int(row["cx"]),"chipy":int(row["cy"]),"energy_ev":float(row["e"]),
                "grade":int(row["grade"]),"original_aperture":bool(original[i]),"common_geometry_supported":bool(accepted[i]),
                "region":"trial_core" if src[i] else "csc_core" if csc[i] else "shared_annulus" if ann[i] else "original_only",
                "trial_offset_px":float(dt[i]),"csc_offset_px":float(dc[i]),"midpoint_offset_px":float(dm[i])})
    def compare(counts):
        c = [[counts["trigger"][k] for k in ("trial_core","csc_core")],
             [sum(counts[p][k] for p in ("before","after")) for k in ("trial_core","csc_core")]]
        ok = min(np.sum(c,axis=0))>0 and min(np.sum(c,axis=1))>0
        return {"trigger_then_pooled_flanks_by_trial_and_csc":c,"fisher_two_sided_descriptive_p":
                float(fisher_exact(c,alternative="two-sided").pvalue) if ok else None,
                "trigger_trial_fraction":fraction(c[0][0],sum(c[0])), "flank_trial_fraction":fraction(c[1][0],sum(c[1]))}
    return by_phase,tags,{"all_recorded":compare(by_phase),
                          "common_geometry_supported":compare({p:v["common_geometry"] for p,v in by_phase.items()})}


def picture(dest, geom, tags):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    fig, axes = plt.subplots(1,3,figsize=(10,3.6),sharex=True,sharey=True,layout="constrained")
    trial = np.array(geom["trial_xy"])
    csc = np.array(geom["csc_xy"])-trial
    mid = csc/2
    extent = max(geom["original_radius_px"]+np.linalg.norm(mid),geom["annulus_px"][1])*1.05
    for ax,phase in zip(axes,PHASES,strict=True):
        p = [r for r in tags if r["phase"]==phase]
        ax.scatter([r["x"]-trial[0] for r in p],[r["y"]-trial[1] for r in p],s=9,c="#343434",alpha=.7)
        for xy,r,color,style in (([0,0],geom["original_radius_px"],"#666666","--"),
                                ([0,0],geom["core_radius_px"],"#1864ab","-"),
                                (csc,geom["core_radius_px"],"#b14d03","-"),
                                *((mid,r,"#888888",":") for r in geom["annulus_px"])):
            ax.add_patch(Circle(xy, r,fill=False,edgecolor=color,linestyle=style,linewidth=1.3))
        ax.plot(0,0,marker="+",color="#1864ab")
        ax.plot(*csc,marker="+",color="#b14d03")
        ax.set(aspect="equal",xlim=(mid[0]-extent,mid[0]+extent),ylim=(mid[1]-extent,mid[1]+extent),
               xlabel="SKY ΔX (pixel)",title=phase.title())
    axes[0].set_ylabel("SKY ΔY (pixel)")
    fig.savefig(dest/"positions.png",dpi=155)
    plt.close(fig)


def page(rec):
    g, p = rec["geometry"],rec["phases"]
    m = rec["csc_master"]
    association_text = ", ".join(f"{r['match_type']} / {r['var_index_b'] or 'unreported'}"
                                 for r in rec["same_observation_csc_rows"])
    src = rec["targeting"]
    edge = g["extraction_circle"]["ellipse_intersection_over_scale_range"]
    overlap = ("overlaps the original extraction aperture" if all(edge) else
               "misses the original extraction aperture" if not any(edge) else "touches the aperture scale uncertainty")
    def how(b):
        return "inside" if b else "outside"
    def conf(stat):
        if not stat["n"]:
            return "—"
        a,b = stat["cp95"]
        return f"{stat['k']}/{stat['n']} ({a:.2f}–{b:.2f})"
    lines=[f"# Two-aperture photometry in ObsID {rec['obsid']}", "",
           f"During the registered trial, the equal disjoint cores contained {p['trigger']['trial_core']} photons "
           f"at the extraction center and {p['trigger']['csc_core']} at {m['name']}. Of the "
           f"{p['trigger']['original_aperture']} original trial photons, {p['trigger']['original_photons_in_csc_core']} "
           f"landed in the CSC core and {p['trigger']['original_photons_outside_both_cores']} fell outside both cores.", "",
           f"The CSC center is {g['extraction_circle']['separation_arcsec']:.3f} arcsec from the extraction center "
           f"and {how(g['extraction_circle']['csc_center_in_circle'])} the original aperture. Its separation from the "
           f"reported X-ray localization is {g['localization_circle']['separation_arcsec']:.3f} arcsec. The "
           f"{g['localization_circle']['radius_arcsec']:.3f}-arcsec X-ray error circle "
           f"{'intersects' if g['localization_circle']['csc_ellipse_intersects_circle'] else 'does not intersect'} "
           f"the CSC positional ellipse ({m['err_ellipse_r0']} × {m['err_ellipse_r1']} arcsec). "
           f"The CSC ellipse {overlap} in the local tangent plane. The named observation target was {src['target_name']}. "
           f"The CSC same-ObsID stack match type and "
           f"broad-band variability index are {association_text}. Code u denotes a unique CSC stack association.", "",
           f"The core radius is {g['core_radius_px']:.3f} event SKY pixels. Each core uses the same local "
           "background annulus around the midpoint. Scaling divides nominal sky areas. Edge or dither effects "
           "can change the effective areas.", "",
           "| Phase | CCD live s | Original aperture | Trial / CSC core | Shared annulus | Expected background / core | "
           "Trial share (exact 95%) | Common supported s |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for k in PHASES:
        v=p[k]
        lines.append(f"| {k} | {v['gt_live_s']:.2f} | {v['original_aperture']} | {v['trial_core']} / {v['csc_core']} | "
                     f"{v['shared_annulus']} | {v['nominal_background_per_core']:.2f} | "
                     f"{conf(v['trial_fraction_among_both_cores'])} | {v['common_geometry']['supported_live_s']:.2f} |")
    test=rec["comparison"]["all_recorded"]
    desc=test["fisher_two_sided_descriptive_p"]
    lines.extend(["", f"The two-sided Fisher comparison of the trial with the pooled flanks gives "
                  + (f"p = {desc:.4g}." if desc is not None else "no estimate because a margin is zero."),
                  "The table conditions on photons in the two cores. Both centers were selected after inspection, "
                  "and the flare center was determined using the discovery data. Sparse counts, diffuse emission, "
                  "PSF wings, astrometric offsets and incomplete geometry limit source attribution. The "
                  "geometry column gives times supported by independent detector witnesses around the entire construction. "
                  "The analysis uses no PSF model or exposure map. Its count test has no discovery false-alarm interpretation.", "",
                  "Blue marks the blind extraction core, orange the CSC core, the dashed line the original aperture, "
                  "and dotted lines the background annulus. Points are individual selected event rows.", "",
                  "![Three phases in matched SKY coordinates](positions.png)", "",
                  "[Counts, associations, detector support and input receipts](measurement.json) · "
                  "[Tagged EVT2 rows and both offsets](photons.json.gz) · "
                  "[CSC position and error ellipse definitions](https://cxc.cfa.harvard.edu/csc2.2/columns/positions.html)", ""])
    return "\n".join(lines)


def run(root, retained):
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO and
         os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH and
         os.getenv("GITHUB_WORKFLOW_REF","").split("@",1)[0]==ingest.REPO+"/"+ingest.FLOW, "Use the Chandra workflow")
    out=root/"csc_deblend"/PRIOR/f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    out.mkdir(parents=True,exist_ok=True)
    inv=root/"inventory/inventory.json"
    need(sha(inv)==INVENTORY,"Frozen pointing metadata changed")
    inventory={int(r["obsid"]):r for r in campaign.read(inv)}
    status={"source_run":ingest.RUN,"archived_pointed_run":PRIOR,"implementation_sha256":sha(Path(__file__)),
            "execution_commit":os.environ["GITHUB_SHA"],"targets":{}}
    dump(out/"status.json",status)
    failed=[]
    for group,artifact,oid,h,png,tag,_ in pointed_visits.COHORTS["high-latitude"]:
        area=out/str(oid)
        area.mkdir(exist_ok=True)
        attempts=[]
        state=status["targets"][str(oid)]={"status":"in progress","stage":"validating source rows"}
        dump(out/"status.json",status)
        try:
            product,event,receipt,member=pointed_visits.original(root,group,artifact,oid,h,png,tag)
            entries=campaign.read_obs(root/ingest.HISTORY/"members"/f"{group}.json.gz")
            original_path=f"photometry/{oid}/{product['filename']}/{tag}/event_rows.csv.gz"
            _,amember=pointed_visits.archived_member(root,entries,original_path)
            need(event["diagnostics"]["input_sha256"]==h and amember["sha256"]==event["diagnostics"]["files"]["event_rows.csv.gz"],
                 "Original event-row manifest differs")
            master,same,catalog=csc_rows(root,oid)
            state["stage"]="replaying pinned EVT2"
            dump(out/"status.json",status)
            with TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as temp:
                info=followup.event_bytes(product,oid,Path(temp),retained,attempts)
                path=Path(info["path"])
                data,lo,hi=diagnostics._get_data(path,event)
                with fits.open(path,memmap=False) as hdul:
                    ge=geometry(_sky_wcs(hdul["EVENTS"]),event,master)
                dist=np.hypot(data["x"]-ge["trial_xy"][0],data["y"]-ge["trial_xy"][1])
                r=ge["original_radius_px"]
                check=Path(temp)/"rows-check.csv.gz"
                t=event["peak_start_met_s"]
                diagnostics._events_csv(check,data,dist<=r,(dist>=2*r)&(dist<=4*r),event["start_met_s"],
                                        event["stop_met_s"],t,t+event["peak_window_s"])
                need(sha(check)==amember["sha256"],"Native EVT2 rows differ from archived original extraction")
                phases,tags,comp=measure(data,ge,event,lo,hi)
            photonfile=area/"photons.json.gz"
            campaign.save_obs(photonfile,tags)
            target={k:inventory[oid][k] for k in ("target_name","proposal_number","proposal_title","catalog_ra_deg",
                                                  "catalog_dec_deg","observation_url","proposal_url")}
            tar=(target["catalog_ra_deg"],target["catalog_dec_deg"])
            if None not in tar:
                target["csc_offset_from_pointing_reference_arcsec"]=separation(*tar,float(master["ra"]),float(master["dec"]))
                target["trial_offset_from_pointing_reference_arcsec"]=separation(*tar,*ge["extraction_circle"]["center_icrs_deg"])
            rec={"obsid":oid,"trigger_key":tag,"input":product,"source_ledger_member":member,"source_row_member":amember,
                 "source_artifact":artifact,"source_zip_digest":receipt["zip_digest"],"same_observation_csc_rows":same,
                 "csc_master":master,"csc_receipts":catalog,"targeting":target,"inventory_sha256":INVENTORY,"geometry":ge,
                 "phases":phases,"comparison":comp,"tagged_photon_count":len(tags),"tagged_photons_sha256":sha(photonfile),
                 "source_stage":{"run":ingest.RUN,"code_commit":ingest.CODE,"plan_commit":ingest.PLAN,"group":group},
                 "trigger_ccd_id":event["ccd_id"],"live_fraction":event["live_fraction"],
                 "original_trigger_vet_status":event["vet"]["status"],"primary_control_replicas":24,"heldout_selection_changed":False}
            dump(area/"measurement.json",rec)
            picture(area,ge,tags)
            (area/"README.md").write_text(page(rec))
            state.update(status="complete",stage="recorded",measurement_sha256=sha(area/"measurement.json"),
                         figure_sha256=sha(area/"positions.png"))
        except Exception as exc:
            state.update(status="incomplete",error=f"{type(exc).__name__}: {exc}")
            failed.append(oid)
        finally:
            dump(area/"acquisition.json",attempts)
            dump(out/"status.json",status)
    if failed:
        raise RuntimeError(f"Incomplete CSC deblending measurements for {failed}")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence",required=True,type=Path)
    parser.add_argument("--retained",action="append",default=[],type=Path)
    args=parser.parse_args()
    run(args.evidence,args.retained)
