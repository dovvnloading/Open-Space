"""Recheck two previously selected Chandra flares from their original Actions evidence."""
import argparse
from pathlib import Path
import csv
import json
import hashlib
import os
import shutil
from tempfile import TemporaryDirectory

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact

from archive import dump, replay
from campaign import read, read_obs
from diagnostics import _get_data
from search import _sky_wcs, PIXEL_ARCSEC

PRIMARY = "44985879443"
REVISION = "4663ab98666b5d9b5b116478f61d6cae5d49e104"
TARGETS = ((31,14744596670,23324,"507a5aaf26b492ef2b5120242d9d33e281803667c8c1edc125fa905e6a09bbd0"),
           (19,14743042773,32380,"3c5f46567187114d7deb741d2a709bbc4e3566801a90a9cdb3f23a7c66a944cf"))


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def seed(staging, out):
    leads = []
    for group, artifact, oid, product_hash in TARGETS:
        origin = staging / f"chandra-evidence-{group}"
        env = read(origin / f"environment-{group}.json")
        if str(env["run_id"]) != PRIMARY or env["git_commit"] != REVISION:
            raise ValueError("Unexpected producer")
        ledger_path = origin / "observations" / f"{oid}.json.gz"
        ledger = read_obs(ledger_path)
        product, = [p for p in ledger["products"] if p["input"]["sha256"] == product_hash]
        rec = product["input"]
        t = product["search"]["triggers"][0]
        if ledger["obsid"] != oid or t["trigger_key"] != "real/0" or not t["instrument_pass"]:
            raise ValueError("Original trigger identity disagrees")
        dest = out / str(oid)
        dest.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ledger_path,dest/ledger_path.name)
        photo = origin/"photometry"/str(oid)/rec["filename"]/"real/0"
        for name,sha in t["diagnostics"]["files"].items():
            if Path(name).name != name or digest(photo/name) != sha:
                raise ValueError("Source diagnostic hash differs")
        shutil.copytree(photo,dest/"photometry",dirs_exist_ok=True)
        original_cache = origin/t["catalog_cache"]
        extracted = dest/"original_catalogs"
        extracted.mkdir(exist_ok=True)
        for query in t["vet"]["queries"].values():
            if "url" not in query:
                continue
            key=hashlib.sha256(query["url"].encode()).hexdigest()[:24]
            for memo in original_cache.glob(f"*-{key}.json"):
                contents=json.loads(memo.read_text())
                if contents.get("url")!=query["url"]:
                    raise ValueError("Mismatched archived catalog request")
                shutil.copyfile(memo,extracted/memo.name)
                for record in [contents,*contents.get("attempts",[])]:
                    name=record.get("body_file")
                    if not name:
                        continue
                    if Path(name).name!=name or digest(original_cache/name)!=record["sha256"]:
                        raise ValueError("Mismatched archived catalog response")
                    shutil.copyfile(original_cache/name,extracted/name)
        dump(dest/"trigger.json", t)
        dump(dest/"input.json", [rec])
        leads.append({k:t[k] for k in ("obsid","trigger_key","ra_deg","dec_deg","event_mjd",
                                        "position_error_arcsec","trial_duration_s","duration_s")})
        leads[-1].update(input_sha256=product_hash,input_name=rec["filename"],input_url=rec["url"],
                         status=t["vet"]["status"],original_artifact_id=artifact,original_shard=group,
                         original_ledger_sha256=digest(ledger_path),
                         peak_flux_erg_cm2_s=t["flux"]["peak_flux_erg_cm2_s"],
                         counting_sigma=t["flux"]["flux_counting_sigma_erg_cm2_s"],
                         gaia_source_ids=[str(m["name"]) for m in t["vet"]["matches"] if m["catalog"]=="gaia"])
    dump(out/"leads.json",leads)
    dump(out/"origin.json",{"primary_run":PRIMARY,"primary_workflow_commit":REVISION,
        "targeted_followup_run":os.environ["GITHUB_RUN_ID"],"selection":"Observed sky-localized flares. "
        "ObsID 23324 is offset from the pointed-at star. ObsID 32380 targeted a ROSAT position. "
        "This follow-up does not modify either trigger or the primary survey false-alarm calibration."})


def counts(source_mask, periods):
    return {name:int(np.count_nonzero(source_mask & mask)) for name,mask in periods.items()}


def comparison(data, t, evt):
    rows, lo, hi = _get_data(evt,t)
    start,stop = t["start_met_s"],t["stop_met_s"]
    periods={"before":rows["t"] < start,"flare":(rows["t"]>=start)&(rows["t"]<stop),"after":rows["t"]>=stop}
    xy = np.column_stack((rows["x"],rows["y"]))
    center = np.array((t["aperture_x"],t["aperture_y"]))
    radius=t["aperture_radius_px"]
    dist=np.linalg.norm(xy-center,axis=1)
    aperture=dist<=radius
    src=counts(aperture,periods)
    if [src[k] for k in ("before","flare","after")] != [t["before_counts"],t["on_counts"],t["after_counts"]]:
        raise ValueError("Original source aperture counts disagree")
    # Select the reference exclusively from the pre- and post-window image.
    # This is a comparison for PSF spillover and correlated variability.
    outside=~periods["flare"]
    bright=xy[outside & (dist>2.2*radius) & (dist<7*radius)]
    measures={"source_counts":src,"source_half_radius_counts":counts(dist<=radius/2,periods),
        "gtis_met_s":[lo,hi],"trial_met_s":[start,stop],"radius_sky_px":radius,
        "source_annulus_counts":counts((dist>=2*radius)&(dist<=4*radius),periods),
        "spectral_slices":[],"comparison":None,"input_sha256":digest(evt)}
    for i,(a,b) in enumerate(zip(np.linspace(start,stop,4)[:-1],np.linspace(start,stop,4)[1:])):
        photons=rows["e"][aperture & (rows["t"]>=a) & (rows["t"]<b)]/1000
        measures["spectral_slices"].append({"third":i+1,"n":len(photons),"n_ge_2kev":int((photons>=2).sum()),
                                           "median_keV":float(np.median(photons)) if len(photons) else None})
    if len(bright)>=5:
        tree=cKDTree(bright)
        count=tree.query_ball_point(bright,max(3,radius*.6),return_length=True)
        chosen=int(np.argmax(count))
        ref=np.mean(bright[tree.query_ball_point(bright[chosen],max(3,radius*.6))],axis=0)
        ref_dist=np.linalg.norm(xy-ref,axis=1)
        if np.linalg.norm(ref-center)>2*radius and count[chosen]>=5:
            ref_ap=ref_dist<=radius
            neighbor=counts(ref_ap,periods)
            with fits.open(evt,memmap=False) as hdul:
                pos=_sky_wcs(hdul["EVENTS"]).all_pix2world([ref],1)[0]
            table=[[src["flare"],neighbor["flare"]],[src["before"]+src["after"],neighbor["before"]+neighbor["after"]]]
            measures["comparison"]={"center_sky_px":ref.tolist(),"ra_deg":float(pos[0]),"dec_deg":float(pos[1]),
                "separation_arcsec":float(np.linalg.norm(ref-center)*PIXEL_ARCSEC),"counts":neighbor,
                "outside_window_cluster_photons":int(count[chosen]),"source_neighbor_table":table,
                "share_of_photons_fisher_two_sided_p":float(fisher_exact(table).pvalue),
                "source_annulus_excluding_neighbor_counts":counts((dist>=2*radius)&(dist<=4*radius)&~ref_ap,periods)}
    else:
        ref_ap=None
    if measures["comparison"] is None:
        ref_ap=np.zeros(len(rows),dtype=bool)
    edges=np.arange(lo,hi,120.).tolist()+[hi]
    columns=[np.histogram(rows["t"][mask],edges)[0] for mask in (aperture,ref_ap,(dist>=2*radius)&(dist<=4*radius))]
    with (data/"lightcurves_120s.csv").open("w",newline="") as fh:
        w=csv.writer(fh);w.writerow(["start_met_s","stop_met_s","source_counts","comparison_counts","annulus_counts"])
        w.writerows(zip(edges[:-1],edges[1:],*columns))
    measures["interpretation"]=("The selected neighboring aperture is a diagnostic of contemporaneous flux changes. "
        "Its off-window selection and spatially varying backgrounds prevent use of this test as a discovery probability. "
        "Empty comparison means fewer than five sufficiently clustered off-window photons or overlapping apertures. "
        "Colors and finite-aperture counts do not establish a distance or distinguish thermal from nonthermal emission.")
    dump(data/"photon_check.json",measures)


def photons(out):
    for lead in read(out/"leads.json"):
        dest=out/str(lead["obsid"])
        with TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as tmp:
            directory=Path(tmp)
            replay(dest/"input.json",directory)
            comparison(dest,read(dest/"trigger.json"),directory/lead["input_name"])


if __name__=="__main__":
    if os.getenv("GITHUB_ACTIONS")!="true" or os.getenv("GITHUB_REF")!="refs/heads/research/chandra-acis-fast-transients":
        raise RuntimeError("Use the Chandra workflow on its research branch")
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="cmd",required=True)
    s=sub.add_parser("seed");s.add_argument("staging",type=Path);s.add_argument("output",type=Path)
    s=sub.add_parser("photons");s.add_argument("output",type=Path)
    args=p.parse_args()
    if args.cmd=="seed":
        seed(args.staging,args.output)
    else:
        photons(args.output)
