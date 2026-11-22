"""Imaging-led EW arm in independently observed Ly-alpha emitters.

Every catalogued repeat-LAE source and every positional match to a good DESI
visual redshift is considered. Reserve independent broadband follow-up for
both populations, including a fixed random subset to avoid targeting only
bright pre-existing optical counterparts. All science runs inside Actions.
"""

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import math
import os
from pathlib import Path
import random

import numpy as np
from scipy.spatial import cKDTree

import spectral as spec
from photometry import HSC_R, hsc_scale, imaging_ew, lan_test, num
from vetting import Ancillary, Cache, atom, chord, cube_test, hsc_photometry, tap_batch


def stage_counts(cat, aux, stats):
    ra = np.asarray(cat.info["RA"])
    dec = np.asarray(cat.info["DEC"])
    z = np.asarray(cat.info["z_hetdex"],dtype=float)
    is_lae = np.char.strip(np.char.lower(np.asarray(cat.info["source_type"],dtype="U8")))=="lae"
    repeated = np.zeros(cat.size, dtype=bool)
    pair = cat.sorted_src_ids[:-1] == cat.sorted_src_ids[1:]
    repeated[cat.src_order[:-1][pair]] = True
    repeated[cat.src_order[1:][pair]] = True
    if getattr(cat,"confirmed_repeat_rows",None):
        repeated[list(cat.confirmed_repeat_rows)]=True
    repeated &= is_lae & (z>=1.88) & (z<=3.52)
    stats["repeat_LAE_source_observations"] = int(repeated.sum())
    priors = {int(i):(float(z[i]), "repeat") for i in np.flatnonzero(repeated)
              if len(cat.repeat_rows(i))}
    if "DESI" in aux.tables:
        table = aux.tables["DESI"][1]
        v = np.asarray(table["VI_Z"],dtype=float)
        secure = ((table["VI_QUALITY"]>=3)&(table["DEX_FLAG"]==1)&(table["SEP"]<=1.25) &
                  (v>=1.88)&(v<=3.52))
        stats["DESI_VAC_highz_secure_input"] = int(np.sum(secure))
        valid = np.isfinite(ra)&np.isfinite(dec)&(ra>=0)&(ra<360)&(abs(dec)<=90)
        idx = np.flatnonzero(valid)
        tree = cKDTree(chord(ra[idx],dec[idx]))
        for j in np.flatnonzero(secure):
            if not np.isfinite(table["TARGET_RA"][j]+table["TARGET_DEC"][j]):
                continue
            rows = tree.query_ball_point(chord(table["TARGET_RA"][j],table["TARGET_DEC"][j]),
                                         np.deg2rad(1.25/3600))
            for i in idx[rows]:
                priors[int(i)] = (float(v[j]),"desi")
    stats["eligible_source_observations"] = len(priors)
    stats["DESI_counterpart_observations"] = sum(t=="desi" for _,t in priors.values())
    return priors


def compact(row):
    return {k:row.get(k) for k in ("row_index","source_id","detectid","shotid","ra","dec",
            "wave","z_lya","flux","flux_err","snr","continuum","continuum_err","ew0",
            "ew0_lower","rescue_population","imaging_rank_A")}


def measure_spectra(cat, priors, aux, stats, ledger, shuffled=False):
    eligible = {}
    rows = np.array(sorted(priors), dtype=int)
    for first in range(0,len(rows),64):
        ix = rows[first:first+64]
        f = np.array(cat.flux[ix],dtype=np.float32,order="C")
        e = np.array(cat.error[ix],dtype=np.float32,order="C")
        bad = (f==-999) | ~np.isfinite(f) | ~np.isfinite(e) | (e<=0)
        e[bad]=np.nan
        if shuffled:
            if len(ix)<2:
                stats["control_singleton_not_shuffled"]+=len(ix)
                continue
            f,e = spec.shuffle_fibres(f,e,seed=572025+first)
        stats["spectra_measured"]+=len(ix)
        all_peaks = defaultdict(list)
        for local,peak in spec.search_batch(cat.wave,f,e,density_unit=1e-17):
            i = int(ix[local])
            z0 = priors[i][0]
            if not (peak["single_pixel"] or peak["sky_near"] or peak.get("badcal",False)
                    or abs(peak["z_lya"]-z0)>.005):
                all_peaks[i].append(peak)
        for i in map(int,ix):
            m = cat.metadata(i)
            st = priors[i][1]
            base = dict(row_index=i,source_id=m["source_id"],detectid=m["detectid"],
                        shotid=m["shotid"],ra=m["ra"],dec=m["dec"],rescue_population=st)
            def record(reason, extra=None):
                stats[reason] += 1
                ledger.write(json.dumps({**base,"status":reason,**(extra or {})},default=atom)+"\n")
            if not all_peaks[i]:
                record("no_clean_peak_at_redshift")
                continue
            p = max(all_peaks[i],key=lambda p:p["snr"])
            if p["ew_extreme"]:
                record("already_in_primary_spectral_EW_arm",p)
                continue
            c,ec = p["continuum"],p["continuum_err"]
            if c is not None and ec is not None and c>3*ec and (
                (p["flux"]+3*p["flux_err"])/((c-3*ec)*(1+p["z_lya"]))<=240):
                record("spectral_continuum_rules_out_extreme_EW",p)
                continue
            if m.get("source_type")=="agn" or m.get("agn_flag") in (0,1):
                record("catalogued_AGN_or_broad_source",p)
                continue
            _,lan = aux.lookup("LAN",m["ra"],m["dec"],3.)
            if lan and abs((num(lan.get("z_hetdex"),-999)-p["z_lya"]))>.015:
                lan=None
            if lan and num(lan.get("combined_eqw_rest_lya"),0)>240:
                record("already_published_LAN_extreme_EW",p)
                continue
            # An HSC aperture with observed flux & uncertainty provides a priority
            # proxy. The published r=30 nondetection flag carries no optical depth.
            # A different line/continuum aperture cannot certify a discovery.
            rank_hint = 0.
            lan_check = lan_test(dict(m,**p),lan)
            if lan_check:
                rank_hint = num((lan_check.get("catalog_aperture_test") or {}).get("ew0_lower_A"),0.)
            rec = {**m,**base,**p,"channel":"ew","stream":"imaging_rescue",
                   "catalog":m,"rescue_prior_z":priors[i][0],"catalog_ew_priority_only_A":rank_hint}
            if lan:
                rec["rescue_LAN"] = {k:lan.get(k) for k in ("ra","dec","detectid","shotid","z_hetdex",
                    "combined_eqw_rest_lya","flag_resolved","sep_arcsec","flux_lya","flux_lya_err",
                    "area_r_iso_circ","hsc-r_mag","hsc-r_mag_err")}
            key = int(m["source_id"])
            old = eligible.get(key)
            if old is None or (st=="desi",p["snr"])>(old["rescue_population"]=="desi",old["snr"]):
                eligible[key]=rec
            record("science_peak_for_imaging",p)
    return list(eligible.values())


def select_scout(rows, population, limit, seed):
    """Priority thirds: available LAN photometry, robust line flux, random remainder."""
    rows = [r for r in rows if r["rescue_population"]==population]
    if len(rows)<=limit:
        return rows
    out, done = [], set()
    orders = [sorted(rows,key=lambda r:-r["catalog_ew_priority_only_A"]),
              sorted(rows,key=lambda r:-(r["flux"]-3*r["flux_err"])),
              random.Random(seed).sample(rows,len(rows))]
    for j,order in enumerate(orders,1):
        for row in order:
            key=row["source_id"]
            if key in done or (j==1 and row["catalog_ew_priority_only_A"]<=0):
                continue
            out.append(row)
            done.add(key)
            if len(out)>=limit*j//3:
                break
    return out


def scout_photometry(rec, cache):
    audit = rec["rescue_audit"]
    try:
        audit["HSC"] = hsc_photometry(rec)
    except Exception as exc:
        audit["HSC"] = dict(status="unavailable",error=str(exc))
    if not (hsc_scale(audit["HSC"]) and num(audit["HSC"].get("flux_err"),0)>0):
        try:
            audit["Legacy_forced"] = hsc_photometry(rec,layer="ls-dr9")
        except Exception as exc:
            audit["Legacy_forced"] = dict(status="unavailable",error=str(exc))
    try:
        audit["cube"] = cube_test(rec,cache.fetch)
    except Exception as exc:
        audit["cube"] = dict(status="unavailable",error=str(exc))
    audit["imaging_EW"] = imaging_ew(rec,audit)
    return rec


def rescue(cat,out):
    out = Path(out)
    cache = Cache(Path(os.environ.get("RUNNER_TEMP","/tmp"))/"hetdex")
    aux = Ancillary(cache)
    stats = Counter()
    try:
        priors = stage_counts(cat,aux,stats)
        observed=[]
        cohort=dict(stats)
        for shuffled in (False,True):
            stats=Counter(cohort)
            prefix="ew-rescue-control" if shuffled else "ew-rescue"
            with gzip.open(out/(prefix+"-spectra.jsonl.gz"),"wt") as ledger:
                leads = measure_spectra(cat,priors,aux,stats,ledger,shuffled)
            stats["distinct_possible_after_deduplication"] = len(leads)
            sample = select_scout(leads,"desi",120,746321)+select_scout(leads,"repeat",120,746322)
            stats["bounded_imaging_positions"] = len(sample)
            for offset in range(0,len(sample),12):
                block = sample[offset:offset+12]
                for rec,ls in zip(block,tap_batch("Legacy",block)):
                    rec["rescue_audit"]={"Legacy":ls}
                    if rec.get("rescue_LAN"):
                        rec["rescue_audit"]["LAN"]=rec["rescue_LAN"]
                    ew=imaging_ew(rec,rec["rescue_audit"])
                    rec["imaging_rank_A"]=num(ew.get("ew0_lower_A"),0.)
            # Among objects without adequate direct Legacy photometry, seek real
            # sky-aperture measurements. No optical catalogue match is not a limit.
            extra=[]
            for j,pop in enumerate(("desi","repeat")):
                untested = [r for r in sample if r["rescue_population"]==pop and r["imaging_rank_A"]==0.]
                extra.extend(select_scout(untested,pop,12,847328+j))
            with ThreadPoolExecutor(max_workers=3) as pool:
                for _ in pool.map(lambda r:scout_photometry(r,cache),extra):
                    pass
            stats["paired_aperture_attempts"]=len(extra)
            with gzip.open(out/(prefix+"-imaging.jsonl.gz"),"wt") as ledger:
                for rec in sample:
                    audit=rec["rescue_audit"]
                    ew=audit.get("imaging_EW") or imaging_ew(rec,audit)
                    rec["imaging_rank_A"]=num(ew.get("ew0_lower_A"),0.)
                    stats["eligible_imaging_EWs"]+=int(ew.get("ew0_lower_A") is not None)
                    stats["imaging_EW_above_240"]+=int(rec["imaging_rank_A"]>240)
                    ledger.write(json.dumps({"source":compact(rec),"imaging":ew,"evidence":audit},
                                             default=atom,allow_nan=False)+"\n")
            summary={"counts":dict(stats),"source_catalog_errors":aux.errors,
                     "inputs":cache.inputs,
                     "scope":"All repeat-LAE observations with same-source IDs or isolated positional repeats and good matched DESI "
                             "VI z=1.88-3.52 considered in this release. At most 120 unique objects per "
                             "population receive direct Tractor checks, prioritized in thirds by LAN "
                             "optical priority, lower line flux and fixed-seed sampling. Up to 12 with "
                             "missing/zero imaging bound per population get paired-aperture attempts "
                             "before final external vetting. No source is promoted by a missing flux."}
            (out/(prefix+"-summary.json")).write_text(json.dumps(summary,indent=2,default=atom)+"\n")
            if not shuffled:
                observed=sample
        return observed
    finally:
        aux.close()
