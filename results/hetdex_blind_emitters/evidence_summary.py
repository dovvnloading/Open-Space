"""Auditable release-level tallies from completed HETDEX search ledgers (stdlib only)."""
from collections import Counter, defaultdict, deque
import gzip
import hashlib
import json
import math
from pathlib import Path


IONS = {"civ": ("civ1", "civ2"), "oiii": ("oiii_uv1", "oiii_uv2"), "ciii": ("ciii1", "ciii2")}
ARCHIVES = ("NED", "SIMBAD", "Legacy", "HSC", "Legacy_forced", "SPARCL", "cube", "DESI", "DESI_lines", "LAN", "HETDEX_HSC")


def num(value, default=None):
    try:
        n = float(value)
        return n if math.isfinite(n) else default
    except (TypeError, ValueError, OverflowError):
        return default


def read(path):
    with path.open() as infile:
        return json.load(infile)


def ledger(path):
    with gzip.open(path, "rt") as infile:
        for number, line in enumerate(infile, 1):
            yield number, json.loads(line)


def stable_key(p):
    return p["row_index"],p["channel"],p.get("stream") or "",float(p["wave"]).hex()


def measurement_id(p):
    """Immutable observed quantities; source and partner IDs distinguish blended lines."""
    fields = ("shotid", "source_id", "detectid", "ra", "dec", "flux", "flux_err", "snr",
              "knot_detectid", "secondary_detectid", "secondary_ra", "secondary_dec",
              "secondary_wave", "repeat_rows", "imaging_rank_A")
    partner = p.get("pair_lya") or {}
    values = [p.get(k) for k in fields] + [partner.get(k) for k in ("wave","flux","flux_err","snr")]
    return stable_key(p), json.dumps(values,separators=(",",":"),allow_nan=False)


def frozen_digest(peak, vet):
    raw = json.dumps([peak,vet],sort_keys=True,separators=(",",":"),allow_nan=False)
    return hashlib.sha256(raw.encode()).digest()


def checked_line(row):
    n = row["first_pass_line"]
    if type(n) is not int or n < 1:
        raise ValueError("Invalid first-pass line in supplementary ledger")
    return n


def legacy_focus(focus, revisions):
    """Recover exact original snapshots saved before the line-number field existed."""
    by_measurement = defaultdict(list)
    for _,r in ledger(focus/"followup.jsonl.gz"):
        if r.get("first_pass_line") is not None:
            continue
        ident = (measurement_id(r["peak"]),r["group"],r["first"]["disposition"],r["first_pass_external"])
        # Historical focus.adjudicate processed this stable sort of the immutable
        # input list. Replay within each equal measurement group only. Distinct
        # Ly-alpha anchors and distinct source/visit geometries never share a queue.
        order = (r["first_pass_external"],not r["independent_redshift_agrees"],
                 not bool(r["same_visit_two_line_repeats"]),r["group"])
        digest = frozen_digest(r["peak"],r["first"])
        by_measurement[ident].append((order,digest))
    ordered = {k:deque(h for _,h in sorted(v,key=lambda pair:pair[0]))
               for k,v in by_measurement.items()}
    by_snapshot = defaultdict(lambda:defaultdict(deque))
    for line,row,path in revisions:
        p = row["peak"]
        ident = (measurement_id(p),row["group"],row["first_pass_disposition"],row["first_pass_external"])
        if not ordered.get(ident):
            raise ValueError(f"Legacy focus row {line} has no matching frozen measurement and companion")
        digest = ordered[ident].popleft()
        by_snapshot[stable_key(p)][digest].append((line,row,path))
    if any(ordered.values()):
        raise ValueError("The frozen focus selection and completed legacy decisions differ in multiplicity")
    return by_snapshot


def adjudications(d, *, include_pair_replay=True):
    """Yield exactly one effective decision per original ledger row, in original order."""
    d=Path(d)
    focus=d/"focus"
    spatial=d/"repeat-blobs"
    register=d/"vetted.jsonl.gz"
    exact={}
    expected_hashes=set()
    old=[]
    if (focus/"vet-summary.json").exists():
        expected_hashes.add(read(focus/"summary.json")["input_sha256"])
        counts=Counter()
        for line,row in ledger(focus/"vetted.jsonl.gz"):
            if tuple(row["key"]) != stable_key(row["peak"]):
                raise ValueError(f"Focus identity disagrees with measured peak on supplemental line {line}")
            counts[row["vet"]["disposition"]]+=1
            item=(line,row,"focus/vetted.jsonl.gz")
            if row.get("first_pass_line") is None:
                old.append(item)
            else:
                n=checked_line(row)
                if n in exact:
                    raise ValueError(f"Multiple focus decisions claim first-pass line {n}")
                exact[n]=item
        summary=read(focus/"vet-summary.json")
        if sum(counts.values())!=summary["rows"] or counts!=Counter(summary["counts"]):
            raise ValueError("Focus vetting counts disagree with completed ledger")
    if (spatial/"summary.json").exists():
        expected_hashes.add(read(spatial/"summary.json")["input_sha256"])
        for line,row in ledger(spatial/"vetted.jsonl.gz"):
            n=checked_line(row)
            if n in exact or row["peak"]["channel"]!="blob":
                raise ValueError(f"Invalid spatial revision for first-pass line {n}")
            exact[n]=(line,row,"repeat-blobs/vetted.jsonl.gz")
    replay={}
    saved=focus/"saved-pair-decisions.jsonl.gz"
    if include_pair_replay and saved.exists():
        manifest=focus/"saved-pair-summary.json"
        expected_hashes.add(read(manifest)["input_sha256"])
        for line,row in ledger(saved):
            n=checked_line(row)
            if n in replay or tuple(row["key"])!=stable_key(row["peak"]):
                raise ValueError(f"Invalid saved pair decision for first-pass line {n}")
            replay[n]=(line,row,"focus/saved-pair-decisions.jsonl.gz")
    if expected_hashes:
        with register.open("rb") as source:
            if expected_hashes!={hashlib.file_digest(source,"sha256").hexdigest()}:
                raise ValueError("A supplemental product belongs to a different frozen first-pass ledger")
    legacy=legacy_focus(focus,old) if old else {}
    for n,entry in ledger(register):
        peak=entry["peak"]
        key=stable_key(peak)
        original=entry["vet"]["disposition"]
        item=exact.pop(n,None)
        if item is None and key in legacy:
            digest=frozen_digest(peak,entry["vet"])
            choices=legacy[key].get(digest)
            if choices:
                item=choices.popleft()
                if not choices:
                    del legacy[key][digest]
                if not legacy[key]:
                    del legacy[key]
        if item is not None:
            line,row,path=item
            if (measurement_id(peak)!=measurement_id(row["peak"]) or
                row.get("first_pass_disposition",original)!=original):
                raise ValueError(f"Supplemental measurement differs from first-pass line {n}")
        else:
            line,row,path=n,entry,"vetted.jsonl.gz"
        if n in replay:
            newer=replay.pop(n)
            prev=newer[1]
            if ((prev["previous_ledger"],prev["previous_ledger_line"],prev["previous_disposition"])!=
                (path,line,row["vet"]["disposition"]) or measurement_id(row["peak"])!=measurement_id(prev["peak"])):
                raise ValueError(f"Replayed pair does not match the effective decision on first-pass line {n}")
            line,row,path=newer
        yield line,row,path,original
    if exact or legacy or replay:
        raise ValueError("A supplementary decision lacks an identical original ledger record")


def keep(obj, keys):
    return {k: obj[k] for k in keys.split() if k in obj}


def ratio(n, d):
    return {"numerator": n, "denominator": d, "fraction": n / d if d else None}


def archives(v, tally):
    for key in ARCHIVES:
        obj = v.get(key)
        if not isinstance(obj, dict):
            tally[key]["not_recorded"] += 1
            continue
        st = obj.get("status", "")
        if key in ("LAN", "DESI", "HETDEX_HSC"):
            outcome = "matched"
        elif key == "DESI_lines" and obj.get("lines"):
            outcome = "spectral_diagnostics_recorded"
        elif st in ("queried", "matched", "no_counterpart"):
            outcome = "matched" if obj.get("matches") else "empty_cone"
        elif st == "unavailable" or obj.get("error"):
            outcome = "error"
        elif st == "measured":
            outcome = "measured"
        else:
            outcome = st or "status_unspecified"
        tally[key][outcome] += 1
        if key == "SPARCL":
            for m in obj.get("matches", []):
                tally["SPARCL_individual_spectra"][m.get("status", "status_unspecified")] += 1


def repeat_sn(v, key, primary_shot):
    """One maximum per distinct recorded other shot; no entries for masked fits."""
    by_shot = {}
    for m in v.get(key) or []:
        shot = m.get("shotid")
        sn = num(m.get("snr"))
        if shot is not None and str(shot) != str(primary_shot) and sn is not None:
            by_shot[shot] = max(sn, by_shot.get(shot, -math.inf))
    return by_shot


def spatial(p, v, counts):
    central = p.get("morphology") == "central"
    arm = "single_center" if central else "two_knots"
    k = counts[arm]
    k["vetted_rows"] += 1
    c = v.get("cube") or {}
    metric = "annular_psf_excess" if central else "bridge"
    passed = "annular_above_5sigma" if central else "bridge_above_5sigma"
    m = c.get(metric) or {}
    measurable = num(m.get("sn")) is not None
    k["cube_attempted"] += int(bool(c))
    k["cube_opened_and_line_window_tested"] += int(c.get("status") == "measured")
    k["spatial_statistic_measurable"] += int(measurable)
    k["central_5sigma"] += int(bool(c.get("central_above_5sigma")))
    excess = bool(c.get(passed))
    clean = excess and c.get("central_above_5sigma") and not c.get("background_suspect")
    k["psf_excess_5sigma"] += int(excess)
    k["central_and_excess_and_clean_background"] += int(bool(clean))
    reps = [x["bridge_cube"] for x in v.get("repeats", []) if x.get("bridge_cube")
            and x.get("shotid") is not None and str(x["shotid"]) != str(p.get("shotid"))]
    second = any(num((r.get(metric) or {}).get("sn")) is not None for r in reps)
    k["different_shot_cube_attempted"] += int(bool(reps))
    k["different_shot_excess_measurable"] += int(second)
    k["first_excess_and_measurable_repeat"] += int(excess and second)
    k["both_shots_excess_5sigma"] += int(excess and any(r.get(passed) for r in reps))
    k["both_shots_central_and_excess_and_clean_background"] += int(bool(clean) and any(
        r.get(passed) and r.get("central_above_5sigma") and not r.get("background_suspect") for r in reps))


def he_ew(p):
    if p.get("heii_rest_ew_A") is not None or p.get("heii_rest_ew_lower_A") is not None:
        return p.get("heii_rest_ew_A"), p.get("heii_rest_ew_lower_A")
    # The single-peak EW and lower bound are tabulated under the Ly-alpha interpretation.
    h, l = num(p.get("z_heii")), num(p.get("z_lya"))
    if h is None or l is None or h <= -1:
        return None, None
    scale = (1 + l) / (1 + h)
    return tuple(None if num(p.get(k)) is None else scale * num(p[k]) for k in ("ew0", "ew0_lower"))


def lya_partner(p):
    # Prefer independently localized blind peaks when provided. Otherwise the stored
    # velocity scan spans -300 to +600 km/s, with fixed-centroid correlated trials.
    for name in ("pair_lya", "lya_velocity_best"):
        obj = p.get(name) or {}
        if num(obj.get("snr")) is not None and num(obj.get("flux"), 0) > 0:
            return name, obj
    return None, {}


def example(p, v, d, line):
    ch = p["channel"]
    ew, elo = he_ew(p) if ch == "heii" else (
        p.get("spectral_EW0_A", p.get("ew0")), p.get("spectral_EW0_lower_A", p.get("ew0_lower")))
    return {
        "ledger": f"evidence/{d.name}/vetted.jsonl.gz", "jsonl_line": line,
        **keep(p, "release row_index source_id detectid shotid ra dec channel stream knot_detectid secondary_detectid"),
        "classification": v["disposition"], "air_wavelength_A": p.get("wave"),
        "flux_cgs": p.get("flux"), "flux_err_cgs": p.get("flux_err"), "snr": p.get("snr"),
        "redshift_hypothesis": p.get("z_heii") if ch == "heii" else p.get("z_lya") if ch != "unknown" else None,
        "rest_EW_A": ew if ch != "unknown" else None, "rest_EW_lower_A": elo if ch != "unknown" else None,
        "imaging_EW_lower_A": (v.get("imaging_EW") or {}).get("ew0_lower_A") if ch == "ew" else None,
        "passed": list(v.get("passed") or []), "open_reasons": list(v.get("reasons") or []),
        "ranked_external_checks_recorded": "SPARCL" in v,
    }


def chemistry(p, v, d, line, stages):
    stages["all_heii_rows_reviewed"] += 1
    if v.get("disposition") != "unresolved":
        return
    stages["unresolved"] += 1
    z, wave, sigma = num(p.get("z_heii")), num(p.get("wave")), num(p.get("sigma_A"))
    flux = num(p.get("flux"), 0)
    if (z is None or not 1.9 <= z <= 2.4 or not wave or not sigma or sigma <= 0 or
        flux <= 0 or num(p.get("snr"), 0) < 4.5 or
        any(p.get(k) for k in ("single_pixel", "sky_near", "badcal"))):
        return
    width = 2.355 * sigma / wave * 299792.458
    if width >= 900:
        return
    stages["clean_narrow_heii"] += 1
    origin, ly = lya_partner(p)
    if num(ly.get("snr"), 0) < 5 or flux / num(ly.get("flux"), math.inf) < 0.10:
        return
    stages["paired_lya_and_flux_ratio"] += 1
    metals = p.get("heii_metal_components") or p.get("pair_metals") or p.get("heii_family") or {}
    doublets = {}
    for ion, keys in IONS.items():
        components = [metals.get(k) or {} for k in keys]
        if any(num(m.get("snr"), -math.inf) >= 3 for m in components):
            return
        measured = all(m.get("status") == "measured" and num(m.get("limit3"), 0) > 0 for m in components)
        doublets[ion] = {"component_status": [m.get("status", "not_recorded") for m in components],
                         "hetdex_ratio_upper": sum(m["limit3"] for m in components)/flux if measured else None}
    limits = dict(v.get("three_ion_limits_relative_to_heii") or {})
    for ion, m in doublets.items():
        u = m["hetdex_ratio_upper"]
        if u is not None and u < num((limits.get(ion) or {}).get("ratio_upper"), math.inf):
            limits[ion] = {"ratio_upper": u, "source": "HETDEX"}
    if not any(doublets[k]["hetdex_ratio_upper"] is not None for k in ("civ", "oiii")):
        return
    stages["no_component_ge3_and_at_least_one_civ_or_oiii_doublet_measured"] += 1
    values = [num((limits.get(k) or {}).get("ratio_upper"), math.inf) for k in IONS]
    ntight = sum(u <= 0.5 for u in values)
    ncoverage = sum(math.isfinite(u) for u in values)
    stages[f"ions_with_ratio_upper_le_0p5_{ntight}"] += 1
    rec = example(p, v, d, line)
    rec.update(partner_measurement=origin, lya=keep(ly,"wave air_wave flux flux_err snr offset_km_s"),
               heii_to_lya=flux/ly["flux"], heii_fwhm_kms=width, metal_doublets=doublets,
               best_ratio_limits=limits, ions_with_ratio_upper_le_0p5=ntight,
               low_z_companion_patterns=v.get("low_z_families") or [],
               single_line_catalog_conflict=bool(v.get("catalog_conflict")))
    independently_seen = (max(repeat_sn(v,"repeats",p.get("shotid")).values(), default=0) >= 4 or
                          "independent instrument confirms line" in rec["passed"])
    rec["independently_reobserved_primary"] = independently_seen
    competing = bool(rec["low_z_companion_patterns"] or rec["single_line_catalog_conflict"])
    rank = (not competing, ntight, independently_seen, ncoverage,
            -sum(x for x in values if math.isfinite(x)), min(num(ly["snr"]),num(p["snr"])), -line)
    return rank, rec


def illustrative(p, v, d, line):
    if v.get("disposition") != "unresolved" or "SPARCL" not in v:
        return
    key = p["channel"]
    cube = v.get("cube") or {}
    repeat = max(repeat_sn(v,"repeats",p.get("shotid")).values(), default=0)
    independent = repeat >= 4 or "independent instrument confirms line" in (v.get("passed") or [])
    probe = {}
    if key == "blob":
        central = p.get("morphology") == "central"
        key = "single_center" if central else "two_knots"
        s = num((cube.get("annular_psf_excess" if central else "bridge") or {}).get("sn"), 0)
        if s < 5 or not cube.get("central_above_5sigma") or cube.get("background_suspect"):
            return
        probe.update(psf_excess_empirical_snr=s, central_empirical_snr=cube["central_aperture"]["sn"],
                     cube=cube.get("cube"),seeing_arcsec=cube.get("seeing_arcsec"))
        rank = (independent, s)
    elif key == "ew":
        spectral = num(p.get("ew0_lower"),0)
        photo = num((v.get("imaging_EW") or {}).get("ew0_lower_A"),0)
        if max(spectral, photo) <= 240 or not (independent or cube.get("central_above_5sigma")):
            return
        rank = (min(spectral,photo)>240, independent, min(spectral,photo), max(spectral,photo))
        probe.update(spectral_lower_A=spectral or None,imaging_lower_A=photo or None)
    elif key == "unknown":
        audit = v.get("unassigned_line_audit") or {}
        if audit.get("status") != "unmatched at independent redshift" or not independent:
            return
        probe.update(host_transition_test=audit["status"])
        rank = (repeat, num(p.get("snr"),0))
    else:
        return
    rec = example(p,v,d,line)
    rec.update(observation=probe,maximum_independent_shot_snr=repeat if repeat else None)
    return key, rank+(-line,), rec


def update_best(pool, key, rank, rec):
    if key not in pool or rank > pool[key][0]:
        pool[key] = (rank,rec)


def summarize_evidence(d):
    d = Path(d)
    archive = defaultdict(Counter)
    arms = {k: Counter() for k in ("single_center", "two_knots")}
    rep = defaultdict(Counter)
    reasons = defaultdict(Counter)
    chem_stages, near, chemicals = Counter(), {}, {}
    for n, entry, path, initial in adjudications(d):
        p, v = entry["peak"], entry["vet"]
        ch = p["channel"]
        archives(v,archive)
        if ch == "blob":
            spatial(p,v,arms)
        r = rep[ch]
        r["all_vetted_rows"] += 1
        visits = repeat_sn(v,"repeats",p.get("shotid"))
        other = repeat_sn(v,"secondary_repeats",p.get("shotid"))
        r["repeat_lookup_recorded"] += int("repeats" in v)
        r["repeat_measurement_available"] += int(bool(visits))
        r["primary_repeats_ge4"] += int(any(sn >= 4 for sn in visits.values()))
        r["both_lines_ge4_in_same_independent_shot"] += int(any(sn >= 4 and other.get(s,0) >= 4 for s,sn in visits.items()))
        if v.get("disposition") == "unresolved":
            reasons[ch].update(set(v.get("reasons") or []))
        if ch == "heii":
            pick = chemistry(p,v,d,n,chem_stages)
            if pick:
                rank, rec = pick
                rec["ledger"]=f"evidence/{d.name}/{path}"
                key = p.get("source_id")
                if key is None:
                    key = (p.get("ra"), p.get("dec"))
                update_best(chemicals,key,rank,rec)
        else:
            pick = illustrative(p,v,d,n)
            if pick:
                key,rank,rec=pick
                rec["ledger"]=f"evidence/{d.name}/{path}"
                update_best(near,key,rank,rec)
    for k in arms.values():
        k["excess_rate_when_measurable"] = ratio(k["psf_excess_5sigma"],k["spatial_statistic_measurable"])
        k["repeat_excess_rate_for_first_excess_with_repeat_measurement"] = ratio(
            k["both_shots_excess_5sigma"],k["first_excess_and_measurable_repeat"])
    chem_sorted = [record for _,record in sorted(chemicals.values(),key=lambda entry:entry[0],reverse=True)[:12]]
    return {"release":d.name, "count_unit":"vetted observed hypothesis rows within this release unless specified; rows and releases can overlap",
            "archive_checks":dict(archive), "archive_note":"Counts are final per-row responses. Missing records can reflect a skipped query; empty_cone requires a successful cone search. Forced cutout measurements may have either sign.",
            "spatial_arms":arms,"repeat_support":{"forced_per_vetted_row":dict(rep), **repeat_controls(d)},
            "open_reasons":{ch:dict(c.most_common()) for ch,c in reasons.items()},
            "near_misses":[x[1] for x in near.values()],
            "near_miss_description":"Externally reviewed unresolved examples with actual positive measurements, at most one per distinct spatial arm, EW, or unmatched-at-host channel. Each remains blocked by the listed reasons.",
            "heii_chemistry":{"selection":"All observed vetted He II hypotheses. Rank only unresolved, unflagged 1.9–2.4 He II peaks with S/N >=4.5, FWHM <900 km/s, Ly-alpha S/N >=5 and He II / Ly-alpha >=0.1. No measured metal component may reach 3 sigma; both components of C IV or O III] must have usable limits. For each ion the HETDEX limit sums the stored positive-flux-plus-3-sigma component bounds. Redshift, chemistry, and physical origin remain unconfirmed.",
                              "stages_rows":dict(chem_stages),"selected_release_source_ids":len(chemicals),
                              "ranked_unresolved":chem_sorted},
            "controls":control_summary(d)}


def repeat_controls(d):
    meta = read(d/"repeat-summary.json")
    rows_by_stream = defaultdict(Counter)
    by_group = defaultdict(set)
    by_peak = defaultdict(set)
    strong = defaultdict(set)
    for _, p in ledger(d/"repeat-leads.jsonl.gz"):
        stream = p["stream"]
        visits = {m["shotid"]:num(m.get("snr"),-math.inf) for m in p["repeat_measurements"]}
        group = tuple(p["repeat_rows"])
        key = (group,p["wave"])
        by_group[stream].add(group)
        by_peak[stream].add(key)
        if sum(s >= 4 for s in visits.values()) >= 2:
            strong[stream].add(key)
        rows_by_stream[stream][p["channel"]] += 1
    stacks = {}
    for stream, label in (("repeat_stack","physical_groups"),("wrong_sky","control_groups")):
        stacks[stream] = {"groups_searched":meta["counts"].get(label,0),
            "distinct_groups_with_a_lead":len(by_group[stream]),"distinct_supported_peaks":len(by_peak[stream]),
            "peaks_with_2_shots_ge4":len(strong[stream]),"hypothesis_rows_by_channel":dict(rows_by_stream[stream])}
    return {"blind_stacks":stacks,
            "definition":"Peaks meeting the spectral channel screen in stacks and independently reaching 3 sigma in at least two distinct shots. Each distinct group and air wavelength is counted once across channel labels. Wrong-sky contains real spectra at different positions.",
            "source":"repeat-leads.jsonl.gz and repeat-summary.json"}


def control_summary(d):
    known = read(d/"positive-controls.json")
    weak = read(d/"weak-heii-summary.json")
    halo = read(d/"halo-injections-summary.json")
    blank, weak_population = Counter(), Counter()
    for _, trial in ledger(d/"offsource.jsonl.gz"):
        m = trial["input"]
        blank["cubes_attempted"] += 1
        blank["cubes_with_blank_apertures"] += int(bool(m.get("coordinates")))
        blank["cubes_with_error"] += int(bool(m.get("error")))
        blank["blank_apertures"] += len(m.get("coordinates") or [])
        for ch in ("heii", "ew", "unknown"):
            blank[ch+"_spectral_preselection_rows"] += trial["counts"].get(ch,0)
    for _, p in ledger(d/"weak-heii-pairs.jsonl.gz"):
        stream = "shuffled" if p.get("row_is_synthetic") else "observed" if p["pair_offset_A"] == 0 else "displaced"
        weak_population[stream] += 1
    inj = defaultdict(Counter)
    for _, p in ledger(d/"injections.jsonl.gz"):
        counts = inj[p["family"]]
        counts["attempts"] += 1
        if p["status"] != "injected":
            counts["unusable"] += 1
        elif p["was_blended"]:
            counts["preexisting_blend"] += 1
        else:
            counts["unblended_trials"] += 1
            counts["recovered_blind_primary"] += int(bool(p.get("recovered")))
            if p["family"] == "lya":
                if num(p.get("input_ew0"),0) > 240:
                    counts["input_ew_gt240"] += 1
                    counts["input_ew_gt240_and_selected"] += int(bool(p.get("selected_ew")))
            else:
                counts["selected_heii_precursor_and_lya_no_civ"] += int(bool(p.get("selected_heii_pair")))
    pos = keep(known,"targets position_or_id_matches target_blind_heii_peaks target_heii_screen_recoveries")
    pos["per_published_source"] = [keep(k,"published_detectid matched detected_near_heii recovered_heii_screen") for k in known["controls"]]
    return {"positive_controls":pos, "weak_pairs_ledger_rows":dict(weak_population),
            "weak_pair_calibration":keep(weak,"observed_clean_pairs unique_sources expected_displaced_pairs_at_true_exposure known_eight_with_rows known_eight_blind_recovered injection_bins"),
            "injections_by_family":dict(inj),"offsource":dict(blank),
            "halo":keep(halo,"cube_attempts usable_position_wavelength_pairs bins"),
            "interpretation":"Spectral injection recoveries are conditional on the input extractions. Halo trials share blank positions and are not independent sky realizations.",
            "sources":["positive-controls.json","weak-heii-summary.json","weak-heii-pairs.jsonl.gz","injections.jsonl.gz","offsource.jsonl.gz","halo-injections-summary.json"]}


if __name__ == "__main__":
    import sys
    json.dump(summarize_evidence(sys.argv[1]),sys.stdout,allow_nan=False,indent=2)
    print()
