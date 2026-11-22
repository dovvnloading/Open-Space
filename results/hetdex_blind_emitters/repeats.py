"""Blind independent-shot coadds and matched wrong-sky controls for HETDEX."""
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import spectral as spec


def vectors(ra, dec):
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(r)*np.cos(d), np.sin(r)*np.cos(d), np.sin(d)))


def safe_groups(cat, stats):
    """Exact-ID groups plus isolated <=1 arcsec pairs among ungrouped rows."""
    ra, dec = np.asarray(cat.info["ra"], float), np.asarray(cat.info["dec"], float)
    good = np.isfinite(ra) & (ra >= 0) & (ra < 360) & np.isfinite(dec) & (np.abs(dec) <= 90)
    good &= cat.shot_ids > 0
    xyz = vectors(ra, dec)
    used = np.zeros(cat.size, dtype=bool)
    # Mark every non-singleton ID, even if later rejected, so a problematic
    # association cannot re-enter through the positional fallback.
    starts = np.r_[0, np.flatnonzero(np.diff(cat.sorted_src_ids))+1]
    ends = np.r_[starts[1:], cat.size]
    for loc in np.flatnonzero(ends-starts >= 2):
        rows = np.sort(cat.src_order[starts[loc]:ends[loc]])
        if cat.src_ids[rows[0]] <= 0:
            continue
        used[rows] = True
        valid = rows[good[rows]]
        shots, indices = np.unique(cat.shot_ids[valid], return_index=True)
        if len(shots) < 2:
            continue
        stats["repeated_ids_ge2_shots"] += 1
        # Distances to centroid <=1 arcsec and maximum pair distance <=2".
        chord = np.linalg.norm(xyz[valid]-xyz[valid].mean(axis=0), axis=1)
        if chord.max() > np.deg2rad(1/3600):
            stats["id_groups_position_rejected"] += 1
            continue
        yield "id", valid[indices].tolist()
    singles = np.flatnonzero(~used & good)
    if len(singles) == 0:
        return
    # Query against all rows so a blended or repeat-ID neighbor can disqualify
    # a positional association. Only globally isolated 2-row pairs are used.
    ids = np.flatnonzero(good)
    tree = cKDTree(xyz[ids])
    chord1 = 2*np.sin(np.deg2rad(1/3600)/2)
    for start in range(0, len(singles), 32768):
        rows = singles[start:start+32768]
        dist, loc = tree.query(xyz[rows], k=3, distance_upper_bound=chord1)
        for row, others, distances in zip(rows, loc, dist):
            if not np.isfinite(distances[1]) or np.isfinite(distances[2]):
                continue
            pair = ids[others[:2]]
            other = int(pair[1] if pair[0] == row else pair[0])
            if row >= other or used[other] or cat.shot_ids[other] == cat.shot_ids[row]:
                continue
            # The other member must also have exactly one other neighbor.
            sd, _ = tree.query(xyz[other], k=3, distance_upper_bound=chord1)
            if np.isfinite(sd[2]):
                continue
            yield "isolated_position_pair", [int(row), other]


def spectra_for_rows(cat, rows):
    # Respect the same public-catalog missing markers as acquisition.Spectra.read.
    f = np.asarray(cat.flux[rows], dtype=float)
    e = np.asarray(cat.error[rows], dtype=float)
    ok = np.isfinite(f) & (f != -999) & np.isfinite(e) & (e > 0)
    return np.where(ok, f, 0), np.where(ok, e, np.inf)


def _stack(f, e):
    weight = 1/e**2
    norm = np.sum(weight, axis=0)
    y = np.divide(np.sum(weight*f, axis=0), norm, out=np.zeros_like(norm), where=norm>0)
    se = np.sqrt(np.divide(1., norm, out=np.full_like(norm, np.inf), where=norm>0))
    return y, se


def stack_for_rows(cat, rows):
    """Reconstruct the exact 1D stack in public FITS numeric units."""
    return _stack(*spectra_for_rows(cat, rows))


def wrong_sky(cat, rows, shot_order, sorted_shots, xyz, rng):
    """Keep the first visit. Draw each other visit from the original shot."""
    fake = [rows[0]]
    selected = {int(cat.src_ids[rows[0]])}
    for i in rows[1:]:
        shot = cat.shot_ids[i]
        lo, hi = np.searchsorted(sorted_shots, shot), np.searchsorted(sorted_shots, shot, "right")
        accepted = False
        for _ in range(64):
            j = int(shot_order[rng.integers(lo, hi)])
            # Reject source-sharing and require >=60 arcsec from real location.
            sep = np.linalg.norm(xyz[j]-xyz[rows[0]])
            if not np.isfinite(sep) or sep < np.deg2rad(60/3600) or int(cat.src_ids[j]) in selected:
                continue
            fake.append(j)
            selected.add(int(cat.src_ids[j]))
            accepted = True
            break
        if not accepted:
            return None
    return fake


def seen_in_single(out):
    matches = defaultdict(list)
    path = out/"screen.jsonl.gz"
    with gzip.open(path, "rt") as inp:
        for line in inp:
            row = json.loads(line)
            matches[row["row_index"]].append((row["channel"], row["wave"]))
    return matches


def repeat_search(cat, out):
    """Return normal channel records and save every lead and reject count."""
    from scan import CAT_KEYS, channels, dump
    out = Path(out)
    single = seen_in_single(out)
    rng = np.random.default_rng(392071 if cat.release=="sc1" else 3297201)
    shot_order = np.argsort(cat.shot_ids, kind="stable")
    sorted_shots = cat.shot_ids[shot_order]
    xyz = vectors(np.asarray(cat.info["ra"], float), np.asarray(cat.info["dec"], float))
    counts = Counter()
    result = []
    groups = safe_groups(cat, counts)
    dx = float(np.median(np.diff(cat.wave)))

    with gzip.open(out/"repeat-leads.jsonl.gz", "wt", compresslevel=3) as output:
        while True:
            block = []
            for association, indices in groups:
                rows = [int(i) for i in indices]
                counts["physical_groups"] += 1
                counts["observed_visits"] += len(rows)
                counts["association_"+association] += 1
                block.append(("repeat_stack", association, rows, rows[0]))
                wrong = wrong_sky(cat, rows, shot_order, sorted_shots, xyz, rng)
                if wrong:
                    block.append(("wrong_sky", association, wrong, rows[0]))
                    counts["control_groups"] += 1
                else:
                    counts["control_no_comparison"] += 1
                if len(block) >= 128:
                    break
            if not block:
                break
            data = [spectra_for_rows(cat, rows) for _, _, rows, _ in block]
            means, errors = zip(*[_stack(f,e) for f,e in data])
            peaks = spec.search_batch(cat.wave, means, errors, density_unit=1e-17)
            # Reuse scatter and mask evaluations for all peaks in one stack.
            calibrated = {}
            for index, peak in peaks:
                stream, assoc, rows, ref = block[index]
                counts[stream+"_peaks"] += 1
                z = float(cat.info["z_hetdex"][ref])
                z = z if np.isfinite(z) and z >= -0.01 else None
                use = channels(cat.wave, means[index], errors[index], peak, z)
                if not use:
                    continue
                counts[stream+"_channel_peaks_before_corroboration"] += 1
                if index not in calibrated:
                    f,e,ok = spec._inputs(cat.wave,*data[index],density_unit=1e-17)
                    _,_,scale = spec._screen(f,e,ok,dx)
                    calibrated[index] = (f, e*scale[:,None])
                f,e = calibrated[index]
                per_visit = []
                for i, y, se in zip(rows, f, e):
                    measurement = spec._force(cat.wave, y, se, dx, peak["wave_vac"], peak["sigma_A"])
                    per_visit.append(dict(measurement, row_index=i, shotid=int(cat.shot_ids[i]),
                                          source_id=int(cat.src_ids[i])))
                n_seen = sum(m.get("snr",0) >= 3 for m in per_visit)
                if n_seen < 2:
                    counts[stream+"_uncorroborated_channel_peaks"] += 1
                    continue
                m = cat.metadata(ref)
                meta = {k:m.get(k) for k in CAT_KEYS}
                for channel in use:
                    previous = any(c==channel and abs(w-peak["wave"]) <= max(4.,peak["sigma_A"])
                                   for i in rows for c,w in single.get(i, ()))
                    rec = {**meta,"catalog":meta,**peak,"channel":channel,
                           "stream":stream,"row_index":ref,"repeat_association":assoc,
                           "repeat_rows":rows,"repeat_measurements":per_visit,
                           "repeat_n_supported":n_seen,"repeat_new_to_single_screen":not previous}
                    counts[stream+"_supported_"+channel] += 1
                    if not previous:
                        counts[stream+"_new_to_single_"+channel] += 1
                    output.write(dump(rec)+"\n")
                    if stream == "repeat_stack":
                        result.append(rec)
    summary = {"release":cat.release,"counts":dict(counts),
               "minimum_coadded_snr":spec.MIN_SNR,"minimum_supporting_shots":2,
               "minimum_forced_per_shot_snr":3.0,"max_position_match_arcsec":1.0,
               "control_minimum_offset_arcsec":60.0,
               "frame":"native topocentric grid; no per-visit barycentric resampling",
               "scope":"within-release only; controls retain shot-specific noise and unrelated real lines"}
    (out/"repeat-summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return result
