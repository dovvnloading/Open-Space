"""Blind ACIS-I/S timed-exposure search in level-2 sky/time event lists.

Coordinates in this module are 1-origin physical sky pixels or CCD pixels.
Time is the EVENTS TIME column in mission seconds (FITS time reference below).
No network access, source list, or target coordinates enter the trigger search.
"""

from dataclasses import dataclass
import math
import re
import warnings

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from astropy.wcs import FITSFixedWarning, WCS
from scipy.spatial import cKDTree
from scipy.stats import binom

iers.conf.auto_download = False

DURATIONS = (60, 180, 600, 1800, 3600)
CELL_WIDTHS = (8, 24)
MIN_COUNTS = 8
SEED_P = 1e-7
PIXEL_ARCSEC = 0.492
EVENT_DTYPE = np.dtype([(key, "f8") for key in ("t", "x", "y", "e")] +
                       [(key, "i4") for key in ("cx", "cy", "ccd")])
EV_TO_ERG = 1.602176634e-12


def _union(intervals):
    merged = []
    for start, stop in sorted(intervals):
        if not np.isfinite([start, stop]).all() or stop <= start:
            raise ValueError("Bad GTI endpoints")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(stop, merged[-1][1]))
        else:
            merged.append((start, stop))
    return merged


def _time_support(spans):
    """Disjoint union weighted by the greatest simultaneous live fraction."""
    changes = {}
    for start, stop, duty in spans:
        for point, delta in ((start, 1), (stop, -1)):
            changes.setdefault(point, []).append((duty, delta))
    active = {}
    out = []
    previous = None
    for point, updates in sorted(changes.items()):
        if active and previous is not None:
            factor = max(active)
            if out and out[-1][1] == previous and out[-1][2] == factor:
                out[-1][1] = point
            else:
                out.append([previous, point, factor])
        for duty, delta in updates:
            active[duty] = active.get(duty, 0) + delta
            if not active[duty]:
                del active[duty]
        previous = point
    return out


def _livetime_union(spans):
    return float(sum((b-a)*d for a,b,d in _time_support(spans)))


def _gti_map(hdul):
    by_chip = {}
    common = []
    for hdu in hdul:
        hdr = hdu.header
        if "GTI" not in (hdu.name + str(hdr.get("HDUCLAS1", ""))).upper():
            continue
        if getattr(hdu, "columns", None) is None:
            continue
        names = {n.upper(): n for n in hdu.columns.names}
        if not {"START", "STOP"} <= names.keys():
            continue
        chip = hdr.get("CCD_ID")
        # Standard Chandra GTI<n> is the chip-specific data-subspace block.
        if chip is None:
            match = next((m for name in (hdu.name, hdr.get("HDUNAME", ""))
                          if (m := re.fullmatch(r"GTI[_-]?([0-9])", name.strip().upper()))), None)
            if match:
                chip = int(match[1])
        if chip is None and "CCD_ID" not in names:
            common.extend(zip(hdu.data[names["START"]], hdu.data[names["STOP"]]))
            continue
        chips = hdu.data[names["CCD_ID"]] if "CCD_ID" in names else [chip]*len(hdu.data)
        for c, a, b in zip(chips, hdu.data[names["START"]], hdu.data[names["STOP"]]):
            by_chip.setdefault(int(c), []).append((float(a), float(b)))
    if by_chip and common:
        # A separately supplied common GTI restricts every chip.
        common = _union(common)
        by_chip = {c: [(max(a, x), min(b, y)) for a, b in intervals for x, y in common
                       if min(b, y) > max(a, x)] for c, intervals in by_chip.items()}
    if not by_chip and common:
        raise ValueError("Cannot attribute common GTI to active CCDs; coverage unverified")
    if not by_chip:
        raise ValueError("No chip-specific GTIs")
    return {c: _union(intervals) for c, intervals in by_chip.items()}


def _read_events(hdu):
    keys = {n.upper(): n for n in hdu.columns.names}
    required = {"X", "Y", "TIME", "CHIPX", "CHIPY", "CCD_ID", "ENERGY", "STATUS", "GRADE"}
    if required - keys.keys():
        raise ValueError("Missing event columns " + str(sorted(required - keys.keys())))
    rows = []
    for left in range(0, len(hdu.data), 150_000):
        block = hdu.data[left:left + 150_000]
        energy = block[keys["ENERGY"]]
        status = block[keys["STATUS"]]
        clean = (status == 0).reshape(len(block), -1).all(axis=1)
        clean &= np.isin(block[keys["GRADE"]], [0, 2, 3, 4, 6])
        clean &= (energy >= 500) & (energy <= 7000)
        for name in ("TIME", "X", "Y"):
            clean &= np.isfinite(block[keys[name]])
        arr = np.empty(clean.sum(), dtype=EVENT_DTYPE)
        for key, source in zip(arr.dtype.names, ("TIME", "X", "Y", "ENERGY", "CHIPX",
                                                  "CHIPY", "CCD_ID")):
            arr[key] = block[keys[source]][clean]
        rows.append(arr)
    return np.concatenate(rows) if rows else np.empty(0, dtype=EVENT_DTYPE)


def _sky_wcs(hdu):
    names = [name.upper() for name in hdu.columns.names]
    inds = [names.index(name) + 1 for name in ("X", "Y")]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        sky = WCS(hdu.header, keysel=["pixel"], colsel=inds).celestial
    if sky.pixel_n_dim != 2 or sky.wcs.lng != 0 or sky.wcs.lat != 1:
        raise ValueError("Unrecognized celestial binary-table X,Y WCS")
    return sky


def _orientation(data):
    # Median-time-neighbour difference removes telescope dither to first order.
    # CHIP->SKY rotation is required for the sky-aligned readout-column test.
    if len(data) < 8:
        return None
    s = np.argsort(data["t"], kind="stable")
    if len(s) > 30000:
        k = np.linspace(0, len(s)-2, 30000, dtype=int)
        left, right = s[k], s[k+1]
    else:
        left, right = s[:-1], s[1:]
    det = np.column_stack([data[k][right]-data[k][left] for k in ("cx", "cy")])
    sky = np.column_stack([data[k][right]-data[k][left] for k in ("x", "y")])
    valid = (np.linalg.norm(det, axis=1) > 8) & (data["t"][right]-data["t"][left] < 60)
    if valid.sum() < 5:
        return None
    matrix, _, rank, _ = np.linalg.lstsq(det[valid], sky[valid], rcond=None)
    scales = np.linalg.svd(matrix, compute_uv=False)
    if rank != 2 or not np.all((scales > .85) & (scales < 1.15)):
        return None
    return matrix


@dataclass
class Segment:
    ccd: int
    start: float
    stop: float
    duty: float
    frame: float
    orientation: object


def _spatial_layout(data):
    """Cache sky tilings across permutations, independent of event times."""
    grids = []
    if len(data) < MIN_COUNTS:
        return grids
    for width in CELL_WIDTHS:
        for dx, dy in ((0, 0), (.5, 0), (0, .5), (.5, .5)):
            gx = np.floor(data["x"] / width + dx).astype("i8")
            gy = np.floor(data["y"] / width + dy).astype("i8")
            keys = (gx-gx.min())*(gy.max()-gy.min()+1) + gy-gy.min()
            # A CCD footprint covers up to ~1500x1500 native sky pixels after
            # roll and dither. Compact labels also tolerate unusual footprints.
            _, cells, counts = np.unique(keys, return_inverse=True, return_counts=True)
            eligible = counts >= MIN_COUNTS
            if not eligible.any():
                continue
            labels = np.full(len(counts), -1, dtype=np.int32)
            labels[eligible] = np.arange(eligible.sum())
            grids.append((width, labels[cells], counts[eligible]))
    return grids


def _proposals(data, seg, grids):
    """Trials start at each photon, in four overlapping sky tilings per width."""
    length = seg.stop - seg.start
    t = data["t"] - seg.start
    time_order = np.argsort(t, kind="stable")
    proposals = []
    for width, cells, size in grids:
        group = cells[time_order]
        eligible = group >= 0
        group = group[eligible]
        # Stable integer sort groups photons, retaining increasing time in
        # each group. This avoids repeating eight two-key sorts per scramble.
        by_cell = np.argsort(group, kind="stable")
        order = time_order[eligible][by_cell]
        group = group[by_cell]
        st = t[order]
        combined = group.astype(np.float64)*(length + max(DURATIONS) + 2) + st
        index = np.arange(len(order))
        for duration in DURATIONS:
            if length < 3*duration:
                continue
            n = np.searchsorted(combined, combined+duration, side="left") - index
            p0 = duration/length
            inside = ((st >= duration) & (st <= length-2*duration) &
                      (n >= MIN_COUNTS) & (n >= 4*p0*size[group]))
            good = np.flatnonzero(inside)
            if not len(good):
                continue
            ps = binom.sf(n[good]-1, size[group[good]], p0)
            # Non-maximum suppression within a given sky cell and duration.
            # Sort by cell then p, and skip overlapping, weaker intervals.
            chosen = np.flatnonzero(ps <= SEED_P)
            chosen = chosen[np.lexsort((st[good[chosen]], ps[chosen], group[good[chosen]]))]
            active_group = -1
            occupied = []
            for i in chosen:
                k = good[i]
                if group[k] != active_group:
                    occupied = []
                    active_group = group[k]
                tt = float(st[k] + seg.start)
                if any(abs(tt-last) < duration for last in occupied):
                    continue
                occupied.append(tt)
                ids = order[k:k+n[k]]
                xs, ys = data["x"][ids], data["y"][ids]
                x, y = np.median(xs), np.median(ys)
                distances = np.hypot(xs-x, ys-y)
                if distances.min() > .70*width:
                    # A split or corner-heavy cell can have an empty median
                    # aperture. Anchor at its nearest contributing event.
                    nearest = distances.argmin()
                    x, y = xs[nearest], ys[nearest]
                proposals.append((float(ps[i]), tt, duration, width, float(x), float(y)))
    return sorted(proposals)


def _upper(n_on, n_off, exp_on, exp_off):
    return float(binom.sf(n_on-1, n_on+n_off, exp_on/(exp_on+exp_off)))


def _assess(proposal, data, seg, tree, meta, wcs):
    p_seed, t0, window, width, x, y = proposal
    t1 = t0 + window
    radius = width * .70
    indices = np.asarray(tree.query_ball_point((x, y), 4*radius), dtype=int)
    dr = np.hypot(data["x"][indices]-x, data["y"][indices]-y)
    it = data["t"][indices]
    inside = (it >= t0) & (it < t1)
    ap = dr <= radius
    idx = indices[inside & ap]
    k = len(idx)
    if not k:
        raise ValueError("Empty aperture on a nonempty seed")
    before = int(np.count_nonzero(ap & (it < t0)))
    after = int(np.count_nonzero(ap & (it >= t1)))
    local_before = int(np.count_nonzero(ap & (it >= t0-window) & (it < t0)))
    local_after = int(np.count_nonzero(ap & (it >= t1) & (it < t1+window)))
    bkg = int(np.count_nonzero(inside & (dr >= 2*radius)))
    # Annulus r=2..4 has 12 aperture areas. Reject truncated geometries below.
    annulus_ratio = 12
    source = data[idx]
    t = source["t"]
    ps = _upper(k, bkg, 1, annulus_ratio)
    pb = _upper(k, before, window, t0-seg.start)
    pa = _upper(k, after, window, seg.stop-t1)
    vetoes = []
    if k < MIN_COUNTS or ps > 1e-5:
        vetoes.append("weak_local_spatial_excess")
    if (max(pb, pa) > .025 or
        max(before/(t0-seg.start), after/(seg.stop-t1)) > .25*k/window or
        max(local_before, local_after) > max(1, .35*k)):
        vetoes.append("active_or_unconstrained_flank")
    whole_on = (data["t"] >= t0) & (data["t"] < t1)
    # Compare CCD-wide rate with the rest of the GTI, masking the source+annulus.
    far_on = int(whole_on.sum() - inside.sum())
    far_off = int(len(data) - len(indices) - far_on)
    outtime = seg.stop-seg.start-window
    pccd = _upper(far_on, far_off, window, outtime)
    ann_off = int(np.count_nonzero((~inside) & (dr >= 2*radius)))
    p_ann = _upper(bkg, ann_off, window, outtime)
    if bkg/window > 2*ann_off/outtime and p_ann < .01:
        vetoes.append("local_background_flare")
    if (far_on/window > 1.8*far_off/outtime and pccd < .001):
        vetoes.append("ccd_background_flare")
    detector_pixel = source["cx"].astype("i8")*2048+source["cy"]
    _, hits = np.unique(detector_pixel, return_counts=True)
    _, frame_hits = np.unique(t, return_counts=True)
    if hits.max() >= max(3, .25*k) or len(hits) < 4:
        vetoes.append("repeated_detector_pixel")
    if np.unique(source["cx"], return_counts=True)[1].max() >= max(5, .70*k):
        vetoes.append("single_detector_column")
    if len(frame_hits) < 4 or frame_hits.max() >= max(4, .5*k):
        vetoes.append("single_or_few_frames")
    if (max(np.unique(source["cy"], return_counts=True)[1]) >= max(5, .65*k)
        and seg.ccd == 8):
        vetoes.append("s4_row_streak")
    if (min(source["cx"].min(),source["cy"].min()) < 4*radius or
        max(source["cx"].max(),source["cy"].max()) > 1025-4*radius):
        vetoes.append("chip_edge_or_truncated_annulus")
    # Subarray footprints can end well inside the 1024x1024 detector.
    for name in ("cx","cy"):
        if (np.median(source[name]) - data[name].min() < 4*radius or
                data[name].max() - np.median(source[name]) < 4*radius):
            if "chip_edge_or_truncated_annulus" not in vetoes:
                vetoes.append("chip_edge_or_truncated_annulus")
    strip_n = None
    strip_side = None
    streak_p = None
    if seg.orientation is None:
        vetoes.append("unresolved_readout_axis")
    else:
        # Test a putative readout streak in sky coordinates: bright points and
        # dispersed events on the same column outside the extraction aperture.
        axis = seg.orientation[1] / np.linalg.norm(seg.orientation[1])
        rx, ry = data["x"]-x, data["y"]-y
        along = np.abs(axis[0]*rx+axis[1]*ry)
        across = np.abs(axis[1]*rx-axis[0]*ry)
        length_mask = (along >= 3*radius) & (along < 800)
        half = max(3, radius*.6)
        core = length_mask & (across < half)
        sides = length_mask & (across >= 2*half) & (across < 4*half)
        strip_n, strip_side = int(core.sum()), int(sides.sum())
        streak_p = _upper(strip_n, strip_side, 1, 2)
        if strip_n >= 15 and strip_n > strip_side and streak_p < .001:
            vetoes.append("readout_column_excess")
    xcen, ycen = (float(np.mean(source[field])) for field in ("x", "y"))
    ra, dec = map(float, wcs.all_pix2world([[xcen, ycen]], 1)[0])
    if not np.isfinite([ra, dec]).all():
        raise ValueError("Nonfinite celestial coordinates")
    spread = np.hypot(source["x"]-xcen, source["y"]-ycen)
    if k >= 12:
        covariance = np.cov(np.stack([source["x"],source["y"]]))
        eig = np.linalg.eigvalsh(covariance)
        if eig[1] > 16*max(eig[0],.1):
            vetoes.append("elongated_event_track")
    r90 = math.hypot(1.2, 1.52*PIXEL_ARCSEC*float(np.sqrt(np.mean(spread**2)/k)))
    peak_length = min(60, window)
    # Maximum raw counts in a fixed 60-second interval fully inside the seed.
    sorted_t = np.sort(t)
    starts = np.unique(np.clip(sorted_t, t0, t1-peak_length))
    peak_n = np.searchsorted(sorted_t, starts+peak_length)-np.searchsorted(sorted_t, starts)
    best = int(np.argmax(peak_n))
    pt = float(starts[best])
    peak_pick = (t >= pt) & (t < pt+peak_length)
    ann_peak = (dr >= 2*radius) & (it >= pt) & (it < pt+peak_length)
    peak_bkg = int(ann_peak.sum())
    energy_sum = source["e"][peak_pick].sum() - data["e"][indices[ann_peak]].sum()/annulus_ratio
    qt = np.quantile(t, [.05,.95])
    if qt[1]-qt[0] < 40:
        vetoes.append("subminute_event_concentration")
    live_peak = peak_length*seg.duty
    mjd = meta["mjdref"] + (meta["timezero"] + pt)/86400
    mjd_utc = float(Time(mjd, format="mjd", scale=meta["timesys"]).utc.mjd)
    report = dict(obsid=meta["obsid"], ccd_id=seg.ccd, ra_deg=ra, dec_deg=dec,
                  sky_x=xcen, sky_y=ycen, aperture_x=x, aperture_y=y, aperture_radius_px=radius,
                  position_error_arcsec=r90, event_mjd=mjd_utc, event_mjd_scale="UTC",
                  start_met_s=t0, stop_met_s=t1, trial_duration_s=window,
                  duration_s=float(qt[1]-qt[0]), duration_definition="aperture photon t95-t05",
                  on_counts=k, annulus_counts=bkg, annulus_area_ratio=annulus_ratio,
                  before_counts=before, after_counts=after, before_s=t0-seg.start,
                  after_s=seg.stop-t1, adjacent_before_counts=local_before,
                  adjacent_after_counts=local_after, seed_local_p=p_seed,
                  seed_score=-math.log10(max(p_seed,1e-300)), live_fraction=seg.duty,
                  spatial_local_p=ps, before_local_p=pb, after_local_p=pa,
                  annulus_off_counts=ann_off, annulus_flare_local_p=p_ann,
                  ccd_on_counts=far_on, ccd_off_counts=far_off, ccd_background_local_p=pccd,
                  max_detector_pixel_counts=int(hits.max()), distinct_detector_pixels=len(hits),
                  distinct_frames=len(frame_hits), max_counts_per_frame=int(frame_hits.max()),
                  readout_strip_counts=strip_n, readout_side_counts=strip_side,
                  readout_local_p=streak_p, peak_start_met_s=pt, peak_window_s=peak_length,
                  peak_live_s=live_peak, peak_counts=int(peak_n[best]),
                  peak_annulus_counts=peak_bkg,
                  peak_rate_ct_s=float((peak_n[best]-peak_bkg/annulus_ratio)/live_peak),
                  peak_detected_energy_erg_s=float(energy_sum*EV_TO_ERG/live_peak),
                  peak_flux_erg_cm2_s=None, instrument_pass=not vetoes, vetoes=vetoes)
    return report, idx


def _same_trigger(a, b):
    r = max(a["aperture_radius_px"], b["aperture_radius_px"])
    return (min(a["stop_met_s"],b["stop_met_s"]) > max(a["start_met_s"],b["start_met_s"])
            and np.hypot(a["sky_x"]-b["sky_x"],a["sky_y"]-b["sky_y"]) <= r)


def _scan(data, seg, meta, wcs, grids=None):
    if grids is None:
        grids = _spatial_layout(data)
    proposals = _proposals(data, seg, grids)
    if not proposals:
        return []
    tree = cKDTree(np.column_stack([data["x"], data["y"]]))
    found = []
    for proposal in proposals:
        event, indices = _assess(proposal, data, seg, tree, meta, wcs)
        existing = next((i for i,(other,_) in enumerate(found) if _same_trigger(other,event)), None)
        if existing is None:
            found.append((event, indices))
        elif event["instrument_pass"] and not found[existing][0]["instrument_pass"]:
            found[existing] = (event, indices)
    return found


def _inject(data, seg, rng, duration, counts, sigma):
    # Uniform interior detector positions; detector counts, not incident photons.
    # Fold the fitted CHIP/SKY orientation through the observed dither offsets.
    if seg.orientation is None or not len(data):
        return None
    first, last = seg.start+duration, seg.stop-2*duration
    if first >= last:
        return None
    t0 = rng.uniform(first,last)
    order = np.argsort(data["t"])
    times = data["t"][order]
    target = rng.uniform(t0,t0+duration,counts)
    def nearest(t):
        right = np.clip(np.searchsorted(times,t),0,len(times)-1)
        left = np.maximum(0,right-1)
        return order[np.where(np.abs(times[right]-t)<np.abs(times[left]-t),right,left)]
    ref = nearest(target)
    # Work from the nearest observed frame to retain local frame-grid phase.
    t = data["t"][ref] + np.rint((target-data["t"][ref])/seg.frame)*seg.frame
    middle = nearest(t0+duration/2)
    chip = np.column_stack([data["cx"],data["cy"]])
    sky = np.column_stack([data["x"],data["y"]])
    low, high = np.min(chip,axis=0)+32, np.max(chip,axis=0)-32
    if (low >= high).any():
        return None
    xy = (rng.uniform(low,high)-chip[middle]) @ seg.orientation + sky[middle]
    coords = xy+rng.normal(0,sigma,(counts,2))
    cxy = np.rint((coords-sky[ref]) @ np.linalg.inv(seg.orientation)+chip[ref]).astype(int)
    use = ((cxy >= chip.min(axis=0)) & (cxy <= chip.max(axis=0))).all(axis=1)
    syn = np.empty(counts,dtype=EVENT_DTYPE)
    syn["t"], syn["x"], syn["y"] = t, coords[:,0], coords[:,1]
    syn["cx"], syn["cy"], syn["ccd"] = cxy[:,0], cxy[:,1], seg.ccd
    syn["e"] = rng.choice(data["e"],counts)
    return np.concatenate([data,syn[use]]), xy, t0, int(use.sum())


def search_file(path, seed=20261121, n_controls=4, n_injections=4):
    """Return JSON-serializable complete per-EVT2 results; fail closed on invalid metadata."""
    if min(n_controls,n_injections,seed) < 0:
        raise ValueError("Counts and seed must be nonnegative")
    with fits.open(path, memmap=False) as hdul:
        evt = hdul["EVENTS"]
        h = evt.header
        if (str(h.get("INSTRUME","")).strip().upper() != "ACIS" or
            str(h.get("GRATING","")).strip().upper() != "NONE" or
            str(h.get("READMODE","")).strip().upper() != "TIMED"):
            raise ValueError("Requires ACIS timed-exposure, non-grating EVT2")
        if u.Unit(h.get("TIMEUNIT","s")) != u.s:
            raise ValueError("Unsupported event time unit")
        gtis = _gti_map(hdul)
        data = _read_events(evt)
        wcs = _sky_wcs(evt)
        mjdref = h.get("MJDREF")
        if mjdref is None:
            mjdref = float(h["MJDREFI"]) + float(h["MJDREFF"])
        meta = dict(obsid=str(h["OBS_ID"]), mjdref=float(mjdref),
                    timesys=h["TIMESYS"].lower().strip(), timezero=float(h.get("TIMEZERO",0)))
        frame = float(h["TIMEDEL"])
        if frame <= 0 or not np.isfinite(frame):
            raise ValueError("Invalid frame duration")
        spans = []
        segments = []
        for c, intervals in sorted(gtis.items()):
            ontime, livetime = h.get(f"ONTIME{c}"), h.get(f"LIVTIME{c}")
            duty = (float(livetime)/float(ontime) if ontime and livetime is not None
                    else float(h["DTCOR"]))
            if not 0 < duty <= 1.00001:
                raise ValueError(f"Invalid live fraction for CCD {c}")
            duty = min(1., duty)
            chip = data[data["ccd"]==c]
            matrix = _orientation(chip)
            for a,b in intervals:
                a,b = float(a),float(b)
                spans.append((a,b,duty))
                if b-a <= 3*min(DURATIONS):
                    continue
                segment_data = chip[(chip["t"] >= a) & (chip["t"] < b)]
                segments.append((Segment(c,a,b,duty,frame,matrix),segment_data))
    n_clean = len(data)
    del data, evt, hdul, chip
    result = dict(obsid=meta["obsid"], seed=int(seed), gti_livetime_s=_livetime_union(spans),
                  cleaned_events=int(n_clean), triggers=[], controls=[], injections=[],
                  start_support_by_duration={}, exposure_by_duration_s={}, ccd_exposure_s={})
    for duration in DURATIONS:
        allowed = [(s.start+duration,s.stop-2*duration,s.duty) for s,_ in segments
                   if s.stop-s.start > 3*duration]
        support = _time_support(allowed)
        result["start_support_by_duration"][str(duration)] = support
        result["exposure_by_duration_s"][str(duration)] = float(sum((b-a)*d for a,b,d in support))
    result["exposure_s"] = result["exposure_by_duration_s"][str(min(DURATIONS))]
    for c in gtis:
        result["ccd_exposure_s"][str(c)] = float(sum((s.stop-s.start-3*min(DURATIONS))*s.duty
                                                     for s,_ in segments if s.ccd==c))
    controls = [dict(replica=i,triggers=[]) for i in range(n_controls)]
    for j,(seg,events) in enumerate(segments):
        grids = _spatial_layout(events)
        result["triggers"].extend(row for row,_ in _scan(events,seg,meta,wcs,grids))
        for i,control in enumerate(controls):
            rng = np.random.default_rng(np.random.SeedSequence([seed,1,j,i]))
            shuffled = events.copy()
            shuffled["t"] = rng.permutation(events["t"])
            control["triggers"].extend(row for row,_ in _scan(shuffled,seg,meta,wcs,grids))
    result["controls"] = controls
    for i in range(n_injections):
        rng = np.random.default_rng(np.random.SeedSequence([seed,2,i]))
        duration = int(rng.choice(DURATIONS))
        counts = int(rng.choice([8,16,32,64]))
        sigma = float(rng.choice([1.2,3.,6.]))
        trial = dict(trial=i,duration_s=duration,counts=counts,sigma_sky_px=sigma,
                     attempted=False,recovered=None)
        result["injections"].append(trial)
        weights = [max(0.,s.stop-s.start-3*duration)*s.duty for s,_ in segments]
        if not sum(weights):
            trial["reason"] = "no eligible CCD time for this duration"
            continue
        seg,events = segments[rng.choice(len(segments),p=np.asarray(weights)/sum(weights))]
        added = _inject(events,seg,rng,duration,counts,sigma)
        if added is None:
            trial["reason"] = "insufficient detector geometry"
            continue
        augmented, xy, start, added_counts = added
        matches = []
        for event, ids in _scan(augmented,seg,meta,wcs):
            if (event["instrument_pass"] and np.count_nonzero(ids >= len(events)) >= 4
                and np.hypot(event["sky_x"]-xy[0],event["sky_y"]-xy[1]) <= event["aperture_radius_px"]
                and min(event["stop_met_s"],start+duration)>max(event["start_met_s"],start)):
                matches.append(event)
        trial.update(attempted=True,recovered=bool(matches),ccd_id=seg.ccd,sky_x=float(xy[0]),
                     sky_y=float(xy[1]),start_met_s=float(start),inserted_counts=added_counts)
    return result
