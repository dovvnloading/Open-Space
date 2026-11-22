"""Compare a stable daily candidate and its one-day aliases on the native orbit products.

Every fit is a follow-up on data also used by the daily search. Improvements and
window correlations diagnose alias ambiguity; they are not independent detections,
false-alarm probabilities or barycentric ephemerides.
"""
import argparse
import gzip
import json
from pathlib import Path
import re
from time import perf_counter
from urllib.parse import urljoin

import numpy as np

from ingest import ascii_curve, bat_fits, fetch, page_links


def product_urls(source, band):
    if band == "bat":
        base = source["bat_page"].rstrip("/")
        yield base + ".orbit.lc.fits"
        yield base + ".orbit.lc.txt"
        pages = [source["bat_page"]]
    else:
        pages = source.get("maxi_pages", [source["maxi_page"]])
        for page in pages:
            if "/star_data/" in page:
                yield urljoin(page, source["maxi_id"] + "_g_lc_1orb_all.dat")
    for page in pages:
        try:
            data, meta = fetch(page)
        except Exception:
            continue
        for _, url, label in page_links(data, meta["url"]):
            if band == "maxi" and re.search(r"1orb.*all\.(dat|txt)(\.gz)?$", url, re.I):
                yield url
            if band == "bat" and re.search(r"orbit\.lc\.(fits|txt)(\.gz)?$", url, re.I):
                yield url


def load_native(source, band):
    attempts, seen = [], set()
    for url in product_urls(source, band):
        if url in seen:
            continue
        seen.add(url)
        try:
            data, receipt = fetch(url)
            if data.startswith(b"\x1f\x8b"):
                data = gzip.decompress(data)
            if data.startswith(b"SIMPLE"):
                t, y, e, flag = bat_fits(data)
            else:
                t, y, e, flag = ascii_curve(data, band)
            good = np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (e > 0) & (t > 50000) & (t < 85000)
            if flag is not None:
                good &= flag == 0
            t, y, e = t[good], y[good], e[good]
            # Exact duplicates cannot supply additional independent weight.
            # The row with the smallest quoted uncertainty is deterministic.
            ix = np.lexsort((e, t))
            t, y, e = t[ix], y[ix], e[ix]
            keep = np.r_[True, np.diff(t) != 0] if len(t) else np.zeros(0, dtype=bool)
            t, y, e = t[keep], y[keep], e[keep]
            if len(t) < 60 or np.median(np.diff(t)) >= 0.9:
                raise ValueError("insufficient subdaily timing coverage")
            receipt.update(rows_raw=len(good), rows_used=len(t), duplicated_epochs=int((~keep).sum()),
                           quality_column=flag is not None, first_mjd=float(t[0]), last_mjd=float(t[-1]),
                           median_separation_days=float(np.median(np.diff(t))), attempts=attempts)
            return (t, y, e), receipt
        except Exception as exc:
            attempts.append({"url":url, "error":str(exc)[:320]})
    return None, {"status":"unavailable", "attempts":attempts}


class BlockFit:
    def __init__(self, t, y, dy, block_days):
        self.t = t - t[0]
        self.w = 1 / np.square(dy)
        self.span = float(self.t[-1])
        self.block = np.floor(self.t / block_days).astype(int)
        self.wsum = np.bincount(self.block, weights=self.w)
        self.invsum = np.divide(1, self.wsum, out=np.zeros_like(self.wsum), where=self.wsum > 0)
        self.y = self.project(y)
        self.wy = self.w * self.y
        self.chi0 = float(self.y @ self.wy)
        # All models share the local frequency grid shape. Cache the offset
        # rotations once per epoch so just two trig evaluations per native
        # timestamp are needed for each predicted alias family centre.
        off = np.arange(-4, 5) / (4 * self.span)
        self.rotations = [(v, np.cos(2*np.pi*v*self.t), np.sin(2*np.pi*v*self.t)) for v in off]

    def project(self, col):
        sums = np.bincount(self.block, weights=self.w * col)
        return col - (sums * self.invsum)[self.block]

    def columns(self, f):
        angle = (2*np.pi*f) * self.t
        return np.stack([self.project(np.cos(angle)), self.project(np.sin(angle))], axis=1)

    def improvement(self, x):
        g = x.T @ (self.w[:, None] * x)
        v = x.T @ self.wy
        b = np.linalg.pinv(g, rcond=1e-10, hermitian=True) @ v
        return float(max(0, v @ b)), b

    def profile(self, group):
        # Include exact prediction. Common ±1/T local search for each model,
        # including the target, is intended to accommodate reference-grid error.
        angle = 2*np.pi*group["frequency"] * self.t
        c0, s0 = np.cos(angle), np.sin(angle)
        best = None
        for offset, cv, sv in self.rotations:
            f = group["frequency"] + offset
            if f <= 0:
                continue
            c, s = c0*cv - s0*sv, s0*cv + c0*sv
            wc, ws = self.w*c, self.w*s
            bc, bs = np.bincount(self.block, wc), np.bincount(self.block, ws)
            mc, ms = self.invsum*bc, self.invsum*bs
            # Block-diagonal intercept elimination through the weighted normal
            # equations avoids constructing an N-by-number-of-blocks matrix.
            cross = wc@s - mc@bs
            g = np.array([[wc@c - mc@bc, cross], [cross, ws@s - ms@bs]])
            v = np.array([c@self.wy, s@self.wy])
            b = np.linalg.pinv(g, rcond=1e-10, hermitian=True) @ v
            improve = float(max(0, v@b))
            if best is None or improve > best["delta_chi2"]:
                best = dict(frequency=float(f), period_days=float(1/f), delta_chi2=improve,
                            semi_amplitude=float(np.hypot(*b)), rank=int(np.linalg.matrix_rank(g)))
        return {**best, "predicted_frequency":group["frequency"], "family":group["family"],
                "descriptions":group["descriptions"]}

    def degeneracy(self, reference, contender):
        x = self.columns(reference["frequency"])
        z = self.columns(contender["frequency"])
        # Canonical correlations of the two sampled sine/cosine subspaces,
        # after removing separate weighted means in every time block.
        sw = np.sqrt(self.w)[:, None]
        xw, zw = x*sw, z*sw
        def invroot(g):
            ev, vec = np.linalg.eigh(g)
            use = ev > np.max(ev) * 1e-10
            return (vec[:, use] / np.sqrt(ev[use])) @ vec[:, use].T
        cc = invroot(xw.T @ xw) @ (xw.T @ zw) @ invroot(zw.T @ zw)
        both = self.improvement(np.column_stack((x, z)))[0]
        f, g = reference["frequency"], contender["frequency"]
        def window(a):
            return float(abs(self.w @ np.exp(2j*np.pi*a*self.t)) / self.w.sum())
        return {"canonical_correlations":np.clip(np.linalg.svd(cc, compute_uv=False), 0, 1).tolist(),
                "window_at_frequency_difference":window(f-g), "window_at_frequency_sum":window(f+g),
                "alias_gain_given_target":max(0, both-reference["delta_chi2"]),
                "target_gain_given_alias":max(0, both-contender["delta_chi2"])}


def predictions(f):
    groups = {}
    def add(freq, family, descr):
        if freq <= 0:
            return
        group = groups.setdefault(round(freq, 13), {"frequency":freq, "family":family, "descriptions":[]})
        group["descriptions"].append(descr)
    add(f, "target", "daily target")
    for h in (1, 2, 3):
        if h > 1:
            add(h*f, "harmonic", f"{h}f")
            add(f/h, "harmonic", f"f/{h}")
        for n in range(1, 7):
            for sign in (-1, 1):
                add(abs(n + sign*h*f), "one_day_alias", f"|{n} {sign:+d} * {h}f|")
                if h > 1:
                    add(abs(n + sign*f) / h, "one_day_alias", f"|{n} {sign:+d} * f|/{h}")
    return list(groups.values())


def consistency(profile):
    epochs = ["full", "early", "late"]
    if not all(profile[t]["status"] == "ok" for t in epochs):
        return {"status":"incomplete_epochs"}
    winner = profile["full"]["best_one_day_alias"]["predicted_frequency"]
    differences = {}
    for tag in epochs:
        chosen = next(r for r in profile[tag]["fits"] if r["predicted_frequency"] == winner)
        differences[tag] = profile[tag]["target"]["delta_chi2"] - chosen["delta_chi2"]
    return {"status":"ok", "full_best_alias_prediction_per_day":winner,
            "target_minus_full_winner_delta_chi2_by_epoch":differences,
            "best_alias_prediction_same_both_halves":all(
                profile[t]["best_one_day_alias"]["predicted_frequency"] == winner for t in ("early", "late")),
            "target_wins_against_all_one_day_aliases_full_and_both_halves":all(
                profile[t]["target_minus_best_alias_delta_chi2"] > 0 for t in epochs)}


def describe_fit(curve, p):
    t, y, dy = curve
    half = (t[0] + t[-1])/2
    results = {}
    for tag, use in (("full", np.ones(len(t), dtype=bool)), ("early", t < half), ("late", t >= half)):
        if use.sum() < 30 or np.ptp(t[use]) < 2*p:
            results[tag] = {"status":"insufficient_span_or_rows", "n":int(use.sum())}
            continue
        started = perf_counter()
        block = max(120., 4*p)
        fit = BlockFit(t[use], y[use], dy[use], block)
        fits = [fit.profile(group) for group in predictions(1/p)]
        ref = fits[0]
        alts = sorted((r for r in fits if r["family"] == "one_day_alias"),
                      key=lambda r:r["delta_chi2"], reverse=True)
        strongest = alts[0]
        window = fit.degeneracy(ref, strongest)
        results[tag] = {"status":"ok", "n":int(use.sum()), "first_mjd":float(t[use][0]), "last_mjd":float(t[use][-1]),
                        "baseline_days":fit.span, "nuisance_block_days":block, "occupied_blocks":int(np.sum(fit.wsum>0)),
                        "null_chi2":fit.chi0, "target":ref, "best_one_day_alias":strongest, "winner_window":window,
                        "target_minus_best_alias_delta_chi2":ref["delta_chi2"]-strongest["delta_chi2"],
                        "runtime_seconds":round(perf_counter()-started, 3), "fits":fits}
    return results


def run(args):
    args.out.mkdir(parents=True, exist_ok=True)
    roster = {r["id"]:r for r in json.loads(args.inventory.read_text())["sources"]}
    selected = []
    for line in args.candidates.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if (not args.ids or r["source_id"] in args.ids) and r.get("stability", {}).get("pass"):
                selected.append(r)
    ids = {r["source_id"] for r in selected}
    if args.ids and set(args.ids) != ids:
        raise ValueError("Each requested ID must have a daily period that passed the disjoint-half tests")
    summary = {"interpretation":"The orbit products rebin the same instrument data used for selection. "
               "Scores compare alias models and do not supply independent false-alarm probabilities. "
               "Canonical correlations near one indicate timing degeneracy. Times are archive-frame, not barycenter-corrected.",
               "native_products":{}, "tests":[]}
    for ident in sorted(ids):
        source = roster[ident]
        curves = {}
        for band in ("bat", "maxi"):
            curves[band], receipt = load_native(source, band)
            summary["native_products"].setdefault(ident, {})[band] = receipt
        for i,r in enumerate(c for c in selected if c["source_id"] == ident):
            period = r["period_days"]
            small = {"source_id":ident, "period_days":period, "instruments":{}}
            diagnostics = {}
            for band, curve in curves.items():
                if curve is None:
                    small["instruments"][band] = {"status":"product_unavailable"}
                    continue
                profile = describe_fit(curve, period)
                diagnostics[band] = profile
                small["instruments"][band] = {"half_consistency":consistency(profile)}
                for epoch, fit in profile.items():
                    small["instruments"][band][epoch] = {k:v for k,v in fit.items() if k != "fits"}
            name = f"{ident}-{i}.json.gz"
            with gzip.open(args.out/name, "wt") as out:
                json.dump(diagnostics, out, separators=(",", ":"), allow_nan=False)
            small["all_predicted_fits"] = name
            summary["tests"].append(small)
            print(ident, period, "native sampling diagnostics written", flush=True)
    (args.out/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True, help="Daily survey survivor JSONL; may be restricted to unresolved periods")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ids", nargs="+", help="Restrict to these BAT identifiers. Default: all supplied stable periods.")
    run(parser.parse_args())
