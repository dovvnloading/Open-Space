"""Compare frozen RVBank orbital predictions to ESO DRS nights after the RVBank release."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ingest import AUTHOR, PRIMARY
from report import control_table
from runtime import require_survey_workflow
from search import blocks_for, nightly, nuisance, orbit_basis, TAU


def key(value):
    # Retain signs, catalog prefixes, and all stellar component letters.
    name = re.sub(r"[\s_*]", "", str(value)).upper()
    match = re.fullmatch(r"(HD|HIP|GJ|GL|HR)0*(\d+)([A-Z]*)", name)
    return (match[1] + str(int(match[2])) + match[3]) if match else name


class Positions:
    def __init__(self, table):
        valid = np.isfinite(table[["ra_deg", "dec_deg"]].to_numpy(dtype=float)).all(axis=1)
        self.names = table.index.to_numpy()[valid]
        self.vectors = np.asarray([self.xyz(r.ra_deg, r.dec_deg) for r in table.loc[valid].itertuples()])
        self.tree = cKDTree(self.vectors) if len(self.names) else None

    @staticmethod
    def xyz(ra, dec):
        a, d = np.radians([ra, dec])
        return np.array([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)])

    def within(self, xyz):
        if self.tree is None:
            return []
        return self.names[self.tree.query_ball_point(xyz, np.radians(3/3600))].tolist()


def identify(row, bank, eso, bank_sky, eso_sky):
    src = bank.loc[row.star]
    if not np.isfinite([src.ra_deg, src.dec_deg]).all():
        return "", np.nan, "missing primary position"
    xyz = Positions.xyz(src.ra_deg, src.dec_deg)
    others = eso_sky.within(xyz)
    if len(others) != 1:
        return "", np.nan, "no unique ESO target within 3 arcsec"
    match = eso.loc[others[0]]
    target_xyz = Positions.xyz(match.ra_deg, match.dec_deg)
    if bank_sky.within(target_xyz) != [row.star]:
        return "", np.nan, "ambiguous reverse position match"
    sep = 3600 * np.degrees(2*np.arcsin(np.clip(np.linalg.norm(target_xyz-xyz)/2, 0, 1)))
    names = {key(n) for n in [row.star, row.identity, *row.aliases.split("|")] if str(n).strip()}
    if key(match.name) not in names:
        return "", sep, "ESO identifier lacks an exact target or vetted alias match"
    return match.name, sep, ""


def unit_shape(t, s):
    basis = orbit_basis(t-s.tperi_bjd, 1/s.period_d, s.ecc, 0.)
    return basis @ [np.cos(s.omega_rad), -np.sin(s.omega_rad)]


def design(nights, others):
    t = nights.bjd.to_numpy()
    cols = [nuisance(t, nights.era.to_numpy())]
    cols.extend(unit_shape(t, s) for s in others)
    return np.column_stack(cols)


def project(a, values):
    return values - a @ np.linalg.lstsq(a, values, rcond=None)[0]


def forecast_uncertainty(train, test, sig, others, nuisance_new, info):
    # Local derivatives use frequency in resolution elements, log K, eccentricity
    # vector and mean longitude at the discovery midpoint (regular for e -> 0).
    t = train.bjd.to_numpy()
    mid, span = np.mean(t), np.ptp(t)
    anchor = (TAU*(mid-sig.tperi_bjd)/sig.period_d + sig.omega_rad) % TAU
    params = np.array([0., np.log(sig.k_ms), sig.ecc*np.cos(sig.omega_rad),
                       sig.ecc*np.sin(sig.omega_rad), anchor])
    def predict(times, par):
        df, logk, h, j, longitude = par
        omega = np.arctan2(j, h)
        basis = orbit_basis(times-mid, 1/sig.period_d + df/span, np.hypot(h, j), longitude-omega)
        return np.exp(logk) * (basis @ [np.cos(omega), -np.sin(omega)])
    def derivatives(times):
        cols = []
        for step in np.eye(5)*1e-5:
            cols.append((predict(times, params+step)-predict(times, params-step))/(2e-5))
        return np.column_stack(cols)

    w = 1/train.err.to_numpy()
    means = design(train, others) * w[:, None]
    d = np.column_stack([means, derivatives(t) * w[:, None]])
    scales = np.linalg.norm(d, axis=0)
    if np.any(scales == 0) or len(t) <= len(scales)+2:
        return np.nan
    d = d / scales
    if np.linalg.cond(d) > 1e7:
        return np.nan
    bread = np.linalg.inv(d.T@d)
    resid = project(means, (train.rv.to_numpy() - predict(t, params))*w)
    white = bread * max(1., resid@resid/(len(t)-len(scales)))
    groups = blocks_for(t, train.era.to_numpy())
    scores = np.zeros((groups.max()+1, len(scales)))
    lev = np.sum((d @ bread)*d, axis=1)
    np.add.at(scores, groups, d*(resid/np.maximum(1-lev, .01))[:, None])
    clustered = (scores @ bread).T @ (scores @ bread)
    dnew = derivatives(test.bjd.to_numpy()) / test.err.to_numpy()[:, None]
    dnew = project(nuisance_new, dnew/scales[-5:])
    # Express uncertainty in the predicted, testable waveform in units of its norm.
    variances = [np.sum((dnew @ c[-5:, -5:])*dnew) for c in (white, clustered)]
    return np.sqrt(max(0., *variances)/(sig.k_ms**2*info))


def randomization_p(scores, seed):
    """One-direction score against 14-day sign-symmetric null blocks."""
    n = len(scores)
    obs = scores.sum()
    tol = 1e-12 * max(1., np.sum(abs(scores)))
    exact = n <= 16
    draws = 2**n if exact else 65535
    rng = np.random.default_rng(int.from_bytes(seed, 'little'))
    exceed = 0
    for start in range(0, draws, 512):
        size = min(512, draws-start)
        if exact:
            signs = ((np.arange(start, start+size)[:, None] >> np.arange(n)) & 1)*2-1
        else:
            signs = rng.integers(0, 2, (size, n), dtype=np.int8)*2-1
        exceed += int(np.count_nonzero(signs @ scores >= obs-tol))
    return exceed/draws if exact else (exceed+1)/(draws+1)


def evaluate(sig, train, held, others, family):
    t, rv, err = (held[k].to_numpy() for k in ["bjd", "rv", "err"])
    result = {"n_later_nights":len(t), "later_span_d":float(np.ptp(t)) if len(t) else 0}
    if len(t) < 12:
        return dict(result, outcome="fewer than 12 later nights")
    z = blocks_for(t, held.era.to_numpy())
    phase = ((t-sig.tperi_bjd)/sig.period_d) % 1
    gap = np.diff(np.sort(np.r_[phase, phase.min()+1])).max()
    result.update(phase_gap=gap, blocks=int(z.max()+1), later_first_bjd=t.min(), later_last_bjd=t.max())
    if np.ptp(t) < sig.period_d or gap>.35 or z.max()+1<8:
        return dict(result, outcome="insufficient independent phase or time coverage")
    a = design(held, others)/err[:, None]
    if np.linalg.cond(a)>1e7 or len(t) < a.shape[1]+6:
        return dict(result, outcome="later nuisance model ill-conditioned")
    x = project(a, unit_shape(t,sig)/err)
    y = project(a, rv/err)
    info = x@x
    if info <= 1e-10*np.sum(1/err**2):
        return dict(result, outcome="orbital template degenerate with nuisances")
    share = np.bincount(z, weights=x*x)/info
    if share.max() > .40:
        return dict(result, outcome="more than 40% of orbital information in one block")
    khat = x@y/info
    resid = y - x*khat
    full = np.column_stack([a, x])
    q = np.linalg.qr(full, mode="reduced")[0]
    leverage = np.clip(np.sum(q*q, axis=1),0,.99)
    sandwich = np.linalg.norm(np.bincount(z,weights=x*resid/(1-leverage)))/info
    se = max(np.sqrt(max(1.,(resid@resid)/(len(t)-full.shape[1]))/info), sandwich)
    p = randomization_p(np.bincount(z, weights=x*y), hashlib.sha256(
        (sig.star+"|"+str(sig.signal_index)).encode()).digest()[:8])
    delta = float(2*sig.k_ms*(x@y)-sig.k_ms**2*info)
    predictive = forecast_uncertainty(train, held, sig, others, a, info)
    result.update(k_later_signed_ms=khat, k_later_se_ms=se, fixed_minus_null_improvement=delta,
                  single_test_p=p, bonferroni_p=min(1.,p*family),
                  later_block_max_share=share.max(), prediction_rel_uncertainty=predictive)
    if not np.isfinite(predictive) or predictive>.35:
        return dict(result, outcome="discovery waveform extrapolation uncertain")
    sigma = np.hypot(se, predictive*sig.k_ms)
    if p*family<=.05 and delta>0 and abs(khat-sig.k_ms)<=2.58*sigma:
        return dict(result, outcome="fixed waveform corroborated on later dates")
    if abs(khat-sig.k_ms)>3*sigma:
        return dict(result, outcome="later signed amplitude discrepant with fixed prediction")
    return dict(result, outcome="later comparison inconclusive")


def read(path):
    return pd.read_csv(path, dtype={"star":str}, keep_default_na=False, low_memory=False)


def execute(work, out, old):
    primary_meta = json.loads((out/"ingest.json").read_text())
    if primary_meta["source"] not in (PRIMARY, AUTHOR):
        (out/"holdout.json").write_text(json.dumps({"status":"requires primary RVBank discovery run"})+"\n")
        return
    bank = read(out/"targets.csv").set_index("star",verify_integrity=True)
    eso = read(old.parent/"targets.csv").set_index("star",verify_integrity=True)
    for table in (bank,eso):
        for col in ("ra_deg","dec_deg"):
            table[col] = pd.to_numeric(table[col], errors="coerce")
    if set(eso.source) != {"ESO Phase 3 HARPS RV catalogue DR1"}:
        raise ValueError("Holdout requires ESO Phase 3 target provenance")
    bank_sky, eso_sky = Positions(bank), Positions(eso)
    train = pd.read_csv(work/"epochs.csv.gz",dtype={"star":str},low_memory=False)
    # Use the release-wide last dated exposure (including rejected RVs). Exclude
    # its entire local night to eliminate timing-convention overlaps at the edge.
    all_bjd = pd.read_csv(out/"epoch-ledger.csv.gz",usecols=["bjd"]).bjd
    good = all_bjd.between(2452700,2465000)
    last = float(all_bjd[good].max())
    if not np.isfinite(last) or last < primary_meta["last_bjd"]-1.e-6:
        raise ValueError("Cannot verify release-wide date boundary")
    future = pd.read_csv(old,dtype={"star":str},low_memory=False)
    future = future[future.bjd.ge(np.floor(last)+1)].groupby("star",sort=False)
    train = train.groupby("star",sort=False)
    signals = read(out/"signals.csv")
    vets = read(out/"vetting.csv")
    all_signals = signals.merge(vets[["star","signal_index","identity","aliases","external_status"]],
                                on=["star","signal_index"],validate="one_to_one")
    controls = control_table(read(out/"stars.csv"), signals, vets)
    labels = defaultdict(list)
    for ref in controls.itertuples():
        if str(ref.outcome).startswith("period recovered"):
            hits = signals[signals.star.eq(ref.matching_targets) & signals.period_d.eq(ref.fit_period_d)]
            if len(hits)==1:
                labels[(ref.matching_targets,int(hits.iloc[0].signal_index))].append(ref.host+ref.companion)
    eligible = (all_signals.rv_status.isin(["pass","inconclusive"]) &
                all_signals.fap.le(.01) & all_signals.delta_bic.ge(10))
    eligible |= [(r.star,int(r.signal_index)) in labels for r in all_signals.itertuples()]
    selected = all_signals[eligible]
    rows, night_tables = [], []
    for star, subset in selected.groupby("star",sort=True):
        identity = subset.iloc[0]
        match, sep, reason = identify(identity,bank,eso,bank_sky,eso_sky)
        held = None
        if match in future.groups:
            exposures = future.get_group(match).copy()
            held = nightly(exposures)
            exposures["day"] = np.floor(exposures.bjd).astype(int)
            ids = exposures.groupby(["day","era"]).epoch_id.agg(lambda g:json.dumps(list(g)))
            packet = held.merge(ids.rename("epoch_ids").reset_index(), on=["day","era"],
                                validate="one_to_one")
            packet.insert(0,"eso_id",match)
            packet.insert(0,"star",star)
            night_tables.append(packet)
        orig = nightly(train.get_group(star))
        components = all_signals[all_signals.star.eq(star) & all_signals.fap.le(.01) & all_signals.delta_bic.ge(10)]
        for sig in subset.itertuples(index=False):
            row = dict(star=star,signal_index=sig.signal_index,period_d=sig.period_d,k_ms=sig.k_ms,
                       fap=sig.fap,rv_status=sig.rv_status,external_status=sig.external_status,
                       discovery_flags=sig.flags,eso_id=match,sep_arcsec=sep,
                       controls="|".join(labels[(star,int(sig.signal_index))]),n_later_nights=0)
            if reason:
                row["outcome"] = reason
            elif held is None:
                row["outcome"] = "no ESO exposure after release-wide cutoff"
            else:
                others = components[components.signal_index.ne(sig.signal_index)].itertuples(index=False)
                try:
                    row.update(evaluate(sig, orig, held, list(others), len(selected)))
                except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                    row.update(outcome="numerical comparison unresolved",detail=str(exc))
            rows.append(row)
    table = pd.DataFrame(rows,columns=None if rows else
                         ["star","signal_index","period_d","k_ms","controls","eso_id","n_later_nights","outcome"])
    table.to_csv(out/"holdout.csv",index=False)
    packets = (pd.concat(night_tables, ignore_index=True) if night_tables else
               pd.DataFrame(columns=["star","eso_id","day","era","bjd","rv","err","epoch_ids"]))
    packets.to_csv(out/"holdout-nights.csv.gz",index=False,
                   compression={"method":"gzip","mtime":0})
    controls["later_result"] = [next((row["outcome"] for row in rows
                               if ref.host+ref.companion in row["controls"].split("|")),
                              "no matching discovery template") for ref in controls.itertuples()]
    controls.to_csv(out/"holdout-controls.csv",index=False)
    meta = dict(status="completed",latest_primary_bjd=last,first_allowed_bjd=np.floor(last)+1,
                tested_family=len(selected),outcomes=table.outcome.value_counts().to_dict(),
                input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                              [old,old.parent/"targets.csv",old.parent/"inputs.csv",old.parent/"ingest.json",
                               work/"epochs.csv.gz",out/"signals.csv",out/"vetting.csv"]})
    (out/"holdout.json").write_text(json.dumps(meta, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--eso-epochs",type=Path,required=True)
    args = parser.parse_args()
    require_survey_workflow()
    execute(args.work,args.out,args.eso_epochs)
