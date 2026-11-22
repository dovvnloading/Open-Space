"""Post-discovery whole-pulsar sign randomization and telescope-split diagnostics."""
import argparse
import csv
import gzip
import itertools
import json
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy.linalg import cho_factor, solve_triangular
from scipy.special import log_ndtr
from scipy.stats import beta, chi2

import array_search as ar
from search import fit_pulsar

DRAWS = 4095
SEED = 21615073
FIELDS = ("family", "shape", "epoch_mjd", "sky_id")


def key(row):
    sky = row["sky_id"]
    return row["family"], row["shape"], float(row["epoch_mjd"]), None if sky in ("", None) else int(sky)


def identifier(event):
    return "-".join(str(x) for x in key(event) if x is not None)


def read(path):
    return json.loads(path.read_text())


def log_score(q, d):
    """Minus log10 of a chi-square tail, including tails below floating-point range."""
    q = np.maximum(q, 0.)
    log_one = np.log(2.) + log_ndtr(-np.sqrt(q))
    if d == 1:
        log_p = log_one
    elif d == 2:
        log_p = -q / 2.
    elif d == 3:
        with np.errstate(divide="ignore"):
            log_p = np.logaddexp(log_one, (np.log(2 * q / np.pi) - q) / 2)
    else:
        raise ValueError("unexpected response dimension")
    return np.maximum(-log_p / np.log(10.), 0.)


class SignBank:
    def __init__(self, signs, targets):
        self.signs = signs
        self.position = {name: i for i, name in enumerate(targets)}
        self.kernels = defaultdict(list)
        self.maxima = np.zeros(signs.shape[0])
        self.tests = 0
        self.observed = 0.

    def add(self, model, cell):
        fit = model.estimate(cell.z, cell.factor * cell.rss)
        # Whitening the low-dimensional coefficient estimate makes its squared
        # norm exactly the archived statistic for every vector of pulsar signs.
        kernel = solve_triangular(np.linalg.cholesky(fit.covariance),
                                  model.matrix * cell.z, lower=True, check_finite=False)
        dim = len(fit.coefficients)
        q = np.sum(kernel.sum(axis=1)**2)
        if not np.isclose(q, fit.statistic, atol=1e-8, rtol=2e-6):
            raise ValueError("sign kernel does not reproduce the observed statistic")
        indices = [self.position[n] for n in cell.names]
        if dim not in self.kernels:
            for signs in self.signs[:3,indices]:
                direct = model.estimate(cell.z*signs, cell.factor*cell.rss).statistic
                if not np.isclose(np.linalg.norm(kernel @ signs)**2,direct,atol=1e-8,rtol=2e-6):
                    raise ValueError("randomized sign kernel disagrees with direct fit")
        full = np.zeros((dim, len(self.position)))
        full[:, indices] = kernel
        self.kernels[dim].append(full)
        self.tests += 1
        self.observed = max(self.observed, float(log_score(q, dim)))
        return fit

    def finish(self):
        # The chi-square tail is monotone within each response dimension, so
        # evaluate that tail only for the largest Q in each group in each draw.
        for dim, matrices in self.kernels.items():
            operator = np.stack(matrices).reshape(-1, len(self.position))
            for start in range(0, len(self.maxima), 64):
                stop = min(start+64, len(self.maxima))
                b = (operator @ self.signs[start:stop].T).reshape(len(matrices), dim, stop-start)
                q = np.einsum("kdi,kdi->ki", b, b).max(axis=0)
                self.maxima[start:stop] = np.maximum(self.maxima[start:stop], log_score(q, dim))
        self.kernels.clear()

    def tail(self, score, comparisons):
        # Count floating-point ties as exceedances. The exact binomial bounds
        # quantify Monte Carlo error conditional on these data and this statistic.
        tolerance = 1e-9 * (1 + abs(score))
        count = int(np.count_nonzero(self.maxima >= score-tolerance))
        tail = .025 / max(comparisons, 1)
        lo = float(beta.ppf(tail, count, DRAWS-count+1)) if count else 0.
        hi = float(beta.isf(tail, count+1, DRAWS-count)) if count < DRAWS else 1.
        return {"exceedances": count, "p_sign_max_bank": (count+1)/(DRAWS+1),
                "mc_low": lo, "mc_high": hi}


def add_fit(bank, cell, ref, event, sky, capture):
    model = ar.Contrast(cell, event["family"], sky, ref)
    fit = bank.add(model, cell)
    if capture:
        return fit, model.matrix
    return fit, None


def bank_checks(data, evidence, summary):
    files = sorted((data/"nb").glob("*.npz"))
    epochs, start, span = ar.grid(files)
    records, errors = ar.read_records(files, epochs, start, span, False)
    names = [p.name for p in records]
    if names != summary["participating_pulsars"]:
        raise ValueError("coherence inputs do not reproduce the original participating pulsars")
    passed = set()
    for name in names:
        path = evidence / "single" / (name + ".json")
        if path.is_file():
            single = read(path)
            if single["status"] == "ok" and single["noise_diagnostics"]["passes"]:
                passed.add(name)
    signs = np.random.default_rng(SEED).choice([-1., 1.], size=(DRAWS, len(names)))
    signs[:, 0] = 1.  # Simultaneously reversing every sign leaves all Q unchanged.
    full, clean = SignBank(signs, names), SignBank(signs, names)
    captured = {key(e): {} for e in summary["survivors"]}
    skies = ar.sky_grid()
    epoch_index = {float(e): k for k, e in enumerate(epochs)}
    max_diff = 0.
    with gzip.open(evidence/"array"/"array_trials.csv.gz", "rt") as table:
        rows = csv.DictReader(table)
        for (shape, epoch), group in itertools.groupby(rows, key=lambda r: (r["shape"], float(r["epoch_mjd"]))):
            c = ar.cell(records, ar.SHAPES.index(shape), epoch_index[epoch], len(records))
            ref = cho_factor(c.common + np.eye(len(c.z)), check_finite=False)
            mask = np.array([name in passed for name in c.names])
            sub = c.subset(mask) if mask.sum() >= ar.MIN_PSR else None
            ref_sub = cho_factor(sub.common + np.eye(len(sub.z)), check_finite=False) if sub else None
            for row in group:
                index = key(row)
                capture = index in captured
                sky = None if index[3] is None else skies[index[3]]
                fit, matrix = add_fit(full, c, ref, row, sky, capture)
                original = float(row["statistic"])
                max_diff = max(max_diff, abs(original-fit.statistic))
                if not np.isclose(original, fit.statistic, atol=1e-7, rtol=2e-5):
                    raise ValueError(f"archived trial mismatch {index}: {original} vs {fit.statistic}")
                if capture:
                    captured[index].update(full=fit, matrix=matrix, cell=c)
                if sub:
                    try:
                        fit_sub, _ = add_fit(clean, sub, ref_sub, row, sky, capture)
                    except (ValueError, np.linalg.LinAlgError) as exc:
                        if capture:
                            captured[index]["subset_status"] = str(exc)
                        continue
                    if capture:
                        captured[index].update(subfit=fit_sub, subset_names=sub.names, subset_status="tested")
                elif capture:
                    captured[index]["subset_status"] = "fewer than eight eligible noise-diagnostic-passing pulsars"
    if full.tests != summary["n_trials"] or any("full" not in row for row in captured.values()):
        raise ValueError("original bank or survivors incompletely reconstructed")
    full.finish()
    clean.finish()
    return full, clean, captured, passed, start, span, max_diff, errors


def describe_events(survivors, fits, full, clean, passed):
    events, pulsars = [], []
    # Missing external tests occupy their place in this family with p=1.
    raw_p = np.asarray([e["epta_check"].get("p_fixed_one_sided", 1.) for e in survivors])
    order = np.argsort(raw_p, kind="stable")
    holm = np.empty(len(raw_p))
    holm[order] = np.minimum(1., np.maximum.accumulate(raw_p[order] * (len(raw_p)-np.arange(len(raw_p)))))
    for ordinal, event in enumerate(survivors, start=1):
        tag = identifier(event)
        row = dict(survivor_ordinal=ordinal, event=tag, **{k:event[k] for k in FIELDS}, amplitude=event["amplitude"],
                   units=event["units"], published_bank_p_joint=event["p_joint"],
                   internal_check=event["internal_check"]["status"],
                   wideband_check=event["wideband_check"]["status"], epta_check=event["epta_check"]["status"])
        saved = fits[key(event)]
        f, c, matrix = saved["full"], saved["cell"], saved["matrix"]
        d = len(f.coefficients)
        direction = f.coefficients / np.linalg.norm(f.coefficients)
        alpha_columns = 2*len(survivors)+2
        row.update(full.tail(float(log_score(f.statistic,d)), alpha_columns))
        row.update(full_n=len(c.names), subset_status=saved["subset_status"],
                   subset_n=sum(n in passed for n in c.names))
        scale, _ = ar.physical(event["shape"],event["family"])
        ref_amp = float(np.linalg.norm(f.coefficients))*scale
        ref_se = float(np.sqrt(direction @ f.covariance @ direction))*scale
        row.update(discovery_direction_amplitude=ref_amp, discovery_direction_se=ref_se,
                   wideband_direction_amplitude=event["wideband_check"].get("signed_amplitude"),
                   wideband_direction_se=event["wideband_check"].get("se"),
                   epta_tested_pulsars=event["epta_check"].get("n_pulsars",0))
        ext = event["epta_check"]
        if "p_fixed_one_sided" in ext:
            if ext["units"] != event["units"]:
                raise ValueError("independent comparison units differ from discovery")
            delta = ext["signed_amplitude"] - ref_amp
            row.update(epta_p_fixed=ext["p_fixed_one_sided"],
                       epta_p_holm=float(holm[ordinal-1]),
                       epta_p_bonferroni=min(1.,len(survivors)*ext["p_fixed_one_sided"]),
                       epta_direction_amplitude=ext["signed_amplitude"], epta_direction_se=ext["se"],
                       epta_minus_nb_amplitude=delta, epta_difference_in_quadrature_se=delta/np.hypot(ext["se"],ref_se))
        row["pass_both_sign_and_subset_checks"] = False
        if "subfit" in saved:
            s = saved["subfit"]
            row["subset_amplitude_in_discovery_direction"] = float(direction @ s.coefficients)*scale
            row["subset_p_parametric_original_bank"] = ar.p_global(s.p, full.tests)
            row.update({"subset_"+k:v for k,v in
                        clean.tail(float(log_score(s.statistic, d)), alpha_columns).items()})
            row["pass_both_sign_and_subset_checks"] = bool(
                row["mc_high"] <= ar.ALPHA and row["subset_mc_high"] <= ar.ALPHA
                and row["subset_amplitude_in_discovery_direction"] > 0
                and row["subset_p_parametric_original_bank"] <= ar.ALPHA)
        events.append(row)
        g = ar.patterns(c.xyz, event["family"], None if key(event)[3] is None else ar.sky_grid()[key(event)[3]])[0]
        weights = (direction @ matrix)*scale
        contributions = weights * c.z
        if not np.isclose(contributions.sum(), np.linalg.norm(f.coefficients)*scale, atol=0., rtol=1e-7):
            raise ValueError(f"{tag}: pulsar contributions do not sum to fitted norm")
        total = np.abs(contributions).sum()
        for i, name in enumerate(c.names):
            pulsars.append(dict(survivor_ordinal=ordinal, event=tag, pulsar=name, passes_noise_diagnostics=name in passed,
                                scalar_score=float(c.z[i]), amplitude_weight=float(weights[i]),
                                signed_contribution=float(contributions[i]), units=event["units"],
                                absolute_contribution_fraction=float(abs(contributions[i])/total),
                                pulsar_delay=float(c.z[i]/c.weight[i]), predicted_delay=float(g[i] @ f.coefficients),
                                delay_units="us/yr" if event["shape"] == "ramp" else "us"))
    return events, pulsars


def system_participation(data, survivors, pulsars):
    by_event = {identifier(row): row for row in survivors}
    by_pulsar = defaultdict(list)
    for r in pulsars:
        by_pulsar[r["pulsar"]].append(r)
    out = []
    for pulsar, rows in by_pulsar.items():
        p = fit_pulsar(data/"nb"/(pulsar+".npz"), fit_soft=False)
        for r in rows:
            event = by_event[r["event"]]
            dt = p.t - event["epoch_mjd"]
            tpl = p.whiten(ar.waveform(event["shape"], dt))
            z = (tpl/np.linalg.norm(tpl))*p.ry
            if not np.isclose(z.sum(), r["scalar_score"], atol=1e-7, rtol=2e-5):
                raise ValueError(f"{pulsar}: backend score partition mismatch")
            window = 3*float(event["shape"][8:]) if event["shape"].startswith("gaussian") else 180.
            for system in np.unique(p.backend):
                mask = p.backend == system
                before, after = dt[mask & (dt <= 0)], dt[mask & (dt >= 0)]
                out.append(dict(survivor_ordinal=r["survivor_ordinal"], event=r["event"], pulsar=pulsar, backend=system,
                                signed_contribution=float(z[mask].sum()*r["amplitude_weight"]),
                                units=r["units"], reduced_rows=int(mask.sum()),
                                days_in_window=len(np.unique(np.floor(p.t[mask & (abs(dt)<=window)]))),
                                support_window_days=window,
                                days_to_nearest_before=float(-before.max()) if len(before) else None,
                                days_to_nearest_after=float(after.min()) if len(after) else None))
    return out


def telescope_records(data, release, epochs, start, span, names):
    """Repeat the original fixed-template projection after isolating each telescope."""
    records, errors = defaultdict(list), []
    row_fields = ("t","y","err","freq","backend","design","ecorr_group","ecorr_us")
    for path in sorted((data/release).glob("*.npz")):
        if path.stem not in names:
            continue
        try:
            records["ALL"].append(ar.reduce_one(path, epochs, start, span, False))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            errors.append({"release":release,"pulsar":path.stem,"site":"ALL","error":str(exc)})
            continue
        if release == "wb":
            continue
        with np.load(path, allow_pickle=False) as source:
            arrays = {k:source[k] for k in row_fields if k in source}
            xyz = {k:source[k] for k in ("ra","dec")}
        # PINT has separate observatory codes for clocks at the same telescope.
        # LEAP remains separate as it coherently combined several observatories.
        aliases = {"jbroach":"jodrell","jbdfb":"jodrell","jbafb":"jodrell",
                   "jodrell_pre_2021":"jodrell","ncyobs":"nancay",
                   "effelsberg_pre_2021":"effelsberg","effelsberg_asterix":"effelsberg",
                   "arecibo_pre_2021":"arecibo","vla_site":"vla",
                   "jb_mkii_dfb":"jb_mkii","jb_mkii_rch":"jb_mkii"}
        raw_sites = [str(system).split("/",1)[0] for system in arrays["backend"]]
        sites = np.array([aliases.get(site,site) for site in raw_sites])
        for site in np.unique(sites):
            select = sites == site
            with TemporaryDirectory() as tmp:
                part = Path(tmp)/path.name
                np.savez(part, **xyz, **{k:v[select] for k,v in arrays.items()})
                try:
                    records[str(site)].append(ar.reduce_one(part, epochs, start, span, False))
                except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                    errors.append({"release":release,"pulsar":path.stem,"site":str(site),"error":str(exc)})
    return records, errors


def telescope_checks(data, survivors, start, span):
    out, errors = [], []
    if not survivors:
        return out, errors
    epochs = np.unique([e["epoch_mjd"] for e in survivors])
    names = set(n for e in survivors for n in e["pulsars"])
    skies = ar.sky_grid()
    for release in ("nb","wb","epta"):
        sites, skipped = telescope_records(data, release, epochs, start, span, names)
        errors.extend(skipped)
        for ordinal, event in enumerate(survivors,start=1):
            j = ar.SHAPES.index(event["shape"])
            k = int(np.searchsorted(epochs,event["epoch_mjd"]))
            sky = None if key(event)[3] is None else skies[key(event)[3]]
            for site in sorted(set(sites) | {"ALL"}):
                records = sites.get(site, [])
                included = [p for p in records if p.name in event["pulsars"]]
                c = ar.cell(included,j,k,max(1,len(records)),minimum=1)
                result = ar.replication(c,event,sky)
                r = {field:value for field,value in result.items() if field != "per_pulsar"}
                r.update(survivor_ordinal=ordinal, event=identifier(event), release=release, telescope=site,
                         units=event["units"], input_pulsars=len(included))
                if "signed_amplitude" in r:
                    r["fraction_of_discovery_amplitude"] = r["signed_amplitude"] / np.linalg.norm(event["components"])
                if c is not None and len(c.names) >= ar.MIN_PSR:
                    try:
                        f = ar.trial(c, "clock", None)
                        r.update(clock_amplitude=float(f.coefficients[0]),
                                 clock_se=float(np.sqrt(f.covariance[0,0])), clock_p_local=f.p,
                                 clock_units="us/yr" if event["shape"] == "ramp" else "us")
                    except (ValueError, np.linalg.LinAlgError) as exc:
                        r["clock_status"] = str(exc)
                out.append(r)
    return out, errors


def write_table(path, rows, fallback):
    fields = list(dict.fromkeys(k for row in rows for k in row)) or fallback
    ar.write_csv(path, rows, fields)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("data","evidence","out"):
        parser.add_argument("--"+arg, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary = read(args.evidence/"array/array_summary.json")
    full, clean, fits, passed, start, span, discrepancy, errors = bank_checks(args.data, args.evidence, summary)
    events, pulsars = describe_events(summary["survivors"], fits, full, clean, passed)
    systems = system_participation(args.data,summary["survivors"],pulsars)
    telescopes, skipped = telescope_checks(args.data,summary["survivors"],start,span)
    maxima = [dict(draw=i+1, full_max_log10_score=float(a),
                   subset_max_log10_score=float(b) if clean.tests else None)
              for i,(a,b) in enumerate(zip(full.maxima,clean.maxima))]
    write_table(args.out/"null-maxima.csv", maxima, ["draw"])
    write_table(args.out/"events.csv", events, ["event"])
    write_table(args.out/"pulsars.csv", pulsars, ["event","pulsar"])
    write_table(args.out/"systems.csv", systems, ["event","pulsar","backend"])
    write_table(args.out/"telescope_checks.csv", telescopes, ["event","release","telescope","status"])
    result = {
        "draws":DRAWS, "seed":SEED, "original_survivors":len(events),
        "n_bank":full.tests, "n_subset_bank":clean.tests, "diagnostic_passing_pulsars":sorted(passed),
        "largest_original_statistic_replay_difference":discrepancy,
        "observed_full_max_log10_score":full.observed,
        "observed_full_max_tail":full.tail(full.observed,2*len(events)+2),
        "observed_subset_max_log10_score":clean.observed if clean.tests else None,
        "observed_subset_max_tail":clean.tail(clean.observed,2*len(events)+2) if clean.tests else None,
        "n_pass_both_sign_and_subset_checks":sum(e["pass_both_sign_and_subset_checks"] for e in events),
        "input_failures":errors, "telescope_subset_failures":skipped,
        "interpretation":(
            "Post-discovery checks of the original survivors. One random sign applies to every residual of each "
            "narrowband pulsar across all epochs, shapes and sky directions. Within each draw the maximum negative "
            "log10 chi-square tail across the bank compares one-, two- and three-parameter responses. "
            "The randomization tests joint invariance under independent pulsar sign changes. A stochastic "
            "Hellings-Downs background, shared clocks or other cross-pulsar noise can violate this assumption. "
            "These checks cannot identify a deterministic burst. Monte Carlo p values use (exceedances+1)/(draws+1); "
            "the exact binomial bounds jointly cover the reported randomization tails at 95% conditional on the data. "
            "The diagnostic subset keeps only original pulsars passing the single-pulsar noise screen; each retained "
            "template still needs the original eight-pulsar, coverage and geometry criteria. Telescope comparisons "
            "freeze the discovery epoch, sky, sign and polarization and report descriptive fixed-template tails."
            " EPTA directional probabilities from the original replication fits are adjusted across every "
            "internally surviving array template, counting unavailable tests as p=1. EPTA minus NANOGrav "
            "amplitudes use the original fitted vector direction and the same physical units. The difference "
            "in quadrature standard errors is descriptive; sharing DE440 and any true common process violates "
            "the assumption of independent errors. Candidate selection also affects discovery amplitudes."
        ),
    }
    (args.out/"summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(f"{len(events)} original array survivors checked with {DRAWS} joint sign draws", flush=True)


if __name__ == "__main__":
    main()
