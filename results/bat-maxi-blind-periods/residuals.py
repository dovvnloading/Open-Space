"""Targeted native-time residual and harmonic-family audit of frozen blind survivors."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import linalg, optimize

from field_checks import load_archives, read_rows
from search import ALPHA, Continuum, prepare, survey_trials

HARMONICS = 6
MAX_FAMILIES = 8


def waves(t, frequency, harmonics=1):
    a = 2*np.pi*t[:, None]*(frequency*np.arange(1, harmonics+1))
    return np.column_stack((np.cos(a),np.sin(a)))


def column_basis(matrix):
    q, r, _ = linalg.qr(matrix, mode="economic", pivoting=True, check_finite=False)
    diag = np.abs(np.diag(r))
    rank = np.count_nonzero(diag > diag.max()*1e-8) if len(diag) else 0
    return q[:, :rank]


class ResidualAudit:
    def __init__(self, data, instrument, n_sources):
        (epochs, self.y, error), _ = prepare(data)
        self.t = epochs-epochs[0]
        self.span = float(np.ptp(epochs))
        w = (error.min()/error)**2
        self.w = w/w.sum()
        self.sw = np.sqrt(self.w)
        self.wy = self.y*self.sw
        self.q0 = self.sw[:,None]
        self.trials = survey_trials(instrument["frequencies"],n_sources)
        if self.trials != instrument["survey_trial_factor"]:
            raise ValueError("Survey multiplicity differs from the frozen search")
        self.noise = Continuum(self.t,self.y,error,self.w)
        if not np.allclose(self.noise.timescales, instrument["red_timescale_days"]):
            raise ValueError("Noise basis differs from the frozen selection")
        self.noise.variances = np.array(instrument["red_variances"])
        self.noise.configure()

    def project(self, x, q):
        xw = self.sw[:, None]*x
        return xw-q@(q.T@xw)

    def measure(self, f, q):
        v = self.project(waves(self.t,f),q)
        gram = v.T@v
        if np.linalg.eigvalsh(gram).min() < 1e-8:
            return {"status":"unresolved against nuisance fit"}
        ginv = linalg.pinvh(gram,check_finite=False)
        z = v.T@self.wy
        beta = ginv@z
        cov = self.noise.bilinear(self.sw[:,None]*v)
        logp = min(0., float(-.5*(z@linalg.pinvh(cov,check_finite=False)@z)))
        beta_cov = ginv@cov@ginv
        amp = float(np.linalg.norm(beta))
        amp_se = float(np.sqrt(max(0.,beta@beta_cov@beta))/amp) if amp else None
        return {"status":"measured", "frequency_per_day":float(f), "period_days":float(1/f),
                "semiamplitude":amp, "semiamplitude_se":amp_se, "single_grid_log_p":logp,
                "survey_model_bound":float(min(1,np.exp(max(-745., logp+np.log(self.trials))))),
                "power":float(z@ginv@z), "phase_cycles":float(np.arctan2(beta[1],beta[0])/(2*np.pi)%1)}

    def improvement(self, f, q, harmonics):
        v = self.project(waves(self.t, f, harmonics),q)
        z = v.T@self.wy
        return float(z@linalg.pinvh(v.T@v, rtol=1e-10,check_finite=False)@z)

    def profile(self, f, q):
        # First locate the correct lobe; a six-harmonic profile can have multiple extrema.
        grid = np.linspace(max(1e-10, f-1/self.span), f+1/self.span, 31)
        k = int(np.argmax([self.improvement(g,q,HARMONICS) for g in grid]))
        lo, hi = grid[max(0,k-1)], grid[min(len(grid)-1,k+1)]
        ans = optimize.minimize_scalar(lambda g: -self.improvement(g,q,HARMONICS),
                bounds=(lo,hi),method="bounded",options={"xatol":1e-5/self.span})
        return float(ans.x)

    def aliases(self, target, parent):
        report = []
        for harmonic in range(1,HARMONICS+1):
            g = parent*harmonic
            diffs = (abs(target-g),target+g)
            powers = [float(abs(self.w@np.exp(2j*np.pi*x*self.t))**2) for x in diffs]
            if max(powers) >= .0025 or diffs[0]*self.span <= 5:
                report.append({"harmonic":harmonic,"frequency_difference_rayleigh":diffs[0]*self.span,
                    "difference_window_power":powers[0],"sum_window_power":powers[1]})
        return report

    def run(self, candidate, instrument, fgrid):
        target = 1/candidate["period_days"]
        original = self.measure(fgrid,self.q0)
        exact = self.measure(target,self.q0)
        neighborhoods, parents = [], []
        for row in instrument["peaks"]:
            if row["index"] == candidate["index"]:
                continue
            f = 1/row["period_days"]
            power = self.improvement(f,self.q0,1)
            if power <= original.get("power",0) and row["score"] <= candidate["score"]:
                continue
            if abs(target-f)*self.span <= 1.25:
                # Adjacent bins of one mode represent the same frequency uncertainty.
                if abs(target-f)*self.span > .5:
                    neighborhoods.append({"period_days":row["period_days"],"score":row["score"]})
                continue
            parents.append((power,f,row))
        parents.sort(key=lambda item:-item[0])
        q = self.q0
        model = [np.ones((len(self.t),1))]
        report, potential, cap = [], [], False
        for power,f,row in parents:
            # A parent already absorbed by previously fitted clocks cannot generate a new family.
            measured = self.measure(f,q)
            if measured["status"] != "measured" or measured["survey_model_bound"] > ALPHA:
                continue
            if len(report) == MAX_FAMILIES:
                cap = True
                break
            # Estimate the nuisance frequency with the tested signal held in the model;
            # this prevents its curvature from pulling the nuisance solution toward the target.
            with_target = column_basis(np.column_stack((q, self.sw[:,None]*waves(self.t,target))))
            refined = self.profile(f, with_target)
            x = waves(self.t,refined,HARMONICS)
            model.append(x)
            q = column_basis(self.sw[:,None]*np.column_stack(model))
            aliases = self.aliases(target,refined)
            relatives = [h for h in range(2,HARMONICS+1) if abs(h*refined-target)*self.span <= 1.25]
            subharmonics = [h for h in range(2,HARMONICS+1) if abs(refined-h*target)*self.span <= 1.25]
            if relatives or subharmonics:
                potential.append({"parent_period_days":float(1/refined),"target_is_harmonic":relatives,
                                  "target_is_subharmonic":subharmonics})
            report.append({"grid_period_days":row["period_days"], "profiled_period_days":float(1/refined),
                           "initial_explained_variance":power,"parent_before_addition":measured,
                           "original_window_flag":row["window_reason"],"native_window_offsets":aliases})
        after = self.measure(fgrid,q)
        refined_after = self.measure(target,q)
        # At the original candidate frequency project the fitted nuisance waveform, including
        # its actual phases/amplitudes and the actual time mask, into a candidate sinusoid.
        y_nuisance_weighted = q@(q.T@self.wy)
        target_basis = self.project(waves(self.t,target),self.q0)
        leakage = linalg.pinvh(target_basis.T@target_basis)@(target_basis.T@y_nuisance_weighted)
        flags = []
        if cap:
            flags.append("additional significant stronger families remain beyond the audit limit")
        if neighborhoods:
            flags.append("stronger neighboring maximum unresolved at this baseline")
        if any(x["target_is_harmonic"] for x in potential):
            flags.append("coincident with a harmonic of a stronger source clock")
        if after.get("survey_model_bound",1.) > ALPHA:
            flags.append("frozen grid frequency fails after conditional nuisance projection")
        if any(x["target_is_subharmonic"] for x in potential):
            flags.append("possible common period family, physical fundamental requires review")
        if exact.get("semiamplitude",0) and np.linalg.norm(leakage) >= .5*exact["semiamplitude"]:
            flags.append("stronger fitted waveform supplies at least half the observed candidate amplitude")
        return {"status":"requires interpretation" if flags else "retained at fixed grid frequency",
                "flags":flags,"grid_before":original,"grid_after":after,
                "refined_before":exact,"refined_after":refined_after,
                "predicted_nuisance_alias_semiamplitude":float(np.linalg.norm(leakage)),
                "unresolved_neighbors":neighborhoods,"frequency_family_matches":potential,
                "nuisance_families":report,"nuisance_design_rank":q.shape[1]}


def surviving_families(records):
    groups = []
    for sid, record in records.items():
        rows = record.get("candidates",[])
        if not rows:
            continue
        span = min(record["instruments"][band]["span_days"] for band in ("bat","maxi"))
        parent = list(range(len(rows)))
        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i
        edges = []
        for i, row in enumerate(rows):
            a = 1/row["period_days"]
            for j in range(i):
                b = 1/rows[j]["period_days"]
                lo, hi = sorted([a,b])
                h = round(hi/lo)
                if 1 <= h <= HARMONICS and abs(hi-h*lo) <= (h+1)*.5/span:
                    edges.append({"periods_days":[1/a,1/b],"multiple":h})
                    parent[root(j)] = root(i)
        for k in sorted(set(root(i) for i in range(len(rows)))):
            periods = [row["period_days"] for i,row in enumerate(rows) if root(i)==k]
            groups.append({"source_id":sid, "periods_days":periods,
                           "relations":[e for e in edges if any(np.isclose(e["periods_days"][0],p,rtol=1e-12) for p in periods)]})
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path,required=True)
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--output", type=Path,required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    records = {r["roster_id"]:r for path in sorted(args.evidence.glob("search-*.jsonl.gz")) for r in read_rows(path)}
    inv = json.loads((args.evidence/"inventory.json").read_text())
    n_sources = len(inv["sources"])
    targets = {s for s,r in records.items() if r.get("candidates")}
    wanted = {f"{s}.{b}" for s in targets for b in ("bat","maxi")}
    archives = {f.name:f for d in (args.snapshot,args.evidence) if d for f in d.glob("curves-*.tar.gz")}
    curves = load_archives(sorted(archives.values()), wanted)
    with (args.output/"residuals.jsonl").open("w") as out:
        for sid in sorted(targets):
            record = records[sid]
            audits = {}
            for band in ("bat","maxi"):
                key = f"{sid}.{band}"
                if key in curves:
                    audits[band] = ResidualAudit(curves[key],record["instruments"][band],n_sources)
            for c in record["candidates"]:
                result = {"source_id":sid,"period_days":c["period_days"],"instruments":{}}
                for band in ("bat","maxi"):
                    if band not in audits:
                        result["instruments"][band] = {"status":"unresolved missing original curve"}
                        continue
                    period = c["stability"][band+"_grid_period_days"]
                    result["instruments"][band] = audits[band].run(c[band],record["instruments"][band],1/period)
                result["both_retained"] = all(r["status"]=="retained at fixed grid frequency"
                                              for r in result["instruments"].values())
                out.write(json.dumps(result,allow_nan=False)+"\n")
    (args.output/"residual-families.json").write_text(json.dumps({
        "scope":"Numerical BAT/MAXI survivors only. A family can represent one physical clock at most unless independent fundamentals are demonstrated. Diagnostic bounds condition on the selected nuisance frequencies and fitted noise spectrum.",
        "max_nuisance_families":MAX_FAMILIES,"harmonics_per_family":HARMONICS,"families":surviving_families(records)
    },indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
