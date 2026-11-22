"""Post-selection checks on known clocks and neighbouring monitored positions.

Run only in Actions. These diagnostics do not participate in blind selection.
"""

import argparse
import gzip
import io
import json
import math
import tarfile
from pathlib import Path

import numpy as np

from search import Continuum, bounded_uncertainties


# Times are rounded literature periods for positive controls. The wider Her X-1
# tolerance accommodates changes in the duration of its precession cycle.
ECLIPSING = ("Falanga et al. 2015, A&A 577, A130",
             "https://ui.adsabs.harvard.edu/abs/2015A%26A...577A.130F/abstract")
HER = ("Staubert et al. 2009, A&A 494, 1025",
       "https://ui.adsabs.harvard.edu/abs/2009A%26A...494.1025S/abstract")
CONTROLS = [
    ("Vela X-1", "orbital", 8.96436, .0001, ECLIPSING),
    ("Cen X-3", "orbital", 2.08714, .0001, ECLIPSING),
    ("SMC X-1", "orbital", 3.89192, .0001, ECLIPSING),
    ("H 1538-522", "orbital", 3.7284, .0001, ECLIPSING),
    ("OAO 1657-415", "orbital", 10.4473, .0001, ECLIPSING),
    ("Her X-1", "superorbital", 34.85, .025, HER),
]
REFERENCE_MJD = 55000.
RADIUS_DEG = 2.


def read_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def key(value):
    return "".join(c.lower() for c in value.replace("−","-") if c.isalnum() or c in "+-")


def relevant(period, trial, span, frac):
    return period > 0 and trial > 0 and abs(1/period-1/trial) <= 1/max(span,1.)+frac/trial


def control_band(instrument, period, tolerance):
    status = instrument.get("status", "missing")
    if status != "searched":
        return dict(stage=status, matching_peaks=[])
    if not instrument["min_period_days"] <= period <= instrument["max_period_days"]:
        return dict(stage="outside search range", matching_peaks=[])
    span = instrument["span_days"]
    hits = [r for r in instrument["peaks"] if relevant(r["period_days"],period,span,tolerance)]
    if not hits:
        return dict(stage="below blind Monte Carlo threshold", matching_peaks=[])
    good = [r for r in hits if r.get("window_reason") is None and r.get("model_ok") and r["survey_fap_bound"]<=.05]
    if good:
        stage = "passed instrument gates"
    else:
        reasons = set()
        for r in hits:
            if r.get("window_reason"):
                reasons.add("window/alias")
            if not r.get("model_ok"):
                reasons.add("continuum")
            if r["survey_fap_bound"]>.05:
                reasons.add("survey FAP")
        stage = "failed "+", ".join(sorted(reasons))
    return dict(stage=stage, matching_peaks=hits)


def controls(roster, records, destination):
    results = []
    for name, cycle, period, tol, (reference, url) in CONTROLS:
        c = dict(source=name, known_cycle=cycle, reference_period_days=period, reference=reference,
                 url=url, relative_tolerance_added_to_rayleigh=tol)
        matched = [b for b in roster if key(name) in {key(a) for a in [b["id"],b["name"],*b["aliases"]]}]
        if len(matched) != 1 or matched[0]["id"] not in records:
            c.update(stage="missing or ambiguous roster/search record", roster_matches=[b["id"] for b in matched])
        else:
            ident = matched[0]["id"]
            r = records[ident]
            instruments = r.get("instruments", {})
            bands = {k:control_band(instruments.get(k,{}),period,tol) for k in ("bat","maxi")}
            def near_band(band, p):
                return relevant(p,period,instruments[band]["span_days"],tol)
            pairs = [p for p in r.get("pairs",[]) if near_band("bat",p["bat_grid_period_days"])
                     and near_band("maxi",p["maxi_grid_period_days"])]
            final = [p["period_days"] for p in r.get("candidates",[]) if near_band("bat",p["period_days"])]
            if final:
                stage = "blind recovery through halves"
            elif pairs:
                stage = "matched; failed stability"
            elif all(v["stage"]=="passed instrument gates" for v in bands.values()):
                stage = "independent peaks did not agree"
            else:
                stage = "did not reach independent match"
            c.update(stage=stage, roster_id=ident, source_status=r["status"], bands=bands,
                     matches=pairs, surviving_periods_days=final)
        results.append(c)
    dump(destination/"controls.json", results)
    lines = ["# Published positive controls", "",
             "These source identities and reference periods entered after blind peak selection. An orbital control matches within one instrumental Rayleigh width plus 0.01% of the reference frequency. Her X-1 allows an extra 2.5% for its variable superorbital cycle. Only a matched fundamental counts as recovery.", "",
             "| Source | Published period (days) | BAT | MAXI | Joint outcome | Reference |",
             "| --- | ---: | --- | --- | --- | --- |"]
    for c in results:
        bands = c.get("bands",{})
        lines.append(f"| {c['source']} | {c['reference_period_days']} | {bands.get('bat',{}).get('stage','unavailable')} | "
                     f"{bands.get('maxi',{}).get('stage','unavailable')} | {c['stage']} | [{c['reference']}]({c['url']}) |")
    (destination/"controls.md").write_text("\n".join(lines)+"\n")


def separation(pos, other):
    ra,dec,r,d = np.deg2rad([*pos,*other])
    a = math.sin((d-dec)/2)**2+math.cos(dec)*math.cos(d)*math.sin((r-ra)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1.,max(a,0.)))))


def field_members(survivor, by_bat, by_maxi, inverse):
    row = by_bat[survivor["source_id"]]
    entries = []
    for band, catalogue, target_id in (("bat",by_bat,row["id"]),("maxi",by_maxi,row["maxi_id"])):
        target = catalogue.get(target_id, {})
        pos = target.get("position")
        if pos is None:
            entries.append(dict(band=band, status="unresolved target position", neighbors=[]))
            continue
        neighbors = []
        for ident, obj in catalogue.items():
            if ident==target_id or obj.get("position") is None:
                continue
            dist = separation(pos,obj["position"])
            if dist<=RADIUS_DEG:
                holders = [ident] if band=="bat" else inverse.get(ident, [])
                neighbors.append(dict(id=ident, name=obj.get("name") or obj.get("aliases",[ident])[0],
                                      position=obj["position"], separation_deg=dist,
                                      curve_key=f"{holders[0]}.{band}" if len(holders)==1 else None))
        entries.append(dict(band=band, target=target_id, position=pos, target_curve=f"{row['id']}.{band}",
                            neighbors=sorted(neighbors,key=lambda v:v["separation_deg"])))
    return entries


def load_archives(archives, wanted):
    curves = {}
    for archive in archives:
        if not wanted-curves.keys():
            break
        with tarfile.open(archive,"r|*") as tar:
            for member in tar:
                path = Path(member.name)
                if path.parent!=Path("curves") or path.suffix!=".npz" or path.stem not in wanted or not member.isfile():
                    continue
                stream = tar.extractfile(member)
                if stream is None:
                    continue
                with np.load(io.BytesIO(stream.read()),allow_pickle=False) as data:
                    curves[path.stem] = tuple(data[k] for k in ("t","y","dy"))
    return curves


def fixed_fit(curve, period, start, stop):
    t,y,e = curve
    m = (t>=start)&(t<=stop)&np.isfinite(t)&np.isfinite(y)&np.isfinite(e)&(e>0)
    t,y,e = t[m],y[m],e[m]
    if len(t)<40 or np.ptp(t)<3*period:
        return {"status":"insufficient common coverage","rows":len(t)}
    e, weight_rules = bounded_uncertainties(t,e)
    w = (e.min()/e)**2
    w /= w.sum()
    if 1/(w@w) < 20:
        return {"status":"insufficient effective sampling","rows":len(t)}
    theta = 2*np.pi*(t-REFERENCE_MJD)/period
    x = np.column_stack((np.ones(len(t)),np.cos(theta),np.sin(theta)))
    v = w[:,None]*x
    normal = x.T@v
    if np.linalg.cond(normal)>1e8:
        return {"status":"ill-conditioned common fit","rows":len(t)}
    h = np.linalg.inv(normal)
    beta = h@(v.T@y)
    model = Continuum(t,y,e,w)
    if not model.enough:
        return {"status":"insufficient covariance lags","rows":len(t)}
    cov = (h@model.bilinear(v)@h)[1:,1:]
    amp = float(np.hypot(*beta[1:]))
    p = math.exp(-.5*max(0.,float(beta[1:]@np.linalg.pinv(cov)@beta[1:])))
    return dict(status="fitted", rows=len(t), weight_rules=weight_rules, start_mjd=float(t[0]), stop_mjd=float(t[-1]),
                semiamplitude=amp, center=float(beta[0]), phase_cycles=float(np.arctan2(beta[2],beta[1])/(2*np.pi)%1),
                conservative_95_upper_semiamplitude=amp+math.sqrt(5.992*max(0.,float(np.linalg.eigvalsh(cov).max()))),
                frequency_local_red_model_p=p)


def inspect_field(candidate, fields, curves):
    verdict = set()
    for field in fields:
        if field.get("status"):
            verdict.add(field["status"])
            continue
        target = curves.get(field["target_curve"])
        if target is None:
            field["status"] = "unresolved missing target curve"
            verdict.add(field["status"])
            continue
        field["target_full"] = fixed_fit(target,candidate["period_days"],min(target[0]),max(target[0]))
        if field["target_full"]["status"]!="fitted":
            verdict.add("unresolved target fit")
        for nbr in field["neighbors"]:
            nearby = curves.get(nbr["curve_key"])
            if nearby is None:
                nbr["status"] = "unresolved missing neighbor curve"
            else:
                tmin,tmax = max(min(target[0]),min(nearby[0])),min(max(target[0]),max(nearby[0]))
                a = fixed_fit(target,candidate["period_days"],tmin,tmax)
                b = fixed_fit(nearby,candidate["period_days"],tmin,tmax)
                nbr.update(target_on_common_span=a, neighbor_on_common_span=b)
                if a["status"]!="fitted" or b["status"]!="fitted" or a["semiamplitude"]<=0:
                    nbr["status"] = "unresolved common fit"
                else:
                    ratio = b["semiamplitude"]/a["semiamplitude"]
                    phase = (b["phase_cycles"]-a["phase_cycles"]+.5)%1-.5
                    nbr.update(amplitude_ratio_neighbor_to_target=ratio, phase_offset_cycles=phase)
                    if a["frequency_local_red_model_p"] >= .01:
                        nbr["status"] = "target modulation weak on overlapping span"
                    elif b["frequency_local_red_model_p"] < .01:
                        nbr["status"] = "neighbor also modulates at this frequency; origin unresolved"
                    elif b["conservative_95_upper_semiamplitude"]>=a["semiamplitude"]:
                        nbr["status"] = "neighbor's fitted upper amplitude admits the target amplitude"
                    else:
                        nbr["status"] = "neighbor has a lower amplitude upper estimate"
            if nbr["status"]!="neighbor has a lower amplitude upper estimate":
                verdict.add(nbr["status"])
    return dict(source=candidate["source"], source_id=candidate["source_id"], period_days=candidate["period_days"],
                radius_deg=RADIUS_DEG, phase_zero_mjd=REFERENCE_MJD, fields=fields,
                status="requires spatial review" if verdict else "no comparable measured signal in roster neighbors",
                flags=sorted(verdict))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    root = args.evidence
    inventory = json.loads((root/"inventory.json").read_text())
    records = {r["roster_id"]:r for p in sorted(root.glob("search-*.jsonl.gz")) for r in read_rows(p)}
    controls(inventory["sources"],records,args.output)
    candidates = list(read_rows(root/"survivors.jsonl"))
    bat = {b["id"]:b for b in inventory["sources"]}
    maxi = {m["id"]:m for m in inventory["maxi_sources"]}
    inverse = {}
    for b in bat.values():
        if b["maxi_id"]:
            inverse.setdefault(b["maxi_id"],[]).append(b["id"])
    fields = [field_members(c,bat,maxi,inverse) for c in candidates]
    wanted = {f.get("target_curve") for fs in fields for f in fs}
    wanted.update(n.get("curve_key") for fs in fields for f in fs for n in f["neighbors"])
    wanted.discard(None)
    archives = {p.name:p for d in [args.snapshot,root] if d for p in d.glob("curves-*.tar.gz")}
    curves = load_archives(sorted(archives.values()),wanted)
    with (args.output/"confusion.jsonl").open("w") as handle:
        for c, fs in zip(candidates,fields):
            handle.write(json.dumps(inspect_field(c,fs,curves),allow_nan=False)+"\n")
    dump(args.output/"confusion_scope.json",dict(candidates=len(candidates),radius_deg=RADIUS_DEG,
         epoch_mjd=REFERENCE_MJD, instrument_references={
             "BAT":"https://ui.adsabs.harvard.edu/abs/2013ApJS..209...14K/abstract",
             "MAXI":"https://ui.adsabs.harvard.edu/abs/2009PASJ...61..999M/abstract"},
         fixed_frequency="1 / refined BAT period", method="Inverse-variance floating-mean sine and cosine in the common time span. Errors use the source's fitted OU-mixture covariance. Each conditional p-value describes one comparison and cannot establish a discovery.",
         limits="Both public monitor rosters are screened independently. The two-degree radius is an angular neighborhood, not a calibrated image response or a source localization. Missing neighbor products, crowding, off-axis artifacts and unlisted transients can prevent attribution."))


if __name__=="__main__":
    main()
