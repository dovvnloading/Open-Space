"""Export complete continuous-fit records, same-window comparisons and injection trials."""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path
import numpy as np

from consolidate import HERE, read_results, write_gz

TABLES = HERE / 'tables'


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as file:
        writer = csv.DictWriter(file, keys, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def compact(data):
    return json.dumps(data, allow_nan=False, sort_keys=True, separators=(',', ':'))


def peak_summary(name, record, expanded):
    models = record['models']
    # Select across the unrestricted fit, each constrained restart and the
    # exact point-source boundary. This prevents a local minimum from producing
    # spurious evidence against the restricted hypothesis.
    fits = [(m['chi2'], m.get('te_d', m.get('te_or_rise_d')), m.get('rho', 0.), k)
            for k, m in models.items() if k in ('ps','fs') or k.startswith('profile_')]
    chi, te, rho, best = min(fits)
    bs = min([models['bin']['chi2']] + ([expanded['models']['bin']['chi2']] if expanded else []))
    flare = min([models['exp']['chi2']] + ([expanded['models']['exp']['chi2']] if expanded else []))
    deltas = {}
    for key, condition in [('long',lambda t,r:t>=.5-1e-8),('small',lambda t,r:r<=1+1e-8),
                           ('short',lambda t,r:t<=.5+1e-8),('large',lambda t,r:r>=1-1e-8)]:
        allowed = [c for c,t,r,_ in fits if condition(t,r)]
        deltas['cost_'+key+'_chi2'] = (round(min(allowed)-chi,3) if allowed
                                      and 'profile_'+key in models else '')
    return dict(event=name, window_halfwidth_d=record.get('window_halfwidth_d',4.),
                bin_width_d=record.get('bin_width_d',0.), peak_samples=record['n_peak_window'],
                calibration_cap=record['cap_factor'], best_single_start=best,
                best_single_te_d=round(te,6), best_single_rho=round(rho,6),
                best_single_chi2=chi, ps_minus_single_chi2=round(models['ps']['chi2']-chi,3),
                binary_minus_single_chi2=round(bs-chi,3),
                flare_minus_single_chi2=round(flare-chi,3),
                expanded_alternatives=bool(expanded), **deltas)


def summaries(stage, injections=None):
    full=[]
    for folder,label in [('continuous','capped KMT grid seed'),('crossings','free grid seed, capped KMT fit')]:
        for n,obj in sorted(read_results(stage/folder).items()):
            row=dict(event=n,seed_pass=label,**{k:v for k,v in obj.items() if k!='event'})
            if 'source_fractions' in row:row['source_fractions']=compact(row['source_fractions'])
            full.append(row)
    write_gz(TABLES/'continuous-fits.csv.gz',full)
    manual=read_results(HERE/'peaks')
    candidates={**read_results(stage/'comparisons'),**manual}
    expanded=read_results(stage/'expanded')
    packed={}; rows=[]
    for n,r in sorted(candidates.items()):
        if 'models' not in r:
            rows.append({'event':n,'error':r.get('error','missing model')})
            continue
        packed[n]={'initial':r,'expanded':expanded.get(n)}
        rows.append(peak_summary(n,r,expanded.get(n)))
    (TABLES/'peak-comparisons.json.gz').write_bytes(gzip.compress(compact(packed).encode(),mtime=0))
    write_csv(TABLES/'peak-comparisons.csv',rows)
    if injections:
        with injections.open() as file: trials=list(csv.DictReader(file))
        trials.sort(key=lambda r:(r['event'],float(r['t0']),float(r['te_d']),float(r['rho'])))
        groups=collections.defaultdict(list)
        for r in trials:
            survey='Gaia' if r['event'].startswith('Gaia') else r['event'].split('-')[0]
            groups[(survey,r['te_d'],r['rho'])].append(r)
        stats=[]
        for (s,t,r),data in sorted(groups.items()):
            total=len(data); sampled=sum(int(x['two_each_side']) for x in data)
            recovered=sum(int(x['paired_pass']) for x in data)
            per_alert=collections.defaultdict(list)
            for trial in data:per_alert[trial['event']].append(int(trial['paired_pass']))
            rates=np.array([np.mean(per_alert[n]) for n in sorted(per_alert)])
            assert len({len(v) for v in per_alert.values()})==1
            # Preserve the 12 windows as a single block for each resampled
            # sightline; shared windows in different geometry rows stay paired.
            rng=np.random.default_rng(64823)
            draws=rates[rng.integers(0,len(rates),size=(5000,len(rates)))].mean(axis=1)
            lo,hi=np.quantile(draws,[.025,.975])
            stats.append(dict(survey=s,te_d=t,rho=r,alerts=len({x['event'] for x in data}),
                              trials=total,two_each_side=sampled,recovered=recovered,
                              controls_passing=sum(int(x['control_pass']) for x in data),
                              fraction_all=round(recovered/total,5),
                              sightline_bootstrap_low=round(lo,5) if rates.max()>0 else '',
                              sightline_bootstrap_high=round(hi,5) if rates.max()>0 else '',
                              fraction_when_sampled=round(recovered/sampled,5) if sampled else ''))
        write_gz(TABLES/'injections.csv.gz',trials)
        write_csv(TABLES/'injection-summary.csv',stats)
    return len(full),len(rows),len(expanded)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path.cwd())
    parser.add_argument('--injections',type=Path)
    args=parser.parse_args()
    print(summaries(args.cache/'staging',args.injections))
