"""Replay the cluster encounter screen and member bootstrap (numpy, pandas, scipy)."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
PC_PER_KMS_MYR = 1.022712165


def phase_space(ra, dec, distance, pmra, pmdec, rv):
    """ICRS unit vectors and heliocentric position (pc) and velocity (km/s)."""
    a, b = np.radians(ra), np.radians(dec)
    radial = np.stack([np.cos(a)*np.cos(b), np.sin(a)*np.cos(b), np.sin(b)], axis=-1)
    east = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    north = np.stack([-np.cos(a)*np.sin(b), -np.sin(a)*np.sin(b), np.cos(b)], axis=-1)
    expand = lambda x: np.asarray(x)[..., None]
    pos = radial * expand(distance)
    vel = radial * expand(rv) + 0.00474047 * expand(distance) * (east * expand(pmra) + north * expand(pmdec))
    return pos, vel


def encounter(first, second):
    pos = first[0] - second[0]
    vel = first[1] - second[1]
    speed = np.linalg.norm(vel, axis=-1)
    time = -np.sum(pos * vel, axis=-1) / np.maximum(speed**2, 1e-20)
    peri = np.linalg.norm(pos + vel * time[..., None], axis=-1)
    return time / PC_PER_KMS_MYR, peri, np.linalg.norm(pos, axis=-1), speed


def screen(hunt):
    vectors = phase_space(*[hunt[c].to_numpy() for c in ['RAdeg','DEdeg','dist50','pmRA','pmDE','RV']])
    i, j = cKDTree(vectors[0]).query_pairs(420, output_type='ndarray').T
    first = [v[i] for v in vectors]
    second = [v[j] for v in vectors]
    t, d, now, speed = encounter(first, second)
    grid = np.hypot(d, (t - 2 * np.round(t/2)) * PC_PER_KMS_MYR * speed)
    radii = hunt.rJpc.to_numpy()[i] + hunt.rJpc.to_numpy()[j]
    use = (np.abs(t)<8) & (now>30) & (d<20) & (grid>20) & (d<0.8*radii)
    return pd.DataFrame(dict(first=hunt.Name.to_numpy()[i], second=hunt.Name.to_numpy()[j],
                             t_myr=t, pericentre_pc=d, current_pc=now, speed_kms=speed,
                             grid_pc=grid, jacobi_sum_pc=radii))[use].sort_values(['first','second'])


def member_draws(group, probability, floor, draws, rng):
    group = group[group.probs >= probability]
    astro = group[group.e_Plx < 0.18]
    rv = group.dropna(subset=['RV','e_RV'])
    med = rv.RV.median()
    spread = (rv.RV-med).abs().median()
    rv = rv[(rv.e_RV<12) & ((rv.RV-med).abs()<max(6,4.5*spread))]
    ix = rng.integers(len(astro), size=(draws, len(astro)))
    jx = rng.integers(len(rv), size=(draws, len(rv)))
    vals = [np.median(astro[c].to_numpy()[ix], axis=1) for c in ['RA_ICRS','DE_ICRS','Plx','pmRA','pmDE']]
    vals[2] = 1000 / (vals[2]+0.017+rng.normal(0, floor, draws))
    rvmed = np.median(rv.RV.to_numpy()[jx], axis=1)
    return phase_space(*vals, rvmed), [len(astro),len(rv)]


def normalise(name):
    return re.sub('[^a-z0-9]', '', name.lower())


def aliases(hr, ucc):
    return {normalise(s) for value in (hr['AllNames'],hr['Name'],ucc['Name(s)']) for s in re.split('[;,]', value)}


def summarise(arr):
    return np.round(np.percentile(arr,[16,50,84]),4).tolist()


def run(draws):
    hunt = pd.read_csv(INPUTS/'hunt_reffert.csv')
    members = pd.read_csv(INPUTS/'ucc_members.csv',dtype={'Source':str})
    centres = pd.read_csv(INPUTS/'ucc_centres.csv').set_index('name')
    previous_rv = pd.read_csv(INPUTS/'tarricq_rv.csv').set_index('Cluster')
    prior = pd.read_csv(INPUTS/'published_pairs.csv')
    pairs = [{normalise(a),normalise(b)} for a,b in zip(prior.name1,prior.name2)]
    other = pd.read_csv(INPUTS/'other_groups.csv')
    systems = {key:{normalise(n) for n in rows.cluster} for key,rows in other.groupby(['catalogue','system'])}
    id_systems = {key:set(group.hr_id) for key,group in other.groupby(['catalogue','system'])}
    leads = screen(hunt)
    indexed = hunt.set_index('Name',drop=False)
    results = []
    for row in leads.itertuples(index=False):
        names = [normalise(n) for n in (row.first,row.second)]
        populations = [members[members.name==n] for n in names]
        rows = [indexed.loc[n] for n in (row.first,row.second)]
        synonyms = [aliases(h,centres.loc[n]) for h,n in zip(rows,names)]
        prior_match = any(pair & synonyms[0] and pair & synonyms[1] for pair in pairs)
        catalogue_matches = [list(key) for key,group in systems.items() if group & synonyms[0] and group & synonyms[1]]
        ids = {r.ID for r in rows}
        exact_groups = [list(key) for key,group in id_systems.items() if ids <= group]
        shared = set(populations[0].Source) & set(populations[1].Source)
        rec = {k:(round(v,5) if isinstance(v,float) else v) for k,v in row._asdict().items()}
        rec.update(dict(shared_sources=len(shared),listed_in_prior=prior_match,alias_group_hits=catalogue_matches,hr_id_group_hits=exact_groups,
                        half_number_sum_pc=round(sum(r.r50pc for r in rows),4)))
        historical = [r.Name for r in rows if r.Name in previous_rv.index]
        if historical:
            states = [phase_space(r.RAdeg,r.DEdeg,r.dist50,r.pmRA,r.pmDE,
                                  previous_rv.loc[r.Name,'RVmean'] if r.Name in historical else r.RV) for r in rows]
            hist_time, hist_distance, _, _ = encounter(*states)
            rec['tarricq_rv_substitution'] = dict(replaced=historical,time_myr=round(float(hist_time),4),
                                                   pericentre_pc=round(float(hist_distance),4))
        scenarios = []
        for p,floor in [(.5,0),(.5,.01),(.8,.01),(.5,.02)]:
            # Reuse the same resamples for each systematic-floor sensitivity.
            rng = np.random.default_rng(24319)
            samples = [member_draws(g,p,floor,draws,rng) for g in populations]
            t,peri,now,speed = encounter(samples[0][0],samples[1][0])
            scenarios.append(dict(min_membership=p,parallax_floor_mas=floor,
                                  astro_rv_counts=[s[1] for s in samples],
                                  time_myr=summarise(t),pericentre_pc=summarise(peri),
                                  speed_kms=summarise(speed),current_pc=summarise(now),
                                  fraction_within_jacobi_sum=round(float(np.mean(peri<row.jacobi_sum_pc)),4),
                                  fraction_within_half_number_sum=round(float(np.mean(peri<rec['half_number_sum_pc'])),4)))
        rec['member_bootstrap']=scenarios
        results.append(rec)
    return {'hr_clusters':len(hunt),'search_pairs_before_distance_cut':len(hunt)*(len(hunt)-1)//2,
            'published_rows':len(prior),'published_epoch_range_myr':[prior.t_myr.min(),prior.t_myr.max()],
            'bootstrap_draws':draws,'candidates':results}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--draws',type=int,default=6000)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    data=json.dumps(run(args.draws),indent=2,allow_nan=False)+'\n'
    if args.output: args.output.write_text(data)
    else: print(data,end='')
