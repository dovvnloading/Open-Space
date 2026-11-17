"""Count unique sources and independent nights across overlapping mass screens and surveys."""
import csv
import gzip
import json
import zipfile
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))


def main():
    masses=defaultdict(dict)
    samples=[('gaia',ROOT),('extension',ROOT/'extension')]
    prior=defaultdict(set)
    exposures=defaultdict(lambda:defaultdict(dict))
    raw_mrs=defaultdict(set)
    diag={}
    for sample,d in samples:
        for r in rows(ROOT/sample/'mass_screen.csv'):
            sid=r['source_id']
            masses[sid][sample]=max(masses[sid].get(sample,0.), float(r['m2_screen_msun']))
        for r in rows(ROOT/sample/'source_diagnostics.csv'):
            diag[r['source_id']]=r
        for r in rows(d/'novelty/prior_matches.csv'):
            prior[r['source_id']].add(r['study'])
        for survey in ('desi','sdss'):
            for r in rows(d/'epochs'/(survey+'_visits.csv')):
                k=(r['targetid'],r['expid']) if survey=='desi' else (r['instrument'],r['spectrum_pk'])
                name=survey if survey=='desi' else r['instrument'].split('-')[0].lower()
                entry=(int(float(r['time_jd'])-2400000.5),r['quality']=='1', r['instrument'])
                exposures[r['source_id']][name][k]=entry
        for r in rows(d/'lamost/exposures.csv'):
            exposures[r['source_id']]['lamost'][(r['obsid'],r['lmjm'])]=(r['obsdate'],True,'LAMOST MRS')
        for receipt in json.loads((d/'lamost/raw/receipts.json').read_text()):
            f=d/'lamost/raw'/receipt['file']
            if f.exists(): data=f.read_bytes()
            else:
                with zipfile.ZipFile(d/'lamost/raw/responses.zip') as zipf: data=zipf.read(receipt['file'])
            for r in csv.DictReader(gzip.decompress(data).decode('utf-8-sig').splitlines()):
                raw_mrs[r['gaia_source_id']].add(r['mobsid'])
    clusters={r['source_id']:r['cluster_name'] for r in rows(ROOT/'environment/cluster_matches.csv')}
    stats={}
    output=[]
    surveys=sorted({s for instruments in exposures.values() for s in instruments})
    for sid in sorted(masses,key=int):
        m=masses[sid]
        r=dict(source_id=sid,screen_msun=round(max(m.values()),5),initial_msun=round(m['gaia'],5) if 'gaia' in m else '',
               flame_msun=round(m['extension'],5) if 'extension' in m else '',published_study_matches=len(prior[sid]),
               cluster_catalog=clusters.get(sid,''),gaia_rv_amplitude=diag[sid].get('rv_amplitude_robust',''),
               lamost_raw_arm_products=len(raw_mrs[sid]))
        for inst in surveys:
            v=exposures[sid][inst].values()
            r[inst+'_visits']=len(v)
            r[inst+'_clean']=sum(ok for _,ok,_ in v)
            r[inst+'_clean_nights']=len({t for t,ok,_ in v if ok})
            by_inst=defaultdict(set)
            for t,ok,site in v:
                if ok: by_inst[site].add(t)
            r[inst+'_max_nights_one_instrument']=max(map(len,by_inst.values()),default=0)
        output.append(r)
    for label,lower in [('above_1p4',1.4),('above_3',3.)]:
        group=[r for r in output if max(masses[r['source_id']].values())>lower]
        stats[label]=dict(sources=len(group),matched_any_prior_table=sum(bool(r['published_study_matches']) for r in group),
                          matched_cluster_catalog=sum(bool(r['cluster_catalog']) for r in group),
                          with_gaia_robust_amplitude=sum(r['gaia_rv_amplitude'] not in ('','null') for r in group),
                          any_clean_survey=sum(any(r[s+'_clean'] for s in surveys) for r in group),
                          with_clean_repeats_same_instrument=sum(any(r[s+'_max_nights_one_instrument']>=2 for s in surveys) for r in group),
                          with_three_clean_nights_same_instrument=sum(any(r[s+'_max_nights_one_instrument']>=3 for s in surveys) for r in group))
        for inst in surveys:
            stats[label][inst]=dict(raw_matched_sources=sum(bool(r['lamost_raw_arm_products'] if inst=='lamost' else r[inst+'_visits']) for r in group),
                clean_sources=sum(bool(r[inst+'_clean']) for r in group),
                clean_epochs=sum(r[inst+'_clean'] for r in group),
                clean_multinight=sum(r[inst+'_clean_nights']>=2 for r in group))
    with gzip.GzipFile(filename='',mode='wb',fileobj=(ROOT/'coverage.csv.gz').open('wb'),mtime=0) as out:
        import io
        with io.TextIOWrapper(out,encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,list(output[0]),lineterminator='\n');writer.writeheader();writer.writerows(output)
    (ROOT/'coverage.json').write_text(json.dumps(stats,indent=2)+'\n')


if __name__=='__main__':main()
