"""Audit published repeat identifications and PN–cluster projections in Fornax3D."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
ARCSEC = 180 * 3600 / np.pi


def load(name, expected):
    p = HERE / name
    if hashlib.sha256(p.read_bytes()).hexdigest() != expected:
        raise ValueError(f'Archived source bytes changed: {name}')
    with p.open() as f:
        return list(csv.DictReader(f))


def sexagesimal(s):
    deg, minutes, seconds = map(float, s.split())
    return np.copysign(abs(deg) + minutes/60 + seconds/3600, deg)


def coordinates(row):
    a = row['RAJ2000']
    if ' ' in a:
        return 15*sexagesimal(a), sexagesimal(row['DEJ2000'])
    return float(a), float(row['DEJ2000'])


def unit_vectors(rows):
    a, d = np.deg2rad([coordinates(r) for r in rows]).T
    return np.column_stack((np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)))


def separation(a,b):
    return 2*ARCSEC*np.arcsin(np.minimum(1,cdist(a,b)/2))


def write(name, rows):
    with (HERE / name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys(),lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def audit(rows):
    links=[]
    gal_names=sorted({r['Galaxy'] for r in rows})
    for gal in gal_names:
        c=[r for r in rows if r['Galaxy']==gal and r['Loc']=='center']
        outside=[r for r in rows if r['Galaxy']==gal and '/' in r['index']]
        if not outside:
            continue
        cs=unit_vectors(c)
        by_index={r['index']:i for i,r in enumerate(c)}
        for field in sorted({r['Loc'] for r in outside}):
            # Use all detections, including entries without a published repeat,
            # in the reciprocal nearest-neighbor and ambiguity tests.
            outer=[r for r in rows if r['Galaxy']==gal and r['Loc']==field]
            sep=separation(unit_vectors(outer),cs)
            for i,r in enumerate(outer):
                if '/' not in r['index']:
                    continue
                own,named=map(str.strip,r['index'].split('/'))
                j=by_index[named]
                best,runner=np.argsort(sep[i])[:2]
                alt=c[best]
                dv=float(r['LOSV'])-float(alt['LOSV'])
                dm=float(r['m5007'])-float(alt['m5007'])
                far=float(sep[i,j])
                near=float(sep[i,best])
                reciprocal=bool(sep[:,best].argmin()==i)
                unique=bool(sep[i,runner]>2 and np.partition(sep[:,best],1)[1]>2)
                # Deliberately strict. No frame registration or fitted shifts.
                secure=bool(j!=best and far>3 and near<0.6 and reciprocal and unique
                            and abs(dv)<30 and abs(dm)<0.4)
                links.append(dict(galaxy=gal,outer=own,outer_id=r['PN'],named=named,
                                  named_id=c[j]['PN'],named_arcsec=round(far,4),
                                  nearest=alt['index'],nearest_id=alt['PN'],
                                  nearest_arcsec=round(near,4),reciprocal=reciprocal,
                                  isolated_2arcsec=unique,dv_named_kms=round(float(r['LOSV'])-float(c[j]['LOSV']),2),
                                  dv_nearest_kms=round(dv,2),dm_named=round(float(r['m5007'])-float(c[j]['m5007']),3),
                                  dm_nearest=round(dm,3),outer_label=r['label'],nearest_label=alt['label'],
                                  meets_relink_criteria=secure))
    return links


def projections(pn, gc):
    samples=[]
    for gal in sorted({r['Galaxy'] for r in pn} & {r['Galaxy'] for r in gc}):
        ps=[r for r in pn if r['Galaxy']==gal and r['label']=='PN' and r['index'].count('/')==0]
        gs=[r for r in gc if r['Galaxy']==gal]
        for row in ps:
            a,d=coordinates(row)
            if row['PN'][5:] != (row['RAJ2000']+row['DEJ2000']).replace(' ',''):
                continue
            samples.append((gal,row,a,d))
    results=[]
    controls=[]
    for angle in [None,*np.arange(24)*np.pi/12]:
        matches=[]
        for gal in sorted({s[0] for s in samples}):
            rows=[(r,a,d) for g,r,a,d in samples if gal==g]
            gs=[r for r in gc if r['Galaxy']==gal]
            aa,dd=np.array([(a,d) for r,a,d in rows]).T
            gg=unit_vectors(gs)
            if angle is not None:
                aa+=5*np.cos(angle)/(3600*np.cos(np.deg2rad(dd)))
                dd+=5*np.sin(angle)/3600
            vv=unit_vectors([dict(RAJ2000=str(a),DEJ2000=str(d)) for a,d in zip(aa,dd)])
            dist=separation(vv,gg)
            for i,j in zip(*np.where(dist<=1)):
                r=rows[i][0]
                dv=float(r['LOSV'])-float(gs[j]['vLOS'])
                matches.append(dict(galaxy=gal,pn=r['PN'],index=r['index'],gc=gs[j]['SimbadName'],
                                    separation_arcsec=round(dist[i,j],4),dv_kms=round(dv,2),
                                    pn_sigma_kms=r['e_LOSV'],gc_sigma_kms=gs[j]['e_vLOS']))
        count=sum(m['separation_arcsec']<=0.5 and abs(m['dv_kms'])<=30 for m in matches)
        controls.append(int(count))
        if angle is None:
            results=matches
    return results,dict(selected_pn_rows=len(samples),near_counts_unshifted=controls[0],
                        counts_at_5arcsec_24_angles=controls[1:])


def main():
    p=load('pn-observations.csv','cb531324169dbbd46d0cc0cda0a1293d7cd2dbcae71c3f7fc50f52bbf0e1b23d')
    g=load('cluster-observations.csv','d6ee0facd1809927600d16ab4b3e32106382804c506d89fad3277dd5f830a96e')
    links=audit(p)
    strong=[r for r in links if r['meets_relink_criteria']]
    write('repeat-audit.csv',links)
    matches,summary=projections(p,g)
    write('projected-pairs-within-1arcsec.csv',matches)
    summary.update(pn_table_rows=len(p),cluster_table_rows=len(g),published_repeat_links=len(links),
                   strong_alternatives=len(strong),galaxies_with_alternatives=sorted({r['galaxy'] for r in strong}),
                   median_named_sep_arcsec=float(np.median([r['named_arcsec'] for r in strong])),
                   median_alternative_sep_arcsec=float(np.median([r['nearest_arcsec'] for r in strong])),
                   name_coordinate_disagreements=[r['PN'] for r in p if
                    r['PN'][5:] != (r['RAJ2000']+r['DEJ2000']).replace(' ','')])
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    main()
