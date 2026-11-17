"""Pair DES and Euclid measurements on identical historic-orbit indices."""
import csv,gzip,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def rows(p):
    if not p.exists() and Path(str(p)+'.gz').exists():p=Path(str(p)+'.gz')
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt') as f:return list(csv.DictReader(f))

def main():
    des=HERE.parent/'archives';euclid=HERE.parent/'orbits'/'euclid'
    e={int(x['grid_id']):x for x in rows(euclid/'fixed_grid_trials.csv')};d={int(x['track']):x for x in rows(des/'difference_tracks.csv.gz')}
    disp=rows(des/'difference_injections.csv.gz');dr=rows(des/'baseline_draws.csv');selected={int(x['track']) for x in dr};assert len(selected)==320
    assert e.keys()==d.keys()==set(range(3617))
    u={k for k in e if int(e[k]['eligible_after_DES_mask'])};v={k for k in d if int(d[k]['valid'])};assert v&selected=={int(x['track']) for x in dr if int(x['valid'])}
    assert not any(int(x['flagged']) for x in d.values()),'Positive DES grid sources need adjudication before interpreting a union.'
    with (HERE/'matched_recovery.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['DES_r_AB','Euclid_VIS_AB','draws','DES_recovered','Euclid_recovered','both','either','neither'])
        for rmag in [22.,23.,23.5,24.]:
            drec={int(x['track']) for x in disp if float(x['r_AB'])==rmag and int(x['recovered'])}
            for vismag in [24.,25.,26.,26.5,27.]:
                erec={k for k in selected & u if int(e[k][f'recovered_VIS{vismag:g}'])};w.writerow([rmag,vismag,len(selected),len(drec),len(erec),len(drec&erec),len(drec|erec),len(selected-drec-erec)])
    (HERE/'matched_geometry.json').write_text(json.dumps(dict(family=len(e),DES=len(v),Euclid=len(u),both=len(v&u),either=len(v|u),neither=len(e)-len(v|u)),indent=2)+'\n')
if __name__=='__main__':main()
