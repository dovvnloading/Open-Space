"""Single-visit ephemeris-directed positive control with same native masks and stationary subtraction."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from link_visits import make_cube, evaluate, deep_score, search


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('control',type=Path)
    p.add_argument('ee',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args()
    c=make_cube(a.control,a.ee,True)
    target=int(c['unmasked_control_body'])
    idx=np.flatnonzero(c['known_id']==target).item()
    path=c['known_path'][idx]
    # Include a local astrometric-offset map and a comparison sample well away from the ephemeris.
    offsets=np.array([(x,y) for y in range(-3,4) for x in range(-3,4)]+
                     [(r*np.cos(t),r*np.sin(t)) for r in [9,13,17,21] for t in np.arange(24)*np.pi/12])
    px=path[None,:,:,0]+offsets[:,0,None,None]
    py=path[None,:,:,1]+offsets[:,1,None,None]
    bank=dict(x_path=px,y_path=py,x=px[:,:,px.shape[-1]//2])
    cuts=search.Cuts(psf_support=.8)
    kernels=search.validate_cube(c,min_visits=1)
    num,info=evaluate(c,bank,kernels,search.prepare(c,kernels,cuts),cuts)
    rank,score,early,late,eligible=deep_score(num,info)
    valid_local=np.flatnonzero(eligible[:49]); valid_off=np.flatnonzero(eligible[49:])+49
    peak=int(valid_local[np.argmax(score[valid_local])]) if len(valid_local) else 24
    def stats(k):
        return dict(dx=offsets[k,0],dy=offsets[k,1],score=score[k],early=early[k],late=late[k],
                    n_exp=int((info[k]>0).sum()),flux_Jy=(float(num[k].sum()/info[k].sum()) if info[k].sum()>0 else None))
    report=dict(body=target,cube_sha256=search.file_hash(a.control),center=stats(24),local_peak=stats(peak),
                eligible_controls=len(valid_off),
                displaced_max=(float(score[valid_off].max()) if len(valid_off) else None),
                displaced_95_percentile=(float(np.percentile(score[valid_off],95)) if len(valid_off) else None),
                percentile_comparison='Offset diagnostics share pixels; not a false-alarm probability.')
    a.out.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
