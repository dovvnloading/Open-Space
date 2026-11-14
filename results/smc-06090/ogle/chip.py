"""Match non-contact binaries on the target's chip before inspecting their timing changes."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from timing import BASE, OBJECT, derived, fit, obtain, signed_phase

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',default='timing-inputs',type=Path)
    p.add_argument('--output',default='chip.json',type=Path)
    a=p.parse_args();obtain(a.cache)
    urls={}
    def read_remote(name):
        path=a.cache/name
        url=BASE+'ogle4/OCVS/smc/ecl/'+name
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            with urlopen(url,timeout=45) as f: path.write_bytes(f.read())
        urls[name]={'url':url,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
        return path
    cat={l.split()[0]:l.split() for l in (a.cache/'ecl.dat').read_text().splitlines()}
    ident={l.split()[0]:l.split() for l in read_remote('ident.dat').read_text().splitlines()}
    chip='.'.join(ident[OBJECT][4].split('.')[:2])
    target_mag,target_p=float(cat[OBJECT][1]),float(cat[OBJECT][3])
    controls=[]
    for name,par in cat.items():
        ent=ident[name]
        if name==OBJECT or ent[1]!='NC' or '.'.join(ent[4].split('.')[:2])!=chip:continue
        if '-' in (par[1],par[3],par[4],par[5],par[6]):continue
        mag,per,zero,d1,d2=map(float,(par[1],par[3],par[4],par[5],par[6]))
        if (abs(mag-target_mag)>.75 or abs(np.log2(per/target_p))>1 or d1<.12 or d2<.1):continue
        path=read_remote('phot/I/'+name+'.dat')
        t,m,e=np.loadtxt(path).T
        rec={'id':name,'comparisons':[]}
        for use in [t<2457000,t>=2459700]:
            g=use&np.isfinite(t+m+e)&(e>0)&(e<=.1)
            if g.sum()<20:rec['comparisons'].append({'n':int(g.sum())});continue
            phase=signed_phase((t[g]-zero)/per)
            par=fit(phase,m[g],e[g],[np.quantile(m[g],.25),d1,0,.03,d2,.5,.03],'gaussian')
            counts=[int(np.sum(abs(signed_phase(phase-c))<.05)) for c in [par[2],par[5]]]
            rec['comparisons'].append({'n':int(g.sum()),'core_obs':counts,'separation':derived(par)[2],
                'center_phases':[par[2],par[5]],'fitted_depths':[par[1],par[4]]})
        controls.append(rec)
    a.output.write_text(json.dumps({'conditions':'Same OGLE field and chip as target, NC classification, |delta I| <= 0.75, 0.5 <= period ratio <= 2, catalog depths >= 0.12 / 0.10 mag. Count all matches. Fits use same OGLE-IV early and late windows, loss, error cut and epoch-specific freely adjusted parameters as target.',
        'inputs':urls,'controls':controls},indent=2)+'\n')
if __name__=='__main__':main()
