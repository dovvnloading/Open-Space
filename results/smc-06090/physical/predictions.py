"""Predict primary/secondary eclipse widths at TESS epochs and July 2027 using only OGLE."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from geometry import HERE, LightCurve,P,TREF,load

def evaluate(curve,q,jd,flux):
    t=(jd-TREF)/365.25
    q=q.copy(); q[11:13]=0
    def depth(phase):
        m=curve.predict_at(q,np.array([t]),np.array([phase]))[0]
        return 1-10**(-.4*m) if flux else m
    pair=[]
    for seed in [0,.5]:
        opt=minimize_scalar(lambda x:-depth(x),bounds=(seed-.114,seed+.13),method='bounded')
        c=opt.x;half=-opt.fun/2
        l=brentq(lambda x:depth(x)-half,c-.19,c)
        r=brentq(lambda x:depth(x)-half,c,c+.19)
        pair.append([(r-l)*P*1440,c])
    return dict(primary_min=pair[0][0],secondary_min=pair[1][0],difference_min=pair[1][0]-pair[0][0],
                separation_hr=(pair[1][1]-pair[0][1])*P*24)
def main():
    inp=json.loads((HERE/'inputs'/'tess_timing.json').read_text())['groups']
    times={label:v['hjd_mid_approx'] for label,v in inp.items()}
    times['2027_July_01']=2461587.5
    data=load(HERE/'inputs')
    out=[]
    for sense,items in json.loads((HERE/'dilution.json').read_text()).items():
        for row in items:
            if row['third']>.301:continue
            for u in (.2,.4,.6):
                c=LightCurve(data,u,row['third'],float(sense),order=48)
                for k,jd in times.items():
                    q=np.array(row['pars'])
                    out.append(dict(sense=int(sense),third=row['third'],u=u,epoch=k,jd=jd,
                                flux=evaluate(c,q,jd,True),mag=evaluate(c,q,jd,False)))
    (HERE/'width_predictions.json').write_text(json.dumps(out,indent=2)+'\n')
    for label in times:
        for sense in (1,-1):
            vals=[v['flux']['difference_min'] for v in out if v['sense']==sense and v['epoch']==label]
            print(label,sense,min(vals),max(vals))
if __name__=='__main__':main()
