"""Refit the eclipse geometry with the leading equilibrium ellipsoidal amplitude imposed."""
import json
import numpy as np
from geometry import HERE,LightCurve,load,fit
from ellipsoidal import expectation

class TideCurve(LightCurve):
    def __init__(self,data,u,third,sense,tau):
        super().__init__(data,u,third,sense,order=32)
        self.tau=tau
    def predict(self,pars):
        q=pars.copy()
        a=expectation(q,self.third,self.u,self.tau)
        q[11]=a['minimum_cos2_mmag']/1000
        return super().predict(q)

def main():
    data=load(HERE/'inputs')
    full=json.loads((HERE/'dilution.json').read_text())['1']
    out=[]
    for row in full:
        if row['third'] not in (0.,.30):continue
        for u,tau in ((.2,0.),(.4,.32),(.4,1.)):
            c=TideCurve(data,u,row['third'],1,tau)
            res=fit(c,start=row['pars'],fixed={11:0.})
            pars=np.array(res['pars'])
            res['coefficient_mmag']=expectation(pars,row['third'],u,tau)['minimum_cos2_mmag']
            ctrl=fit(LightCurve(data,u,row['third'],1.,order=32),start=row['pars'])
            res.update(u=u,tau=tau,assumed_dilution=row['third'],free_harmonic_loss=ctrl['loss'],
                       increase_in_loss=res['loss']-ctrl['loss'])
            print(row['third'],u,tau,res['coefficient_mmag'],res['increase_in_loss'],flush=True)
            out.append(res)
    (HERE/'joint_tide.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
