"""Score the two auxiliary visits offline using the date-split rules frozen on GitHub."""
import csv
import gzip
from pathlib import Path
import numpy as np
from auxiliary import VISITS, clean, phase_score
from locked.score import depth

ROOT = Path(__file__).resolve().parent

def calculate(training, testing, parameters):
    obj,sector,split,period,duration,spin=parameters
    if min(len(training),len(testing))<40:
        return [obj,sector,len(training),len(testing),"insufficient"]
    t,r=clean(training,spin,split)
    tt,rr=clean(testing,spin,split)
    candidate=np.arange(100)/100
    phase=candidate[np.argmax([phase_score(t,r,period,duration,c) for c in candidate])]
    events=depth(tt,rr,period,duration,phase)
    if not len(events):
        return [obj,sector,len(t),len(tt),phase,0]
    d=events[:,1]
    low,high=np.quantile(np.median(np.random.default_rng(3839).choice(d,(10000,len(d))),axis=1),[.025,.975])
    return [obj,sector,len(t),len(tt),phase,len(d),int(sum(d>0)),np.median(d),low,high,1.4826*np.median(abs(rr-np.median(rr)))]

def main():
    with gzip.open(ROOT/"aux_difference_photometry.csv.gz","rt") as f:
        records=list(csv.DictReader(f))
    with (ROOT/"other_visits.csv").open() as f:
        expected=list(csv.reader(f))[1:]
    for key,want in zip(VISITS,expected,strict=True):
        parameters=VISITS[key]
        rows=[r for r in records if int(r["mpc"])==parameters[0] and int(r["sector"])==parameters[1]]
        pair=[]
        for subset in ("train","test"):
            run=[r for r in rows if r["subset"]==subset]
            array=np.zeros((len(run),9))
            for i,r in enumerate(run):
                array[i,[0,7,8]]=[float(r[k]) for k in ("jd_utc","reference_noise","prf_flux")]
            pair.append(array)
        result=calculate(*pair,parameters)
        if len(result)!=len(want) or not np.allclose(result,np.array(want,dtype=float),rtol=1e-8,atol=1e-8):
            raise ValueError(f"Mismatch in {key}: {result}")
        print(key,"verified",result[5:8])

if __name__=="__main__":
    main()
