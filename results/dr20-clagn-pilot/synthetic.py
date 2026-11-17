"""Synthetic-only response of the fixed-wing proxy; not an observational test."""
import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(20261116);x=np.arange(4750.,5100.,1.);left=(x>=4750)&(x<4770);right=(x>=5080)&(x<5100);wing=((x>=4805)&(x<4840))|((x>=4885)&(x<4920));shape=np.exp(-.5*((x-4862.68)/20.)**2)
def proxy(y):
 a=np.mean(y[left]);b=np.mean(y[right]);line=a+(b-a)*(x[wing]-np.median(x[left]))/(np.median(x[right])-np.median(x[left]));return float(np.sum(y[wing]-line))
rows=[]
for amp in (0,1,3):
 d=[]
 for i in range(200):
  a=10+rng.normal(size=len(x));b=10+amp*shape+rng.normal(size=len(x));d.append(proxy(b)-proxy(a))
 q=np.percentile(d,[5,50,95]);rows.append(dict(injected_peak=amp,n=200,p05=float(q[0]),median=float(q[1]),p95=float(q[2])))
Path(__file__).with_name('synthetic.json').write_text(json.dumps(rows,indent=2)+'\n');print(rows)
