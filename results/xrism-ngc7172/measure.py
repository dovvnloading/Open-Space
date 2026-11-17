"""Frozen split-count test; requires numpy and astropy; see PROTOCOL.md."""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
from astropy.io import fits
EVENT_SHA='74230a4c915f9fa7d7ff54bb1d09e242766e1dd1f5af67db3b6d0bf0ad26224c'
GTI_SHA='6cf9aedbb903fb8ff2289f835f70a1dfcdd45d2ed0a825f19f82547d44b23d81'
def verify(p,sha):
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha

def fisher(a,b,c,d):
 n=a+b+c+d; k=a+c; row=a+b
 def logcomb(n,k): return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
 def prob(x): return math.exp(logcomb(k,x)+logcomb(n-k,row-x)-logcomb(n,row))
 p=prob(a);return min(1,sum(prob(x) for x in range(max(0,row-(n-k)),min(row,k)+1) if prob(x)<=p*(1+1e-10)))

def main(root,out):
 ep=root/'xa202074010rsl_p0px1000_cl.evt.gz';gp=root/'xa202074010rsl_px1000_exp.gti.gz';verify(ep,EVENT_SHA);verify(gp,GTI_SHA)
 with fits.open(ep) as h, fits.open(gp) as gh:
  e=h['EVENTS'].data; g=h['GTI'].data; pg=gh['GTIPIXEL'].data
  starts,stops=g['START'],g['STOP']; durations=stops-starts; target=durations.sum()/2; idx=np.searchsorted(np.cumsum(durations),target);cut=float(starts[idx]+target-durations[:idx].sum())
  eligible=(e['ITYPE']==0)&(e['PIXEL']!=12)&(e['PIXEL']!=27)
  pixelok=np.zeros(len(e),dtype=bool); lives=[]
  for pix in range(36):
   if pix in (12,27): continue
   rows=pg[pg['PIXEL']==pix]; mask=e['PIXEL']==pix; tt=e['TIME'][mask]; ok=np.zeros(len(tt),dtype=bool)
   for x in rows:ok|=(tt>=x['START'])&(tt<x['STOP'])
   pixelok[mask]=ok
   for side,(lo,hi) in enumerate([(float(starts[0]),cut),(cut,float(stops[-1]))]):
    lives.append(dict(pixel=pix,half=side+1,exposure_s=float(np.maximum(0,np.minimum(rows['STOP'],hi)-np.maximum(rows['START'],lo)).sum())))
  mask=eligible&pixelok; t=e['TIME']; pi=e['PI']; counts=[]
  for lo,hi in [(starts[0],cut),(cut,stops[-1])]:
   m=mask&(t>=lo)&(t<hi);line=int(np.sum(m&(pi>=12620)&(pi<12760)));side=int(np.sum(m&(((pi>=12400)&(pi<12540))|((pi>=12840)&(pi<12980)))));counts.append([line,side])
  result=dict(cut_met_tt=cut,global_gti_s=float(durations.sum()),eligible=int(eligible.sum()),pixel_gti_rejected=int((eligible&~pixelok).sum()),counts=counts,p_value=fisher(*counts[0],*counts[1]),pixel_exposures=lives)
  out.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps({k:v for k,v in result.items() if k!='pixel_exposures'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('inputs',type=Path);p.add_argument('output',type=Path);a=p.parse_args();main(a.inputs,a.output)
