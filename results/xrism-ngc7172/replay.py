"""Offline replay of predeclared primary count test from compact event/GTI subset."""
import csv,gzip,json,math
from pathlib import Path
def fisher(a,b,c,d):
 n=a+b+c+d; k=a+c; row=a+b
 def lc(n,k): return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
 def prob(x): return math.exp(lc(k,x)+lc(n-k,row-x)-lc(n,row))
 p=prob(a);return min(1,sum(prob(x) for x in range(max(0,row-(n-k)),min(row,k)+1) if prob(x)<=p*(1+1e-10)))
p=Path(__file__).resolve().parent; g=json.loads((p/'gti.json').read_text());glob=g['global_gti']; pg=g['pixel_gti']; total=sum(b-a for a,b in glob);remaining=total/2
for a,b in glob:
 if remaining<=b-a: cut=a+remaining;break
 remaining-=b-a
by={i:[] for i in range(36)}
for pix,a,b in pg:by[pix].append((a,b))
counts=[[0,0],[0,0]]; n=0
with gzip.open(p/'events.csv.gz','rt',newline='') as f:
 for r in csv.DictReader(f):
  t=float(r['TIME']);pix=int(r['PIXEL']);itype=int(r['ITYPE']);pi=int(r['PI'])
  if pix in (12,27) or itype!=0 or not any(a<=t<b for a,b in by[pix]):continue
  half=int(t>=cut)
  if 12620<=pi<12760:counts[half][0]+=1
  if 12400<=pi<12540 or 12840<=pi<12980:counts[half][1]+=1
  n+=1
ref=json.loads((p/'result.json').read_text());assert counts==ref['counts'] and n==ref['eligible']-ref['pixel_gti_rejected'];assert math.isclose(fisher(*counts[0],*counts[1]),ref['p_value'],abs_tol=1e-12)
print('Verified offline:',counts,'eligible',n,'p',ref['p_value'])
