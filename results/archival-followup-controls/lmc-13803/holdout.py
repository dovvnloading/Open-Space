"""Measure coverage and, where eligible, depths in three disjoint epochs."""
import json
import numpy as np
import screen

inputs = screen.ROOT / 'inputs'
data = np.loadtxt(inputs / 'ogle4-I.dat')
row = (inputs / 'catalog-row.txt').read_text().split()
period, epoch = map(float, row[3:5])
records={}
for label,limits in [('early',screen.EARLY),('training',screen.TRAIN),('holdout',screen.HOLD)]:
    t,m,e=data.T
    keep=(t>=limits[0])&(t<limits[1])&np.isfinite(t+m+e)&(e>0)&(e<=.10)
    ph=np.abs(((t-epoch)/period+.5)%1-.5)
    record={'observations':int(keep.sum())}
    for feature,w in [('quadrature',(ph>.16)&(ph<.34)),('primary',ph<.03),('secondary',ph>.47)]:
        record[feature+'_nights']=len(np.unique(np.floor(t[keep&w])))
    depths=screen.per_era(data,limits,period,epoch)
    record['core_minus_quadrature_mag']=None if depths is None else [round(s[0]-depths[0][0],6) for s in depths[1:]]
    records[label]=record
(screen.ROOT/'holdout.json').write_text(json.dumps(records,indent=2)+'\n')
