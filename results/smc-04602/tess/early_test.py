"""Historical sector test, initialized at interpolated OGLE-era predictions.

Use the same wide free-center search window as the 2025 analysis. 144-bin
coverage was designed for 200-s data; 30-min data use a bandpass coverage cut.
"""
from pathlib import Path
import sys,json
import numpy as np
import test_tess as t
if __name__=='__main__':
    out={}
    for year,sectors,centers in [(2018,[1],[-.00750,.60793]),(2020,[27,28],[-.00588,.61035]),(2023,[67,68],[.00084,.61597])]:
        fs=[];metas=[]
        for s in sectors:
            meta,folds=t.reduce(Path(f's{s:03}.fits'))
            fs+=folds;metas.append(meta)
        if fs:
            ar=[np.stack(a) for a in zip(*fs)]
            fit=t.evalf(ar,centers=np.array(centers));rng=np.random.default_rng(4602)
            bs=np.array([t.evalf([a[idx] for a in ar],centers=np.array(centers))[5:7] for idx in rng.integers(len(fs),size=(399,len(fs)))])
            out[year]={'n':len(fs),'fixed_baseline':centers,'centers':fit[5:7].tolist(),'centers95':np.quantile(bs,[.025,.975],axis=0).tolist(), 'separation':fit[6]-fit[5],'sep95':np.quantile(bs[:,1]-bs[:,0],[.025,.975]).tolist(),'meta':metas}
        else: out[year]={'n':0,'meta':metas}
    print(json.dumps(out,indent=2))
