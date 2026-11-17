"""Check every indexed QLP visit to TIC 316937197 for repeat dimmings."""
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
import numpy as np
from acquire import fetch, product_url, sha256
from confirm import unpack
from scan import Binning, Curve, detect


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    tic=316937197
    params=dict(columns='obs_id,obs_collection,provenance_name,sequence_number,dataURL',
                filters=[dict(paramName='target_name',values=[str(tic)]),
                         dict(paramName='provenance_name',values=['QLP'])])
    request=dict(service='Mast.Caom.Filtered',params=params,format='json',pagesize=200,page=1)
    url='https://mast.stsci.edu/api/v0/invoke?request='+quote(json.dumps(request))
    response=fetch(url)
    data=json.loads(response)
    if data.get('status')!='COMPLETE':
        raise ValueError(data)
    visits=list(data['data'])
    if not any(v['sequence_number']==105 for v in visits):
        visits.append(dict(sequence_number=105,dataURL='mast:HLSP/'+
                           product_url(105,tic).split('/hlsps/')[1]))
    def examine(r):
        sector=r['sequence_number']
        url='https://archive.stsci.edu/hlsps/'+r['dataURL'].split('mast:HLSP/')[1]
        got=dict(sector=sector,url=url)
        try:
            raw=fetch(url)
            hd,table=unpack(raw)
            if int(hd['TICID'])!=tic or int(hd['SECTOR'])!=sector:
                raise ValueError('identity mismatch')
            got.update(sha256=sha256(raw),bytes=len(raw))
            names=table.dtype.names
            col=next(k for k in ['DET_FLUX','KSPSAP_FLUX'] if k in names)
            c=Curve(table['TIME'].astype(float),table[col].astype(float),table['QUALITY']==0,{}, {})
            b=Binning(c); y=b.bin_flux(b.flux)
            e=detect(y,b.centers,require_shoulders=True)
            # Independent event-shape check for the deep, nearly flat-bottomed
            # S105 eclipse. Require an excursion >12% in >=2 consecutive
            # half-hour bins against each five-hour shoulder in 0.35–0.65 d.
            lows=[]
            for i in range(32,len(y)-32):
                if not np.isfinite(y[i:i+2]).all():
                    continue
                left=y[i-31:i-17];right=y[i+17:i+31]
                if np.isfinite(left).sum()<5 or np.isfinite(right).sum()<5:
                    continue
                depth=min(np.nanmedian(left),np.nanmedian(right))-np.max(y[i:i+2])
                if depth>.12:
                    if lows and b.centers[i]-lows[-1]['last']<.6:
                        lows[-1]['last']=b.centers[i]
                        if depth>lows[-1]['depth']:
                            lows[-1].update(peak=b.centers[i],depth=depth)
                    else:
                        lows.append(dict(peak=b.centers[i],depth=depth,last=b.centers[i]))
            return dict(**got,flux=col,first=float(b.lo),last=float(b.hi),good_bins=int(b.valid.sum()),
                        dips=e,deep_epochs=lows, status='ok'),(b.centers,y)
        except Exception as exc:
            return dict(**got,status='failed',detail=str(exc)),None
    results=[];series=[]
    with ThreadPoolExecutor(5) as pool:
        for r,b in pool.map(examine,visits):
            results.append(r)
            if b is not None:
                series.append((r['sector'],*b))
    out=dict(tic=tic,query_url=url,response_sha256=sha256(response),
             query_response=data,products=sorted(results,key=lambda r:r['sector']))
    # Test the apparent separation between the two S105 detector triggers.
    epoch=4217.63199
    p=10.54167
    phase=[]
    for s,t,y in series:
        for k in range(int(np.ceil((t[0]-epoch)/p)), int(np.floor((t[-1]-epoch)/p))+1):
            center=epoch+k*p
            w=np.abs(t-center)<.12
            sh=(np.abs(t-center)>.35)&(np.abs(t-center)<.65)
            if np.isfinite(y[w]).sum()>=5 and np.isfinite(y[sh]).sum()>=12:
                phase.append(dict(sector=s,k=k,center=center,core_median=float(np.nanmedian(y[w])),
                             shoulder_median=float(np.nanmedian(y[sh]))))
    out.update(same_sector_pair_period_trial=p,predicted_epochs=phase)
    # Four consecutive minima near 4161/4180/4199/4218 and a 2018 minimum
    # near 1327.48 imply ~154 cycles in the intervening baseline. Estimate an
    # ephemeris using the timing of each covered return to that phase.
    ref, trial = 4217.572, 18.76683
    timings = []
    for sector,t,y in series:
        for k in range(int(np.ceil((t[0]-ref)/trial)),int(np.floor((t[-1]-ref)/trial))+1):
            pred=ref+k*trial
            a=np.isfinite(y)&(t>=pred-.55)&(t<pred-.30)
            b=np.isfinite(y)&(t>pred+.30)&(t<=pred+.55)
            c=np.isfinite(y)&(np.abs(t-pred)<=.22)
            if min(a.sum(),b.sum())<6 or c.sum()<16:
                continue
            baseline=np.interp(t[c],[np.median(t[a]),np.median(t[b])],[np.median(y[a]),np.median(y[b])])
            deficit=baseline-y[c]
            depth=np.median(np.sort(deficit)[-3:])
            noise=1.4826*np.median(abs(y[a|b]-np.median(y[a|b])))
            w=np.maximum(deficit-.2*depth,0.)
            center=np.average(t[c],weights=w) if np.sum(w)>0 else pred
            timings.append(dict(sector=sector,cycle=k,btjd=float(center),empirical_depth=float(depth),
                                shoulder_mad=float(noise),accepted=bool(depth>max(.002,4*noise))))
    clock=[r for r in timings if r['accepted']]
    if len(clock)>=3:
        p,e=np.polyfit([r['cycle'] for r in clock],[r['btjd'] for r in clock],1)
        residuals=np.array([r['btjd']-(e+p*r['cycle']) for r in clock])
        out['refined_clock']=dict(period_d=p,epoch_btjd=e,n_timings=len(clock),
                         baseline_d=max(r['btjd'] for r in clock)-min(r['btjd'] for r in clock),
                         rms_residual_d=float(np.sqrt(np.mean(residuals**2))))
    out['timings']=sorted(timings,key=lambda r:r['cycle'])
    args.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
