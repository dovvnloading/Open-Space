"""Aggregate color evidence by scan direction; keep dropout sensitivity separate from measured color."""
import argparse
import csv
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

def mid(x):
    return None if not len(x) else round(float(np.median(x)),3)
def main():
    p=argparse.ArgumentParser();p.add_argument('tracks',type=Path);p.add_argument('folder',type=Path);args=p.parse_args()
    trackinfo={r['source_name']:r for r in csv.DictReader(args.tracks.open())}
    bytrack=defaultdict(list)
    for r in csv.DictReader((args.folder/'paired_epochs.csv').open()):bytrack[r['track']].append(r)
    out=[]
    for name,rr in bytrack.items():
        track=trackinfo[name]
        def arr(k):return np.array([float(r[k]) for r in rr])
        f1,s1,f2,s2=map(arr,['w1_flux','w1_sigma','w2_flux','w2_sigma'])
        sep=arr('nearest_w1_arcsec');parity=arr('parity').astype(int);local=arr('local_sigma_p90')
        good=(sep<=2)&(arr('w1_flags')==0)&(arr('w1_qf')>=.85)&(s1>0)&(f1>=5*s1)
        color=2.5*np.log10(f2/f1)
        upper=2.5*np.log10((f2+2*s2)/(f1-2*s1),where=(f1>2*s1),out=np.full(len(f1),np.nan))
        lower=2.5*np.log10((f2-2*s2)/(f1+2*s1))
        proxy=2.5*np.log10((f2-2*s2)/(6*local))
        clean_empty=(sep>6)&(arr('n_local')>=10)&(arr('w2_fracflux')>=.70)&np.isfinite(proxy)
        result=dict(track=name,w2_epochs=int(track['n']),paired_epochs=len(rr),w2_catalog_mag=round(float(track['w2']),3),
                    w2_median_fracflux=mid(arr('w2_fracflux')),w2_p90_centroid_sigma_arcsec=round(float(np.quantile(
                    2.75*np.hypot(arr('w2_dx_pix'),arr('w2_dy_pix')),.9)),3),
                    w1_matches=int(good.sum()),median_measured_w1_w2=mid(color[good]))
        isblue=True;isred=True;issensitive=True
        for h in [0,1]:
            g=good&(parity==h); clean=clean_empty&(parity==h)
            result.update({f'n_w1_p{h}':int(g.sum()),f'median_color_p{h}':mid(color[g]),
                           f'n_clean_empty_p{h}':int(clean.sum()),f'median_local_depth_color_p{h}':mid(proxy[clean])})
            isblue &= g.sum()>=1 and np.median(upper[g])<2
            isred &= g.sum()>=1 and np.median(lower[g])>=2
            issensitive &= clean.sum()>=2 and np.median(proxy[clean])>=2
        if good.sum()>=3 and isblue:
            status='repeated_W1_below_red_cut'
        elif good.sum()>=3 and isred:
            status='repeat_red_measured'
        elif issensitive and np.sum(good & np.isfinite(upper) & (upper<2))==0:
            status='sensitive_catalog_dropout'
        else:
            status='inconclusive'
        result['color_status']=status
        out.append(result)
    out.sort(key=lambda x:(-x['w2_epochs'],x['track']))
    with (args.folder/'track_colors.csv').open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(out[0]));wr.writeheader();wr.writerows(out)
    print(Counter(x['color_status'] for x in out))
    print('non-blue with color support:',[(x['track'],x['w2_epochs'],x['color_status']) for x in out if x['color_status'] in ['repeat_red_measured','sensitive_catalog_dropout']])
if __name__=='__main__':main()
