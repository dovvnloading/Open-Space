"""Replay the blinded spectral ratio and physically normalized flux injections."""
import json
import warnings
import numpy as np
from catalogs import HERE, rows
from screen import load,save_csv

BANDS={'H':(1.51,1.72),'M':(3.16,3.62),'L':(4.42,4.90)}
def spread(v):
    return 1.4826*np.nanmedian(abs(v-np.nanmedian(v)))


def medians(y,wl,target):
    bands={}
    for name,(lo,hi) in BANDS.items():
        sel=(wl>=lo)&(wl<=hi)
        channels=np.isfinite(target)&sel
        need=max(2,int(np.ceil(sel.sum()/2)))
        if channels.sum()<need:
            bands[name]=np.full(y.shape[0],np.nan)
            continue
        masked=y[:,channels]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            med=np.nanmedian(masked,axis=1)
        med[np.isfinite(masked).sum(axis=1)<need]=np.nan
        bands[name]=med
    return bands


def brightness(wl,t):
    return wl**-3/np.expm1(14387.76877/(wl*t))


def measure(y,wl):
    d=medians(y,wl,y[0])
    h,m,l=(d[k][0] for k in ('H','M','L'))
    bh=d['H'][1:];bl=d['L'][1:]
    good=np.isfinite(bh)&np.isfinite(bl)
    if good.sum()<8 or not np.isfinite(h) or h<=0 or not np.isfinite(l):
        return None
    eh,el=spread(bh[good]),spread(bl[good])
    if not (eh>0 and el>0 and h/eh>=5):
        return None
    return dict(H=h,M=m,L=l,H_snr=h/eh,L_snr=l/el,q=l/h,bh=bh[good],bl=bl[good],blank_n=int(good.sum()))


def peers(ids,i,measurements,color,g):
    ix=[k for k,m in enumerate(measurements) if k!=i and m is not None
        and abs(color[k]-color[i])<=.5 and abs(g[k]-g[i])<=1.5]
    ix.sort(key=lambda k:(((color[k]-color[i])/.5)**2+((g[k]-g[i])/1.5)**2,ids[k]))
    ix=ix[:15]
    if len(ix)<4: return None
    ratios=np.array([measurements[k]['q'] for k in ix])
    center=np.median(ratios)
    if center<=0: return None
    return dict(ix=ix,center=center,scatter=max(.05*center,spread(ratios)))


def score(m,peer):
    center,scatter,n=peer['center'],peer['scatter'],len(peer['ix'])
    blank=m['bl']-center*m['bh']
    sigma=np.sqrt(spread(blank)**2+(m['H']*scatter)**2*(1+1.253**2/n))
    z=(m['L']-center*m['H']-np.median(blank))/sigma
    ratio=m['L']/m['H']/center
    return z,ratio,bool(ratio>=1.25 and z>=5 and m['L_snr']>=5)


def main():
    cohort=rows('dr1_positions.csv')
    ids=[r['source_id'] for r in cohort]
    color=np.array([float(r['bp_rp']) for r in cohort])
    g=np.array([float(r['phot_g_mean_mag']) for r in cohort])
    actual={s['source_id']:s for s in load()}
    dataset=np.load(HERE/'dr1_spectra_and_blanks.npz')
    wl=dataset['wavelength_um']
    full={a:dataset[a].astype(float) for a in ids}
    results=[];injection=[];summary={}
    flags={name:set() for name in ('strict','outlier_included')}
    for mask_index,name in enumerate(flags):
        spectra=[np.vstack([full[a][mask_index],full[a][2:]]) for a in ids]
        measurements=[measure(y,wl) for y in spectra]
        matches=[peers(ids,i,measurements,color,g) for i in range(len(ids))]
        controls=[];tested=0
        for i,(m,peer) in enumerate(zip(measurements,matches)):
            rec=dict(source_id=ids[i],extraction=name,status='',h_mJy='',m_mJy='',l_mJy='',
                     blank_n='',h_snr='',l_snr='',ratio_over_peers='',z='',comparison_ids='')
            if m is None: rec['status']='insufficient channels, sky, or H signal'
            elif peer is None: rec['status']='insufficient peer model'
            else:
                z,ratio,flag=score(m,peer)
                tested+=1
                rec.update(status='trigger' if flag else 'screened',h_mJy=m['H'],m_mJy=m['M'],l_mJy=m['L'],
                           blank_n=m['blank_n'],h_snr=m['H_snr'],l_snr=m['L_snr'],ratio_over_peers=ratio,z=z,
                           comparison_ids=' '.join(ids[j] for j in peer['ix']))
                if flag: flags[name].add(ids[i])
                t=float(actual[ids[i]].get('teff_gspphot') or np.nan)
                if abs(z)<3 and 2500<=t<=8000: controls.append((i,t))
            results.append(rec)
        for td in [100,200,300,400,500,600]:
            for gamma in [0.,.001,.01,.03,.1]:
                hits=0
                for i,t in controls:
                    y=spectra[i];m=measurements[i]
                    hmask=(wl>=1.51)&(wl<=1.72)&np.isfinite(y[0])
                    norm=m['H']/np.median(brightness(wl[hmask],t))
                    added=gamma*(t/td)**4*norm*brightness(wl,td)
                    altered=y.copy();altered[0]+=added
                    test=measure(altered,wl)
                    if test is not None: hits+=int(score(test,matches[i])[2])
                injection.append(dict(extraction=name,dust_K=td,L_excess_over_L_star=gamma,
                                      injected=len(controls),recovered=hits))
        summary[name]=dict(attempted=len(ids),valid_target_photometry=sum(m is not None for m in measurements),
                           tested=tested,triggers=sorted(flags[name]),injection_controls=len(controls))
    summary['triggers_both_masks']=sorted(set.intersection(*flags.values()))
    save_csv('dr1_screen.csv',results);save_csv('dr1_injections.csv',injection)
    (HERE/'dr1_counts.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
