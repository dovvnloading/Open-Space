"""Measure contaminating flux from KIC 6367628 in the two KIC 6367663 stamps."""
from pathlib import Path
import hashlib
from urllib.request import urlopen
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
PERIOD = 3.779903762546505

def product(name):
    cache = HERE/'.cache'
    cache.mkdir(exist_ok=True)
    path=cache/name
    row=pd.read_csv(HERE/'input_manifest.csv').set_index('filename').loc[name]
    if not path.exists():
        path.write_bytes(urlopen(row.url,timeout=120).read())
    if hashlib.sha256(path.read_bytes()).hexdigest()!=row.sha256:
        raise ValueError(f"Input checksum mismatch for {name}")
    return path

def lightcurve(name):
    with fits.open(product(name)) as h:
        d=h[1].data
        good=(d['SAP_QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d['PDCSAP_FLUX'])
        return d['CADENCENO'][good].astype(int),d['TIME'][good].astype(float),d['PDCSAP_FLUX'][good].astype(float)

def orbit_fit(t, predictor, measurements, degree, rng):
    orbit=np.floor((t-t.min())/PERIOD).astype(int)
    numerator,denominator,folds=[],[],[]
    for number in np.unique(orbit):
        use=orbit==number
        if use.sum()<10:continue
        local=t[use]-t[use].mean()
        basis=np.vander(local,degree+1,increasing=True)
        template=predictor[use]-basis@np.linalg.lstsq(basis,predictor[use],rcond=None)[0]
        detrended=measurements[use]-basis@np.linalg.lstsq(basis,measurements[use],rcond=None)[0]
        numerator.append(template@detrended)
        denominator.append(template@template)
        folds.append(number%2)
    num,den=np.asarray(numerator),np.asarray(denominator)
    sample=rng.integers(len(den),size=(1000,len(den)))
    boot=num[sample].sum(axis=1)/den[sample].sum(axis=1)[:,None]
    slope=num.sum(axis=0)/den.sum()
    halves=[num[np.array(folds)==f].sum(axis=0)/den[np.array(folds)==f].sum() for f in range(2)]
    return slope,boot,halves

def main():
    records, pixels=[],[]
    rng=np.random.default_rng(5275)
    for quarter,stamp_id in [(14,'2012277125453'),(15,'2013011073258')]:
        tpf=f'kplr006367663-{stamp_id}_lpd-targ.fits.gz'
        sn,st,s=lightcurve(f'kplr006367628-{stamp_id}_llc.fits')
        tn,_,tf=lightcurve(f'kplr006367663-{stamp_id}_llc.fits')
        _,i,j=np.intersect1d(tn,sn,return_indices=True)
        source=s[j]/np.median(s[j])-1
        target=tf[i]/np.median(tf[i])-1
        b,draws,halves=orbit_fit(st[j],source,target[:,None],0,rng)
        low,high=np.quantile(draws[:,0],[.025,.975])
        rec=dict(quarter=quarter,n_common_lc=len(i),transfer_ppm=b[0]*1e6,lower=low*1e6,upper=high*1e6,
                 odd_even_A=halves[0][0]*1e6,odd_even_B=halves[1][0]*1e6)
        with fits.open(product(tpf)) as h:
            data=h[1].data
            good=(data['QUALITY']==0)&np.isfinite(data['TIME'])
            n=data['CADENCENO'][good].astype(int)
            stamp=data['FLUX'][good].astype(float)
            finite=np.all(np.isfinite(stamp),axis=0)
            y,x=np.where(finite)
            optimal=(h[2].data[finite]&2)>0
            w=WCS(h[2].header)
            target_xy=w.world_to_pixel_values(h[0].header['RA_OBJ'],h[0].header['DEC_OBJ'])
            source_xy=w.world_to_pixel_values(292.98524259893,41.71687270018)
        _,i,j=np.intersect1d(n,sn,return_indices=True)
        source=s[j]/np.median(s[j])-1
        b,draws,_=orbit_fit(st[j],source,stamp[i][:,finite],2,rng)
        e=draws.std(axis=0)
        best=np.argmax(b/e)
        rec.update(common_pixel_cadences=len(i),highest_snr_x=x[best],highest_snr_y=y[best],
                   corner_slope=b[best],corner_bootstrap_se=e[best],target_x=target_xy[0],target_y=target_xy[1],
                   source_x=source_xy[0],source_y=source_xy[1])
        for col,row,sl,std,ap in zip(x,y,b,e,optimal):
            pixels.append(dict(quarter=quarter,x=col,y=row,slope_e_s=sl,bootstrap_se=std,optimal_aperture=ap))
        records.append(rec)
    pd.DataFrame(records).to_csv(HERE/'quarter_results.csv',index=False)
    pd.DataFrame(pixels).to_csv(HERE/'pixel_transfer.csv',index=False)
    print(pd.DataFrame(records)[['quarter','n_common_lc','transfer_ppm','lower','upper','corner_slope','corner_bootstrap_se']].round(2).to_string(index=False))

if __name__=='__main__':main()
