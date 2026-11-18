"""Estimate the total baseline in pySIS differential-flux units from the supplied magnitudes."""
import tarfile
import numpy as np


def stream_flux_scale(data):
    """Return (baseline, robust fractional residual) from dF, mag, mag_err columns.

    pySIS uses fdiff = -Delta_flux = alpha*10**(-0.4*(I-Iref)) + offset.
    Infer alpha from robustly separated quantiles to avoid rounding at one epoch.
    This estimates total light on the survey's calibration, including blends.
    """
    delta, mag, mag_err = data.T
    good = np.isfinite(data).all(axis=1) & (mag>8) & (mag<27) & (mag_err>0) & (mag_err<0.5)
    if np.count_nonzero(good)<15:
        return np.nan, np.nan
    mag=mag[good]
    f=-delta[good]
    scale=10**(-0.4*(mag-np.median(mag)))
    lower,upper=np.quantile(scale,[0.15,0.85])
    lo,hi=scale<=lower,scale>=upper
    span=np.median(scale[hi])-np.median(scale[lo])
    if span<=0:
        return np.nan,np.nan
    alpha=(np.median(f[hi])-np.median(f[lo]))/span
    offset=np.median(f-alpha*scale)
    baseline=np.median(f)-offset
    scatter=np.median(abs(f-offset-alpha*scale))/max(baseline,1e-12)
    return (float(baseline),float(scatter)) if baseline>0 and scatter<0.01 else (np.nan,float(scatter))


def baseline_caps(path):
    """Map I-band tar member name to total baseline flux; failed calibrations have NaN."""
    caps={}
    with tarfile.open(path) as tar:
        for m in tar:
            if m.isfile() and m.name.endswith('_I.pysis'):
                a=np.loadtxt(tar.extractfile(m),usecols=(1,3,4),ndmin=2)
                if len(a):
                    caps[m.name]=stream_flux_scale(a)[0]
    return caps
