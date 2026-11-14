"""Search MAXI timing bins for the pulse frequency and its time derivative."""
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
START = 61132
FREQUENCIES = np.linspace(1890,1910,801)  # microhertz
SLOPES = np.linspace(-6,3,73)            # microhertz per Julian year


def search(table,lo,hi,slopes):
    epoch = (lo+hi)/2
    times = START+(16*table[:,0]+8)/86400
    mask = (times>=lo)&(times<hi)
    seconds = (times[mask]-epoch)*86400
    n,m = table[mask,1:].sum(axis=0)
    weights = table[mask,1]/n-table[mask,2]/m
    scale = 2/(1/n+1/m)
    powers = np.empty((2,len(slopes),len(FREQUENCIES)))
    for h in (1,2):
        drift = np.exp(1j*np.pi*h*seconds[None,:]**2*slopes[:,None]*1e-6/(365.25*86400))
        drift = (drift*weights[None,:]).astype(np.complex64)
        for a in range(0,len(FREQUENCIES),64):
            band = FREQUENCIES[a:a+64]
            phase = np.exp(2j*np.pi*h*seconds[:,None]*band[None,:]*1e-6).astype(np.complex64)
            powers[h-1,:,a:a+len(band)] = scale*abs(drift@phase)**2

    def pick(score):
        i,j = np.unravel_index(np.argmax(score),score.shape)
        return {"frequency_uHz":round(float(FREQUENCIES[j]),4),"slope_uHz_per_year":float(slopes[i]),
                "fundamental_power":round(float(powers[0,i,j]),2),
                "second_harmonic_power":round(float(powers[1,i,j]),2)}

    best = pick(powers[0])
    line = np.flatnonzero(slopes==0)[0]
    j = int(powers[0,line].argmax())
    return {"epoch_mjd":epoch,"source_photons":int(n),"annulus_photons":int(m),
            "best_fundamental":best,"best_harmonic_geometric_mean":pick(np.sqrt(powers[0]*powers[1])),
            "best_sum_two_harmonics":pick(powers.sum(axis=0)),
            "best_constant_frequency_uHz":round(float(FREQUENCIES[j]),4),
            "best_constant_fundamental_power":round(float(powers[0,line,j]),2)}


def daily_coherence(table):
    # This test specifies the phase evolution using GBM alone. Its frequency
    # and slope are not selected by maximizing a MAXI periodogram.
    import csv
    with (ROOT/"gbm.csv").open() as stream:
        history = [(float(r["barytime"]),1e6*float(r["frequency"])) for r in csv.DictReader(stream)
                   if float(r["barytime"])>=61131]
    t,f = np.array(history).T
    slope,nu = np.polyfit((t-61228.5)/365.25,f,1)
    seconds = (table[:,0]+.5)*16
    day = (seconds/86400).astype(int)
    dt = seconds-96.5*86400
    phase = 2*np.pi*(dt*nu*1e-6 + .5*dt**2*slope*1e-6/(365.25*86400))
    n,m = table[:,1:].sum(axis=0)
    weights = table[:,1]/n-table[:,2]/m
    vectors = weights*np.exp(1j*phase)
    by_day = np.bincount(day,vectors.real,minlength=193)+1j*np.bincount(day,vectors.imag,minlength=193)
    scale = 2/(1/n+1/m)
    observed = abs(by_day.sum())**2*scale
    rng = np.random.default_rng(74382)
    exceed = 0
    for _ in range(100):
        randomized = np.exp(2j*np.pi*rng.random((1000,193)))@by_day
        exceed += np.count_nonzero(abs(randomized)**2*scale>=observed)
    return {"gbm_frequency_uHz":round(float(nu),6),"gbm_slope_uHz_per_year":round(float(slope),6),
            "observed_fundamental_power":round(float(observed),2),"randomizations":100000,
            "as_large":int(exceed),"seed":74382}


if __name__ == "__main__":
    with gzip.open(ROOT/"maxi_bins.csv.gz","rt") as stream:
        data = np.loadtxt(stream,delimiter=",",skiprows=1,dtype=int)
    results = {}
    for lo,hi,slopes in [(61132,61325,SLOPES),(61132,61228,np.array([0.])),
                         (61228,61325,np.array([0.]))]:
        results[f"{lo}-{hi}"] = search(data,lo,hi,slopes)
    results["fixed_gbm_track"]=daily_coherence(data)
    print(json.dumps(results,indent=2))
