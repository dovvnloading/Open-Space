"""Extrapolate each telescope's astrometry of WISE J0615+1526 independently."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forecast as shared

HERE = Path(__file__).resolve().parent


def main():
    obs = pd.read_csv(HERE / "astrometry.csv")
    origin = pd.Series(dict(ra=93.988341, dec=15.439648, pmra=-29.1, pmdec=-532.1, parallax=52.8))
    bg = pd.read_csv(HERE/"gaia.csv", dtype={"source_id":str}).iloc[0]
    state = shared.params(bg, origin, 2016.0)
    state[2] = max(0.0, state[2])
    source_cov = shared.covariance(bg)
    times = shared.grid(2037, 2044, .5)
    factors = -shared.unit_vectors(origin.ra, origin.dec)[:2] @ shared.earth_position(times).T
    results = {}
    for label, train in [("Spitzer and WISE", np.ones(len(obs), bool)), ("Spitzer only", obs.band=="ch2"),
                         ("WISE only", obs.band=="W2")]:
        fit, cov, resid = shared.fit_spitzer(obs, origin, train)
        rel = state-fit
        d = np.hypot(*shared.tracks(rel,times,factors))
        i = np.argmin(d)
        draw = np.random.default_rng(12345)
        lens = draw.multivariate_normal(fit,cov,size=5000)
        # Gaia published solution has negative parallax. Retain covariances, then clip to the distant-source boundary.
        raw_source = shared.params(bg,origin,2016.0)
        source = draw.multivariate_normal(raw_source,source_cov,size=len(lens))
        source[:,2] = np.maximum(0,source[:,2])
        minima=[]
        for chunk in np.array_split(source-lens,50):
            sep=np.hypot(*np.moveaxis(shared.tracks(chunk,times,factors),1,0))
            minima.extend(sep.min(axis=1).tolist())
        results[label]=dict(epochs=int(sum(train)), lens_parameters=fit.tolist(), parallax_pm_errors=np.sqrt(np.diag(cov))[2:].tolist(),
                            closest_date=shared.iso(times[i]), closest_mas=float(d[i]),
                            peak_shift_at_0p05_solar_masses_mas=float(shared.centroid_shift(d,-rel[2],.05).max()),
                            minimum_mas_5_50_95_percentile=np.percentile(minima,[5,50,95]).tolist(),
                            rms_fit_mas=np.sqrt(np.mean(resid[train]**2,axis=0)).tolist())
        if fit[2] <= state[2]:
            results[label]["peak_shift_at_0p05_solar_masses_mas"] = None
            results[label]["status"] = "Negative relative parallax; geometry only, no physical lensing prediction."
        if not np.all(train):
            results[label]["other_telescope_prediction_rms_mas"]=np.sqrt(np.mean(resid[~train]**2,axis=0)).tolist()
    (HERE/"fits.json").write_text(json.dumps(results,indent=2)+"\n")
    for label,r in results.items():
        print(label,r['closest_date'],round(r['closest_mas']),r['peak_shift_at_0p05_solar_masses_mas'], np.round(r['parallax_pm_errors'],2),np.round(r['minimum_mas_5_50_95_percentile']))
if __name__=="__main__":
    main()
