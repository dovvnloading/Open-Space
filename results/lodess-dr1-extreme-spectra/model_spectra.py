"""Model the six positive LoLSS in-band slopes and retain the full high-band cohort."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import log_ndtr

# One datum per survey in each fit. LoLSS channels and the wideband point share photons.
BANDS = [("LoDeSS", 23, .39), ("LoLSS", 54, .10), ("VLSSr", 73.8, .15),
         ("LoTSS", 144, .07), ("TGSS", 150, .10), ("WENSS", 325, .08), ("NVSS", 1400, .05)]
PREFIX = {"VLSSr": "vlssr", "LoTSS": "lotss", "TGSS": "gmrtas150m", "WENSS": "wenss", "NVSS": "nvss"}


def photometry(row):
    rec = []
    for survey, frequency, scale in BANDS:
        if survey == "LoDeSS":
            flux, error = row.psf_jy, row.psf_noise_jy
        elif survey == "LoLSS":
            flux, error = row.Total_flux / 1000, row.E_Total_flux / 1000
        else:
            pre = PREFIX[survey]
            flux, error = row[pre + "_flux_jy"], row[pre + "_err_jy"]
        if np.isfinite(flux + error) and error > 0:
            limit = max(0, flux) + 1.645 * error if survey == "LoDeSS" and flux < 3 * error else np.nan
            rec.append((row.Source_name, survey, frequency, flux, error, scale, limit))
    return rec


def spectrum(parameters, x, model):
    intercept, alpha = parameters[:2]
    log_flux = intercept + alpha * x
    if model == "curved":
        log_flux += parameters[2] * x**2
    if model == "absorbed":
        log_flux -= parameters[2] * np.exp(-2.1 * x)
    return np.exp(np.clip(log_flux, -100, 30))


def fit(frame, model):
    x = np.log(frame.freq_mhz.to_numpy() / 144)
    flux, error, scale = (frame[col].to_numpy() for col in ["flux_jy", "stat_jy", "frac_scale"])
    # Independent Gaussian gain per survey, analytically marginalized. Each gain is shared
    # within that survey. One wideband point per survey avoids counting correlated subbands.
    def statistic(parameters):
        prediction = spectrum(parameters, x, model)
        variance = error**2 + (scale * prediction)**2
        return np.sum((flux - prediction)**2 / variance + np.log(variance))
    init, bounds = [np.log(np.median(flux[flux > 0])), -.7], [(-10, 8), (-5, 4)]
    if model == "curved":
        init.append(-.15)
        bounds.append((-3, 0))
    if model == "absorbed":
        init.append(.025)
        bounds.append((0, 5))
    result = minimize(statistic, init, method="Nelder-Mead", bounds=bounds,
                      options={"maxiter": 2000, "xatol": 1e-8, "fatol": 1e-8})
    if not result.success:
        raise RuntimeError(result.message)
    par = result.x
    peak = np.nan
    if model == "curved" and par[2] < -1e-6:
        xp = -par[1] / (2 * par[2])
        if -5 < xp < 6:
            peak = 144 * np.exp(xp)
    if model == "absorbed" and par[1] < 0 and par[2] > 0:
        peak = 144 * (-2.1 * par[2] / par[1])**(1 / 2.1)
    return par, result.fun, peak


def run(parent, matches, out):
    frame = pd.read_csv(parent)
    for prefix in PREFIX.values():
        ref = pd.read_csv(matches / f"match_{prefix}.csv.gz", index_col="row")
        if not ref.index.equals(frame.index):
            raise ValueError("Input row order differs from LoLSS selection")
        frame = frame.join(ref)
    frame["alpha_hi"] = np.log(frame.nvss_flux_jy / frame.lotss_flux_jy) / np.log(1400 / 144)
    frame["pred23"] = frame.lotss_flux_jy * (23 / 144)**frame.alpha_hi
    frame["response23"] = frame.psf_jy / frame.pred23
    frame["response54"] = (frame.Total_flux / 1000) / (frame.lotss_flux_jy * (54 / 144)**frame.alpha_hi)
    cols = ["lolrow","Source_name","RA","DEC","field","alpha_hi","pred23","response23","response54",
            "lotss_name","lotss_flux_jy","lotss_sum75_jy","lotss_sep_as","lotss_dc_maj",
            "nvss_sep_as","nvss_flux_jy","vlssr_flux_jy","wenss_flux_jy","gmrtas150m_flux_jy"]
    frame[cols].to_csv(out / "parent_comparison.csv", index=False, float_format="%.9g")
    shortlist = frame[frame.slope / frame.eslope > 3].copy()
    phot = pd.DataFrame([point for _, r in shortlist.iterrows() for point in photometry(r)],
                        columns=["source","survey","freq_mhz","flux_jy","stat_jy","frac_scale","limit95_jy"])
    phot.to_csv(out / "photometry.csv", index=False, float_format="%.9g")
    results = []
    for name, group in phot.groupby("source"):
        for selection in ("all","independent_of_23","high_freq"):
            data = group if selection == "all" else group[group.freq_mhz >= (54 if selection == "independent_of_23" else 144)]
            for model in ("power","curved","absorbed"):
                if selection == "high_freq" and model != "power":
                    continue
                if len(data) <= (2 if model=="power" else 3):
                    continue
                par, stat, peak = fit(data, model)
                pred23 = spectrum(par, np.log(23 / 144), model)
                limit = group.loc[group.survey == "LoDeSS", "limit95_jy"].iloc[0]
                # For a threshold-censored observation the likelihood factor is Phi((L-mu)/sigma).
                # The fits themselves use the more informative signed forced flux.
                noise23 = group.loc[group.survey == "LoDeSS", "stat_jy"].iloc[0]
                logp = log_ndtr((limit-pred23)/np.hypot(noise23,.39*pred23)) if np.isfinite(limit) else np.nan
                results.append([name,selection,model,len(data),*par[:2],par[2] if len(par)==3 else np.nan,
                                stat,stat+2*len(par),peak,pred23,logp])
    fits = pd.DataFrame(results, columns=["source","selection","model","n","log_A","alpha","shape",
                                          "minus2loglike_const","aic_const","peak_mhz","pred23_jy","logp_below_reported_limit"])
    fits.to_csv(out / "fits.csv", index=False, float_format="%.9g")
    name = "LOL1J112022.1+502936"
    lead = frame[frame.Source_name==name].iloc[0]
    distance = np.degrees(2*np.arcsin(np.sqrt(
        np.sin(np.radians(frame.DEC-lead.DEC)/2)**2 + np.cos(np.radians(frame.DEC))*np.cos(np.radians(lead.DEC))*
        np.sin(np.radians(frame.RA-lead.RA)/2)**2)))
    peers=frame[(frame.Source_name!=name)&(frame.field==lead.field)&(distance<=2)&
                (frame.lotss_flux_jy/frame.lotss_sum75_jy>.9)&(frame.lotss_sep_as<10)&
                (frame.nvss_sep_as<12)&(frame.lotss_dc_maj<12)].copy()
    peers[["Source_name","RA","DEC","response23","response54","slope"]].to_csv(out/"local_peers.csv",index=False)
    fig, ax = plt.subplots(figsize=(7.4,4.8),layout="constrained")
    g = phot[phot.source==name]
    ax.errorbar(g.freq_mhz,g.flux_jy,yerr=np.hypot(g.stat_jy,g.frac_scale*g.flux_jy),fmt='o',
                color='#124a6b',label='Independent survey measurements',capsize=3)
    freq = np.arange(44,65,4)
    ax.errorbar(freq,[lead[f's{f}']/1000 for f in freq],yerr=[lead[f'e{f}']/1000 for f in freq],
                fmt='s',color='#b35631',markersize=4,label='LoLSS subbands (fit errors)',capsize=2)
    anchors = fits[(fits.source==name)&(fits.model=='power')&(fits.selection=='high_freq')].iloc[0]
    grid=np.geomspace(23,1400,200)
    ax.plot(grid,spectrum([anchors.log_A,anchors.alpha],np.log(grid/144),'power'),ls='--',color='#777',
            label='Power law fitted at 144–1400 MHz')
    curved=fits[(fits.source==name)&(fits.model=='curved')&(fits.selection=='independent_of_23')].iloc[0]
    ax.plot(grid,spectrum([curved.log_A,curved.alpha,curved['shape']],np.log(grid/144),'curved'),
            color='#208068',label='Curved fit with 23 MHz held out')
    ax.set(xscale='log',yscale='log',xlabel='Observed frequency (MHz)',ylabel='Flux density (Jy)',
           title='NVSS J112022+502936')
    ax.legend(fontsize=8);ax.grid(alpha=.16);fig.savefig(out/"J112022-spectrum.png",dpi=180);plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["parent", "matches", "output"]:
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    run(args.parent, args.matches, args.output)
