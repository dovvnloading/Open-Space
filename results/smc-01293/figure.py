"""Regenerate the phase profile and independent-survey projections."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from replay import ROOT, EPOCH, get_photometry, fit_periods, clean, waveforms, season, project

data = get_photometry()
train = clean(data["OGLE-IV-I"])
periods, coeffs, _ = fit_periods(train)
fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), constrained_layout=True)
ax = axes[0]
for name,color in (("OGLE-IV-I", "#15447a"),("EROS-sm0010l7230-B","#ab4d25")):
    d = clean(data[name])
    t,m,e = d.T
    a,b = waveforms(t,periods,coeffs)
    result = project(data[name],periods,coeffs)
    res = m - result["primary_scale"]*a
    for y in np.unique(season(t)):
        k = season(t)==y
        res[k] -= np.median(res[k]-result["secondary_scale"]*b[k])
    ph = ((t-EPOCH)/periods[1])%1
    xx=[]; yy=[]; se=[]
    for j in range(32):
        mask=(ph>=j/32)&(ph<(j+1)/32)
        samples=res[mask]
        if samples.size<2:continue
        xx.append((j+.5)/32); yy.append(np.median(samples)*1000)
        se.append(1.86*np.median(abs(samples-np.median(samples)))/np.sqrt(samples.size)*1000)
    label="OGLE-IV I" if name=="OGLE-IV-I" else "EROS-2 blue"
    ax.errorbar(xx,yy,yerr=se,fmt="o",ms=3.5,c=color,lw=.8,label=label,alpha=.85)
x=np.linspace(0,1,501)
ax.plot(x,1000*waveforms(EPOCH+x*periods[1],periods,coeffs)[1],color="#111111",lw=1.5,label="OGLE-IV fit")
ax.invert_yaxis();ax.set(xlabel="Phase at 0.529917 days",ylabel="Residual magnitude (mmag)",title="Second signal after subtracting the binary")
ax.legend(loc="lower right",fontsize=8,frameon=False)
s=json.loads((ROOT/"results.json").read_text())
rows=[("OGLE-II I",s["survey_projections"]["OGLE-II-I"]),("OGLE-III I (disjoint)",s["ogle_iii_i_excluding_ii_epochs"]),
      ("OGLE-IV V",s["survey_projections"]["OGLE-IV-V"]),("EROS-2 red",s["survey_projections"]["EROS-sm0010l7230-R"]),
      ("EROS-2 blue",s["survey_projections"]["EROS-sm0010l7230-B"]),("DECam griz",s["decam"]["all_bands"])]
ax=axes[1];y=np.arange(len(rows));ax.axvline(0,color="#777",lw=.8);ax.axvline(1,color="#777",lw=.8,ls="--")
ax.errorbar([v["secondary_scale"] for _,v in rows],y,xerr=[1.96*v["secondary_se"] for _,v in rows],fmt="o",c="#155a52",capsize=3)
ax.set(yticks=y,yticklabels=[k for k,_ in rows],xlabel="Amplitude relative to OGLE-IV I",title="Fixed-period, fixed-phase tests",xlim=(-.08,1.45))
ax.invert_yaxis();fig.savefig(ROOT/"secondary-signal.png",dpi=170)
