"""Phase curve, signed nightly residuals, and companion-mass likelihood sensitivity."""
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fit_controls import HERE,rv_model


def main():
    all_fits=json.loads((HERE/"lamost_candidate_fits.json").read_text())["models"]
    label,nominal=next(iter(all_fits.items()))
    with (HERE/"lamost_candidate_residuals.csv").open() as f:
        rows=[r for r in csv.DictReader(f) if r["fit"]==label]
    t,v,e,r=np.array([[float(row[k]) for k in ["jd_rel_reference","rv_kms","err_kms","observed_minus_model_kms"]]
                     for row in rows]).T
    x=nominal["best_fit"]
    grid=np.linspace(0,1,201)
    phase=((t-x[2])/x[0]) % 1
    fig=plt.figure(figsize=(9,4.3),layout="constrained")
    layout=fig.add_gridspec(2,2,height_ratios=[3,1])
    ax=fig.add_subplot(layout[0,0])
    res=fig.add_subplot(layout[1,0],sharex=ax)
    mass=fig.add_subplot(layout[:,1])
    ax.plot(grid,rv_model(x,x[2]+grid*x[0]),color="#126663",label="Joint orbit")
    for p,y in [(ax,v),(res,r)]:
        p.errorbar(phase,y,yerr=e,fmt="o",ms=3,color="#223e5e",elinewidth=.8)
        p.set_xlim(0,1)
    ax.set_title("HD 264291",loc="left")
    ax.set_ylabel("Velocity (km/s)")
    ax.tick_params(labelbottom=False)
    res.set(xlabel="Orbital phase from periastron",ylabel="O−C\n(km/s)")
    res.axhline(0,color="#126663",lw=1)
    lbls=["Catalog covariance","Inflated errors","FLAME primary mass"]
    for info,col,lbl in zip(all_fits.values(),["#126663","#b35b35","#506f97"],lbls):
        p=info["mass_posterior_density"]
        mass.stairs(p["pdf"],p["edges"],color=col,label=lbl,linewidth=1.3)
    mass.set(xlabel="Companion mass (solar masses)",ylabel="Posterior density",
             title="Conditional on zero companion light")
    mass.legend(frameon=False,fontsize=8,loc="upper right")
    for a in [ax,res,mass]:a.spines[["top","right"]].set_visible(False)
    fig.savefig(HERE/"hd264291_orbit.png",dpi=140)
if __name__=="__main__":
    main()
