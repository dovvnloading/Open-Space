"""Plot the measured orbit and the sensitivity of the SB1 minimum mass to M1."""
import csv
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fit_controls import HERE
from sdss_sb1 import FIELDS,curve

def main():
    fit=json.loads((HERE/"giant_sb1_fits.json").read_text())["results"]["nominal"]
    p=np.array([fit["best"][k] for k in FIELDS])
    with (HERE/"giant_sb1_residuals.csv").open() as src:
        rows=[r for r in csv.DictReader(src) if r["model"]=="nominal"]
    ts,rv,err,res=np.array([[float(r[k]) for k in ["jd_rel_DR3","rv_kms","err_kms","residual_kms"]]
                            for r in rows]).T
    phase=((ts-p[5])/p[0])%1
    grid=np.linspace(0,1,301)
    fig=plt.figure(figsize=(9,4.3),layout="constrained")
    layout=fig.add_gridspec(2,2,height_ratios=[3,1])
    ax=fig.add_subplot(layout[0,0])
    oc=fig.add_subplot(layout[1,0],sharex=ax)
    mass=fig.add_subplot(layout[:,1])
    ax.plot(grid,curve(p,p[5]+p[0]*grid).ravel()+fit["best"]["lamost_offset_kms"],color="#126663")
    for a,y in [(ax,rv),(oc,res)]:
        a.errorbar(phase,y,yerr=err,fmt="o",ms=3,color="#223e5e",elinewidth=.8)
        a.set_xlim(0,1)
    ax.set_title("UCAC4 720-033289",loc="left")
    ax.set_ylabel("Velocity (km/s)")
    ax.tick_params(labelbottom=False)
    oc.axhline(0,lw=1,color="#126663")
    oc.set(xlabel="Orbital phase from periastron",ylabel="O−C\n(km/s)")
    for tag,lbl,col in [("catalog","Catalog uncertainty","#126663"),
                         ("30pct_uncertainty","30% primary-mass uncertainty","#b35b35"),
                         ("30pct_lower_primary","Primary mass reduced 30%","#506f97")]:
        pdf=fit["mass_models"][tag]["density"]
        mass.stairs(pdf["pdf"],pdf["edges"],color=col,linewidth=1.4,label=lbl)
    mass.axvline(3,color="0.35",ls=":",lw=1)
    mass.legend(frameon=False,fontsize=7.5)
    mass.set(title="Sensitivity to primary mass",xlabel="Edge-on companion mass (solar masses)",
             ylabel="Posterior density")
    for a in [ax,oc,mass]:a.spines[["top","right"]].set_visible(False)
    fig.savefig(HERE/"ucac4_720_033289_orbit.png",dpi=140)
if __name__=="__main__":
    main()
