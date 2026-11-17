"""Plot the fitted Keplerians and signed RV residuals from archived fits."""
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fit_controls import HERE, rv_model


def main():
    fit = json.loads((HERE/"control_fits.json").read_text())
    with (HERE/"control_residuals.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    fig = plt.figure(layout="constrained", figsize=(9,6.6))
    panels = fig.subfigures(2,2)
    for name, figpart in zip(["BH1","BH2","NS1","BH3"],panels.ravel()):
        ax, resid = figpart.subplots(2,1,sharex=True,height_ratios=[3,1])
        x = fit[name]["best_fit"]
        phase = np.linspace(0,1,301)
        time = x[2]+phase*x[0]
        ax.plot(phase,rv_model(x,time),color="#23827e",lw=1.4)
        phase, v, sig, dv = np.array([
            [((float(r["jd_rel_reference"])-x[2])/x[0]) % 1,float(r["rv_kms"]),
             float(r["err_kms"]),float(r["observed_minus_model_kms"])]
            for r in rows if r["control"]==name]).T
        for panel, y in [(ax,v),(resid,dv)]:
            panel.errorbar(phase,y,yerr=sig,fmt=".",ms=3,color="#213851",
                           ecolor="#7d93a6",elinewidth=.7)
            panel.spines[["top","right"]].set_visible(False)
        resid.axhline(0,color="#23827e",lw=.8)
        ax.set_title("Gaia "+name+("  ·  pre-DR4" if name=="BH3" else ""),loc="left")
        ax.set_ylabel("RV (km/s)")
        resid.set(xlabel="Orbital phase from periastron",ylabel="O−C\n(km/s)",xlim=(0,1))
    plt.savefig(HERE/"control_orbits.png",dpi=140)
if __name__=="__main__":
    main()
