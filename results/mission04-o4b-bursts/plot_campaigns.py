"""Draw the time-shift ranks and real-data detection efficiencies for both surveys."""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).parent/"results")
    directory = parser.parse_args().results
    first = json.loads((directory/"survey.json").read_text())
    second = json.loads((directory/"holdout"/"campaign.json").read_text())
    sensitivity = json.loads((directory/"sensitivity.json").read_text())["cases"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.3), gridspec_kw={"width_ratios":[1.25, 1]},
                               layout="constrained")
    colors = ("#427085", "#b1482e")
    for survey, name, color in ((first, "10 frames", colors[0]), (second,"125 frames",colors[1])):
        bg = np.sort([r[-1] for r in survey["campaign_slides"] if r[0]])
        observed = survey["foreground_loudest"]
        unique, counts = np.unique(bg,return_counts=True)
        survival = np.cumsum(counts[::-1])[::-1] / len(bg)
        a.step(np.r_[1.,unique],np.r_[1.,survival],where="pre",color=color,label=name)
        rank = np.mean(bg>=observed)
        a.scatter(observed, max(rank,0.5/len(bg)),s=50,color=color,zorder=5,edgecolors="white",linewidth=.6)
    a.axvline(20.,color="#555",linestyle=":",lw=1)
    a.text(21,.18,"Rate threshold",fontsize=8,color="#464646")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(left=2.,right=max(30,1.12*max(first["campaign_slides"][i][-1] for i in range(1,129)),
                                    1.12*max(second["campaign_slides"][i][-1] for i in range(1,129))))
    a.set_ylim(.004,1.2)
    a.set(xlabel="Maximum coherent score per campaign",ylabel="Fraction of shifted campaigns reaching score",
          title="Observed maxima and time shifts")
    a.set_yticks([.01,.05,.2,.5,1],labels=["1%","5%","20%","50%","100%"])
    a.legend(title="Dots mark zero lag",fontsize=8,title_fontsize=8,loc="lower left")
    fields = ["lower","middle","primary","higher_frequency"]
    names = ["150 Hz   0.2","150 Hz   0.5","150 Hz   1.0","350 Hz   1.0","Holdout   150 Hz   1.0"]
    for i, name in enumerate(fields+["holdout"]):
        if name=="holdout":
            row = second["rate"]
            k, n = row["recovered"], row["draws"]
        else:
            row = sensitivity[name]
            k,n = row["recovered"],row["n"]
        p = k/n
        lo = beta.ppf(.025,k,n+1-k) if k else 0.
        hi = beta.ppf(.975,k+1,n-k) if k<n else 1.
        b.errorbar(p,i,xerr=[[p-lo],[hi-p]],color=colors[int(name=="holdout")],fmt="o",capsize=3,ms=6)
        b.text(p+.045,i,f"{100*p:.1f}%",va="center",fontsize=8)
    b.axhline(3.5,color="#aaa",lw=.8,ls=":")
    b.set_yticks(range(5),labels=names)
    b.invert_yaxis()
    b.set_xlim(.5,1.)
    b.set_xticks([.5,.6,.7,.8,.9,1],labels=["50%","60%","70%","80%","90%","100%"])
    b.set(xlabel="Fraction recovered (95% binomial interval)", title="Injections into observed strain")
    b.set_ylabel("Frequency and incident strain\n" + r"($10^{-21}$ Hz$^{-1/2}$)")
    for ax in (a,b):
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="x" if ax is b else "y",alpha=.2)
        ax.set_axisbelow(True)
    fig.savefig(directory/"campaigns.png",dpi=160)
    plt.close(fig)


if __name__=="__main__":
    main()
