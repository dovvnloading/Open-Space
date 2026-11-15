"""Plot native radio-power changes across adjoining seven-second antenna commands."""
import argparse
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from bias_check import HERE, BANDS, AEB_URL, AEB_SHA, seconds, native_power, get_cdf


def main():
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aeb", type=Path, required=True)
    parser.add_argument("--hfr", type=Path, required=True)
    args = parser.parse_args()
    event = next(e for e in json.loads((HERE/"manifest.json").read_text()) if e["date"] == "2024-11-06")
    aeb = get_cdf(args.aeb,AEB_URL,AEB_SHA)
    hfr = get_cdf(args.hfr,event["cdf"],event["sha256"])
    bt = seconds(aeb,"epoch_aeb2")
    current = aeb.varget("psp_fld_l2_aeb2_V3_WHIP_CURR_DAC").astype(float)
    ts,power = native_power(hfr,"ch1_V3V4")
    from bias_check import step_differences
    dc,dp = step_differences(bt,current,ts,power,0)
    fig,axs=plt.subplots(1,2,figsize=(8.4,3.65),layout="constrained",sharex=True)
    for ax,bb in zip(axs,[0,2]):
        ax.axhline(0,color="#b7bbbc",lw=.8)
        ax.axvline(0,color="#b7bbbc",lw=.8)
        ax.scatter(dc,dp[:,bb],s=22,color="#176e78",alpha=.7,edgecolors="none")
        line=np.polyfit(dc,dp[:,bb],1)
        bounds=np.array([dc.min(),dc.max()])
        ax.plot(bounds,np.polyval(line,bounds),color="#8f365b",lw=1.6)
        ax.set_title(f"{BANDS[bb][0]}–{BANDS[bb][1]} MHz",loc="left",fontsize=11,weight="bold")
        ax.set_xlabel("Change in commanded current (µA)")
        for side in ["top","right"]:ax.spines[side].set_visible(False)
    axs[0].set_ylabel("Change in V3–V4 power (dB)")
    fig.suptitle("Parker radio response between adjacent antenna settings",weight="bold",fontsize=14)
    fig.savefig(HERE/"bias-response.png",dpi=135)
    plt.close(fig)


if __name__=="__main__":
    main()
