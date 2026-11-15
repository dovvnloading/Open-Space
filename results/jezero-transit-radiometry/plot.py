"""Plot the eight simultaneous calibrated photodiode profiles on the two highlighted sols."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from analysis import AXIS, CHANNELS
here=Path(__file__).resolve().parent
raw=pd.read_csv(here/"selected_windows.csv")
fig,axes=plt.subplots(2,1,figsize=(8.4,6.2),sharex=True,layout="constrained")
for ax,(sol,pair) in zip(axes,[(403,[6,7]),(1062,[2,3])]):
    block=raw[raw.sol==sol][CHANNELS].to_numpy()
    flank=(abs(AXIS)>=45)&np.isfinite(block).all(axis=1)
    m,b=np.polyfit(AXIS[flank],block[flank],1)
    values=100*(block/(AXIS[:,None]*m+b)-1)
    for i,name in enumerate(CHANNELS):
        key=int(name[-1])
        emphasized=(i==0 or key in pair)
        color="#12181b" if i==0 else ("#bb4d24" if key==pair[0] else "#057e86") if key in pair else "#a3afa9"
        name="Hemispheric TOP7" if i==0 else f"LAT{key}"
        smooth=pd.Series(values[:,i]).rolling(5,center=True).mean()
        ax.plot(AXIS,smooth,color=color,lw=1.9 if emphasized else 0.9,alpha=1 if emphasized else 0.65,
                zorder=3 if emphasized else 1,label=name if emphasized else None)
    ax.axvline(0,color="#9e9e9e",ls="--",lw=0.9)
    ax.set(xlim=(-55,55),ylim=(-2.2,0.55),ylabel="Flux change (%)",title=f"Sol {sol}  ·  {'8 April 2022' if sol==403 else '14 February 2024'}")
    ax.legend(loc="lower left",ncol=3,fontsize=8,frameon=False)
    ax.spines[["top","right"]].set_visible(False)
axes[-1].set_xlabel("Seconds from the predicted closest approach of Phobos to the Sun")
fig.savefig(here/"radiometry.png",dpi=155)
