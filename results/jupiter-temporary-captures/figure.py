"""Plot the projected Jovian paths and the near-term astrometric tests."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from replay import AU, GM, load_tracks

ROOT = Path(__file__).resolve().parent
tracks = load_tracks()
fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.9), layout="constrained")
for ax, key, title in zip(axes, ["uq146-past","uq146-future","oe31-future"],
                          ["2025 UQ146  ·  1984–87","2025 UQ146  ·  2069–71",
                           "2019 OE31  ·  2094–2107"]):
    states = tracks[key]
    coords = states[:,1:4]/AU
    radius = np.linalg.norm(coords,axis=1)
    energy = np.sum(states[:,4:]**2,axis=1)/2-GM/(AU*radius)
    selected = energy < 0
    xs,ys = coords[:,:2].T
    ax.plot(xs,ys,color="#b7bdc4",lw=1)
    ax.plot(xs[selected],ys[selected],color="#185779",lw=1.9)
    ax.scatter([0],[0],s=38,color="#bf653d",zorder=4,label="Jupiter")
    ax.scatter(xs[selected][0],ys[selected][0],s=35,facecolor="white",
               edgecolor="#185779",zorder=3)
    ax.set_title(title,fontsize=11)
    ax.set(xlim=(-.75,.75),ylim=(-.75,.75),aspect="equal",xlabel="Jupiter-centered x (au)")
    ax.set_xticks([-.5,0,.5]); ax.set_yticks([-.5,0,.5])
    ax.grid(alpha=.17); ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("Jupiter-centered y (au)")
fig.suptitle("Projected paths from JPL Horizons",fontsize=14)
fig.savefig(ROOT/"paths.png",dpi=140)
