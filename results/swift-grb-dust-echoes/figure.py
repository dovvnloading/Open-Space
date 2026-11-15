"""Draw the GRB 180728A time–radius track and radial profile."""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from search import ROOT

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("photons")
p.add_argument("output")
args=p.parse_args()
t,x,y,e=np.load(args.photons)
soft=(e>=70)&(e<=300)
r=np.hypot(x,y)
fig,(ax,bx)=plt.subplots(1,2,figsize=(9.8,4.2),layout="constrained",
                         gridspec_kw={"width_ratios":[1.3,1]})
# Plot radial counts per unit solid angle in each occupied snapshot; normalize each column
# to its own exposure so the interrupted observing schedule does not erase late visits.
times=np.sort(t)
windows=np.split(times,np.flatnonzero(np.diff(times)>500)+1)
edges=np.arange(120,551,8)
for w in windows:
    if w.mean()/3600>14:continue
    keep=soft&(t>=w[0])&(t<=w[-1])
    h=np.histogram(r[keep],edges)[0]/np.diff(edges**2)
    scale=np.percentile(h,90)
    if scale>0:
        xmid=w.mean()/3600
        ax.pcolormesh([xmid-.36,xmid+.36],edges,(h/scale)[:,None],
                      cmap="Blues",vmin=0,vmax=1.65,shading="flat")
time=np.linspace(1.1,14,200)
ax.plot(time,np.sqrt(826.6*3600*time/157.3),color="#b4491f",lw=1.7)
ax.text(.04,.96,"157 pc prediction",transform=ax.transAxes,color="#9d3c1a",va="top")
ax.set(xlim=(1.1,14),ylim=(120,550),xlabel="Time since burst (hours)",ylabel="Radius (arcsec)",
       title="The excess moves outward")
data=np.genfromtxt(ROOT/"profile.csv",delimiter=",",names=True)
bx.step(data["u"],data["photons"],where="mid",color="#244965",lw=1.4,label="XRT events")
bx.plot(data["u"],data["background"],"--",color="#777",label="Smooth background")
bx.plot(data["u"],data["background"]+data["line"],color="#b4491f",label="Background + ring")
bx.set(xlabel="Squared radius / predicted squared radius",ylabel="Events per bin",
       title="The expanding track aligns at unity")
bx.legend(frameon=False,fontsize=8)
for a in (ax,bx):
    a.spines[["top","right"]].set_visible(False)
fig.savefig(args.output,dpi=145)
