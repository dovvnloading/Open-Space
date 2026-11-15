"""Plot derived measurements. Run after replay.py."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

HERE = Path(__file__).resolve().parent
with (HERE / "profiles.csv").open() as fh:
    data = [r for r in csv.DictReader(fh) if r["status"] == "measured"]
with (HERE / "simultaneous_quarters.csv").open() as fh:
    quarters = list(csv.DictReader(fh))
plt.rcParams.update({"svg.fonttype":"none","svg.hashsalt":"gleam-j1627","axes.spines.top":False,
                     "axes.spines.right":False,"font.size":10})
fig,(ax,sub)=plt.subplots(1,2,figsize=(9.5,3.85),gridspec_kw={"width_ratios":[1.05,1]})
band_ids={r["obsid"] for r in quarters}
for label,marker,color,predicate in [
    ("13 March, five published spectra","o","#006783",lambda r:r["obsid"] in band_ids),
    ("13 March, six other rotations","s","#c46124",lambda r:r["utc_start"].startswith("2018-03-13") and r["obsid"] not in band_ids),
    ("11, 15 and 16 March","x","#44444c",lambda r:not r["utc_start"].startswith("2018-03-13"))]:
    selected=[r for r in data if predicate(r)]
    ax.scatter([float(r["MHz"]) for r in selected],[float(r["spacing_s"]) for r in selected],
               marker=marker,color=color,s=35,label=label,zorder=3,alpha=.85)
ax.set(title="Spacing across separate rotations",xlabel="Observing band center (MHz)",
       ylabel="Fitted component spacing (s)",ylim=(5.9,9.5))
ax.legend(loc="upper right",frameon=False,fontsize=8.3)
for j,(band,color) in enumerate(zip([85,119,154,185,215],["#006783","#6a4692","#c46124","#36874f","#5f6672"])):
    rs=[r for r in quarters if int(r["band_file"])==band]
    sub.plot([float(r["MHz"]) for r in rs],[float(r["spacing_s"]) for r in rs],
             color=color,marker="o",markersize=4,lw=1.6)
sub.set(title="Simultaneous channels within each rotation",
        xlabel="Channel-group center (MHz)",ylim=(5.9,9.5))
for a in (ax,sub):
    a.grid(axis="y",lw=.5,color="#cccccc")
    a.xaxis.set_major_locator(FixedLocator([88,119,155,185,216]))
fig.tight_layout(w_pad=2)
fig.savefig(HERE / "spacing.svg",metadata={"Date":None})
