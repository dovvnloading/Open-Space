"""Draw the two observing seasons from the deposited measurements."""
import csv
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = Path(__file__).resolve().parent
EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)
with (ROOT / "lsxps.csv").open() as src:
    dates = {r["obsid"]: EPOCH + timedelta(days=float(r["mjd"])) for r in csv.DictReader(src) if r["band"] == "Total"}
packet = json.loads((ROOT / "swift-products.json").read_text())
result = json.loads((ROOT / "result.json").read_text())
fig, axes = plt.subplots(1,2,figsize=(9,3.9),sharey=True)
for ax,year,color in zip(axes,(2025,2026),("#355b83","#a3432f")):
    for product in packet["products"]:
        stamp = dates[product["obs"]]
        if stamp.year != year or product["mode"] != "UVOT":
            continue
        for row in product["rows"]:
            if row["band"]=="UVW1" and row["flags"]==0 and row["sep"]<=2 and row["magerr"]<.2:
                ax.errorbar(stamp,row["mag"],yerr=row["magerr"],fmt="o",ms=4,color=color,capsize=2)
    bounds = result["uvot"][str(year)]["UVW1"]
    ax.axhline(bounds["median_mag"],color=color,ls="--",lw=.8)
    ax.set_title(f"{year} observing season",loc="left",fontsize=11)
    ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3,maxticks=5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(axis="y",alpha=.18)
    ax.spines[["top","right"]].set_visible(False)
axes[0].invert_yaxis()
axes[0].set_ylabel("UVW1 magnitude (brighter upward)")
fig.suptitle("AT 2021aeuk was ultraviolet-bright during the X-ray low state",x=.07,ha="left",fontsize=14)
plt.subplots_adjust(left=.09,right=.98,top=.78,bottom=.20,wspace=.14)
fig.text(.09,.04,"One point per Swift/UVOT visit. Error bars show catalogue uncertainties. Dashed lines mark season medians.",fontsize=8)
plt.savefig(ROOT / "uv-seasons.png",dpi=160)
