"""Plot source-time L1 measurements and lunar/OMNI measurements on separate UTC axes."""
from pathlib import Path
from io import BytesIO

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

root = Path(__file__).resolve().parent
data = pd.read_csv(root / "minute.csv.gz", parse_dates=["utc"], index_col="utc")
panels = [
    ("a   Solar wind monitors near L1", "17:30", "18:30", [
        ("ace_bz", "ACE", "#00816d"), ("dscovr_bz", "DSCOVR", "#1268aa"), ("wind_bz", "Wind", "#a14812")]),
    ("b   Near the Moon and at the modeled bow-shock nose", "18:30", "19:30", [
        ("thb_bz", "ARTEMIS P1", "#00816d"), ("thc_bz", "ARTEMIS P2", "#1268aa"), ("omni_bz", "OMNI", "#a14812")]),
]
fig, axes = plt.subplots(2, 1, figsize=(8.1, 5.8), sharex=True, sharey=True)
for ax, (name, start, stop, channels) in zip(axes, panels):
    ax.axvspan(pd.Timestamp("2019-01-31 " + start, tz="UTC"), pd.Timestamp("2019-01-31 " + stop, tz="UTC"),
               facecolor="#eef0e9", linewidth=0)
    ax.axhline(0, color="#616161", linewidth=0.7)
    for field, label, color in channels:
        ax.plot(data.index, data[field], label=label, color=color, linewidth=1)
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=8)
    ax.set_title(name, loc="left", fontsize=11, fontweight="bold")
    ax.set_ylabel("GSM Bz (nT)")
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_ylim(-24, 26)
ax.set_xlim(pd.Timestamp("2019-01-31 16:30",tz="UTC"), pd.Timestamp("2019-01-31 20:30",tz="UTC"))
ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0,30]))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_xlabel("31 January 2019 (UTC)")
fig.tight_layout()
buffer = BytesIO()
fig.savefig(buffer, format="png", dpi=100)
buffer.seek(0)
Image.open(buffer).convert("RGB").quantize(colors=64).save(root / "polarities.png", optimize=True)
