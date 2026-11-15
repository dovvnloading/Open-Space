"""Draw the direct-light predictions from the archived inputs."""
from datetime import datetime
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from replay import calculate


def plot(output):
    vg = calculate()
    ca = calculate(cassini=True)
    t = [datetime.fromisoformat(r[0]) for r in vg]
    fig, ax = plt.subplots(figsize=(10, 4.7), constrained_layout=True)
    ax.plot(t, [r[2] for r in vg], color="#707478", ls="--", lw=1.7, label="Saturn's globe alone")
    ax.plot(t, [r[1] for r in ca], color="#b95720", lw=2, label="Globe and rings, Cassini VIMS profile")
    ax.plot(t, [r[1] for r in vg], color="#14577a", lw=1.4, label="Globe and rings, Voyager 2 PPS profile")
    ax.set_ylim(-.045, 1.10)
    ax.set_xlim(datetime(2032, 3, 1, 10), datetime(2032, 3, 3, 21))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 6, 12, 18]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M"))
    ax.set_xlabel("Earth reception time (UTC), 2032")
    ax.set_ylabel("Fraction of direct sunlight reaching Phoebe")
    ax.set_title("Saturn’s main rings cast a separate deep shadow on each side of Phoebe’s eclipse", loc="left", fontsize=12)
    ax.legend(loc="center", bbox_to_anchor=(.54,.53), frameon=True, fontsize=9)
    ax.grid(axis="y", color=".92")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.savefig(output, dpi=165)


if __name__ == "__main__":
    plot(Path(__file__).with_name("illumination.png"))
