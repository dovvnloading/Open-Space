"""Draw the measured S Vir and V Del peaks and the post-gap comparison."""

from datetime import date
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from replay import HERE, load_rows


def main():
    peaks = load_rows("comparison-maxima.csv")
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.8), constrained_layout=True)
    for ax, (star, first, last, stop) in zip(
        axes, [("S Vir", -6, 18, 35), ("V Del", -1, 25, 36)]
    ):
        data = [
            r
            for r in peaks
            if r["star"] == star and first - 2 <= int(r["cycle"]) <= stop
        ]
        cycles = {int(r["cycle"]): r for r in data}
        for e, row in cycles.items():
            if e + 1 in cycles:
                nxt = cycles[e + 1]
                ax.plot(
                    [date.fromisoformat(row["date"]), date.fromisoformat(nxt["date"])],
                    [float(row["magnitude"]), float(nxt["magnitude"])],
                    c="#9ba9ac",
                    lw=1,
                    zorder=1,
                )
        ax.axvspan(
            date.fromisoformat(cycles[first]["date"]),
            date.fromisoformat(cycles[last]["date"]),
            color="#c5dfd9",
            alpha=0.45,
            lw=0,
        )
        for parity, color, label in [
            (0, "#185b68", "Even cycle"),
            (1, "#ab4b29", "Odd cycle"),
        ]:
            band = [r for r in data if int(r["cycle"]) % 2 == parity]
            ax.scatter(
                [date.fromisoformat(r["date"]) for r in band],
                [float(r["magnitude"]) for r in band],
                color=color,
                s=22,
                label=label,
                zorder=2,
            )
        ax.invert_yaxis()
        ax.set_ylabel("Peak visual magnitude")
        ax.set_title(star, loc="left", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.2)
    axes[0].legend(loc="upper right", frameon=False, ncols=2, fontsize=9)
    axes[1].set_xlabel("Year of maximum")
    fig.savefig(HERE / "alternating-maxima.png", dpi=150)


if __name__ == "__main__":
    main()
