"""Draw the measured pulses and published 68% localization contours as an SVG."""

import argparse
from pathlib import Path

import h5py
import numpy as np

NAMES = ("FRB20200511C", "FRB20211128E")
COLORS = ("#136273", "#ad4b28")
PALETTE = ["#173451", "#365774", "#596e82", "#78838a", "#92948f",
           "#aca18d", "#c4ae78", "#d9bf67", "#eccc50", "#ffe335"]


def heatmap(path):
    with h5py.File(path, "r") as f:
        good = f["good_freq"][:].astype(bool)
        valid = f["flag"][:].astype(bool) & good[:, None]
        times = (f["index_map/times"][:] - f.attrs["center_time"]) * 1000
        data = np.ma.masked_array(f["data"][:], mask=~valid)
    off = np.abs(times) > 20
    data = (data - data[:, off].mean(axis=1)[:, None]) / data[:, off].std(axis=1, ddof=1)[:, None]
    rebinned = data.reshape(16, 1024, len(times)//2, 2).mean(axis=(1,3))
    centers = times.reshape(-1,2).mean(axis=1)
    return centers, rebinned / rebinned[:, np.abs(centers) > 20].std(axis=1, ddof=1)[:, None]


def draw(directory):
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="980" height="645" viewBox="0 0 980 645">',
           '<title>Measured CHIME burst intensities and localization regions</title>',
           '<rect width="980" height="645" fill="white"/>',
           '<g font-family="sans-serif" fill="#172633">']

    def label(x, y, text, size=15, anchor="middle", color="#172633", bold=False):
        weight = "bold" if bold else "normal"
        svg.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}"'
                   f' text-anchor="{anchor}" fill="{color}">{text}</text>')

    label(44,36,'Two CHIME bursts at a dispersion measure near 1,835 pc cm⁻³',25,"start",bold=True)
    for n, name in enumerate(NAMES):
        left, top, width, height = 65 + n*475, 97, 391, 178
        label(left + width/2,80,name.replace("FRB","FRB "),19,color=COLORS[n],bold=True)
        times, intensity = heatmap(directory / f"{name}_stokesi_dynamic_spectrum.h5")
        dt = times[1] - times[0]
        for t, center in enumerate(times):
            l, r = max(-30, center-dt/2), min(30, center+dt/2)
            if l >= r:
                continue
            for k,value in enumerate(intensity[:,t]):
                index = int(np.clip(np.floor((float(value)+2.5)/6.5*len(PALETTE)),0,len(PALETTE)-1))
                x, y = left+(l+30)*width/60, top+(15-k)*height/16
                svg.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{(r-l)*width/60+.05:.2f}"'
                           f' height="{height/16+.05:.2f}" fill="{PALETTE[index]}"/>')
        for tick in range(-30,31,10):
            label(left+(tick+30)*width/60,top+height+22,tick,12)
        for tick in (400,500,600,700,800):
            label(left-7,top+height-(tick-400)*height/400+4,tick,12,"end")
        label(left+width/2,top+height+44,"Time from CHIME fit (ms)",14)
    label(16,97,"MHz",12,"start")
    label(44,375,"CHIME’s 68% localization contours",21,"start",bold=True)
    left, top, w, h = 75, 410, 848, 130
    def sky(x,y):
        return left+(24-x)*w/8, top+(51.52-y)*h/.82
    for tick in range(16,25):
        x, _ = sky(tick,51.0)
        label(x,top+h+20,tick,12)
    for d in (50.8,51.0,51.2,51.4):
        _, y = sky(20,d)
        label(left-9,y+4,f"{d:.1f}°",12,"end")
    label(44,398,"Dec.",12,"start")
    svg.append(f'<rect x="{left}" y="{top}" width="{w}" height="{h}" fill="none" stroke="#9aa2a5"/>')
    for n,name in enumerate(NAMES):
        color=COLORS[n]
        with h5py.File(directory / f"{name}_localization.h5") as f:
            for g in f["contours/68"].values():
                points = " ".join(f"{x:.2f},{y:.2f}" for x,y in (sky(*p) for p in g[:].T))
                svg.append(f'<polygon points="{points}" fill="{color}" fill-opacity=".23"'
                           f' stroke="{color}" stroke-width="1.6"/>')
        label(345+290*n,616,name.replace("FRB","FRB "),16,color=color,bold=True)
    label(500,584,"Right ascension (degrees, ICRS)",14)
    svg.append("</g></svg>")
    return "\n".join(svg)+"\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing the four original HDF5 products")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(draw(args.input_dir))
