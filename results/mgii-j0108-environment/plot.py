"""Plot projected positions and independent archival spectra."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def rows(filename):
    with (HERE / filename).open() as handle:
        return list(csv.DictReader(handle))


absorbers = rows("j0108-absorbers.csv")
galaxies = [r for r in rows("foreground-galaxies.csv") if float(r["center_arcmin"]) <= 6]
cluster, = [r for r in rows("cluster-crossmatches.csv") if r["wen_han_id"] == "78385"]
fig = plt.figure(figsize=(9, 4.65))
grid = fig.add_gridspec(3, 2, width_ratios=[1.03, 1.4], wspace=.38, hspace=.16)
sky = fig.add_subplot(grid[:, 0])


def offset(ra, dec):
    return (float(ra)-16.999238)*60*np.cos(np.deg2rad(-1.130688)), (float(dec)+1.130688)*60


for i, g in enumerate(galaxies):
    sky.scatter(*offset(g["ra"], g["dec"]), marker="s", s=28, c="#818c7b",
                label="Galaxy within ±500 km/s" if i == 0 else None)
sky.scatter(*offset(cluster["ra_deg"],cluster["dec_deg"]), marker="*", s=125, color="#ab782d",
            label="WH cluster center", zorder=3)
for letter, a in zip("ABC",absorbers):
    px,py = offset(a["RA"],a["DEC"])
    sky.scatter(px,py,facecolors="none",edgecolors="#2e4c79",s=78,lw=1.8)
    sky.annotate(letter,(px,py),xytext=(6,3),textcoords="offset points",weight="bold",color="#2e4c79")
sky.set(xlim=(4.1,-4.1), ylim=(-6.8,5.8), xlabel="RA offset (arcmin)", ylabel="Dec offset (arcmin)", aspect="equal")
sky.grid(alpha=.16)
sky.legend(loc="upper left",fontsize=8,frameon=False)
sky.set_title("Foreground field",loc="left",weight="bold")
for i,(a,fiber) in enumerate(zip(absorbers,[228,230,232])):
    ax=fig.add_subplot(grid[i,1])
    z=float(a["Z_MGII"])
    for s,c,label in [(f'spec-4226-55475-{fiber:04}',"#ca7945","SDSS 2010"),("DESI-"+a["TARGETID"],"#254c7d","DESI 2021")]:
        d=np.genfromtxt(HERE/"spectra"/"pixels"/(s+".csv"),delimiter=",",names=True,dtype=None,encoding="utf-8")
        w=d["wave_A"]/(1+z)
        x=w[(w>2781)&(w<2820)]
        keep=(w>2781)&(w<2820)
        y=d["flux"][keep]; iv=d["ivar"][keep]
        side=(abs(x-2796.3543)>4.1)&(abs(x-2803.5315)>4.1)&(iv>0)&(d["mask"][keep]==0)
        continuum=np.polyval(np.polyfit(x[side]-2800,y[side],1,w=np.sqrt(iv[side])),x-2800)
        ax.step(x,np.where((iv>0)&(d["mask"][keep]==0),y/continuum,np.nan),where="mid",color=c,lw=1.15,label=label)
    ax.axvline(2796.3543,ls=":",lw=.8,color=".45");ax.axvline(2803.5315,ls=":",lw=.8,color=".45")
    ax.set(xlim=(2789,2812),ylim=(-.15,1.65),yticks=[0,.5,1,1.5])
    ax.text(.025,.12,"ABC"[i],transform=ax.transAxes,weight="bold")
    if i<2:ax.tick_params(labelbottom=False)
    else:ax.set_xlabel("Absorber rest wavelength (Å)")
    if i==0:
        ax.set_title("Mg II doublet in two surveys",loc="left",weight="bold")
        ax.legend(loc="upper right",frameon=False,ncol=2,fontsize=8)
    if i==1:ax.set_ylabel("Normalized flux")
fig.subplots_adjust(left=.09,right=.98,bottom=.14,top=.87)
fig.savefig(HERE/"field-and-spectra.png",dpi=135)
