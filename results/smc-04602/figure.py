"""Plot survey timings and the two endpoint eclipse profiles."""
import json
from datetime import datetime,timedelta
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from analyze import BASE,PERIOD,EPOCH,ogle,wrap
from cadence import model

def year(jd):
    date=datetime(2000,1,1,12)+timedelta(days=jd-2451545)
    return date.year+(date-datetime(date.year,1,1)).total_seconds()/(86400*365.25)

def figure():
    results=json.loads((BASE/"results/timing.json").read_text())["measurements"]
    space=json.loads((BASE/"tess/tess_result.json").read_text())["combined"]
    colors={"OGLE":"#113F58","EROS":"#B6542A","TESS":"#267C70"}
    plt.rcParams.update({"font.size":9,"font.family":"DejaVu Sans","axes.spines.top":False,
                         "axes.spines.right":False,"axes.spines.left":False})
    fig,(ax,light)=plt.subplots(1,2,figsize=(10.2,3.5),gridspec_kw={"width_ratios":[1.4,1]})
    for a in (ax,light):
        a.grid(axis="y",color="#E2E8E9",linewidth=.65);a.set_axisbelow(True)
        a.spines["bottom"].set_color("#9DB0B2");a.tick_params(length=0,pad=5)
    g=[(k,v) for k,v in results.items() if k.startswith("OGLE")]
    for _,r in g:
        ax.plot([year(x) for x in r["hjd_range"]],[r["hours"]]*2,color=colors["OGLE"],lw=1,alpha=.6)
    def points(x,y,low,high,col,mark,label=None):
        ax.errorbar(x,y,yerr=[np.subtract(y,low),np.subtract(high,y)],fmt=mark,color=col,
                    ms=6,elinewidth=1.25,capsize=2,label=label,zorder=3)
    x=[year(r["median_hjd"]) for _,r in g]; y=[r["hours"] for _,r in g]
    points(x,y,[r["interval95"][0] for _,r in g],[r["interval95"][1] for _,r in g],
           colors["OGLE"],"o","OGLE I")
    e=[r for k,r in results.items() if k.startswith("sm0057l21161-B")]
    points([year(r["median_hjd"]) for r in e],[r["hours"] for r in e],
           [r["interval95"][0] for r in e],[r["interval95"][1] for r in e],
           colors["EROS"],"s","EROS-2 blue")
    ts=space["separation"]*PERIOD*24
    li,ui=np.array(space["separation95"])*PERIOD*24
    points([2025.55],[ts],[li],[ui],colors["TESS"],"D","TESS")
    ax.set(ylim=(17.6,26.5),xlim=(1996,2027),xticks=[2000,2005,2010,2015,2020,2025],
           yticks=[18,20,22,24,26],ylabel="Primary to secondary eclipse (hours)",xlabel="Observation year")
    ax.set_title("Eclipse spacing across three surveys",loc="left",fontweight="bold",pad=12)
    ax.legend(loc="lower right",frameon=False,fontsize=8,handletextpad=.2)
    for series,key,clr,label in [(ogle('ii'),'OGLE-II','#968579','1997–2000'),(ogle('iv'),'OGLE-IV-late',colors['OGLE'],'2022–2024')]:
        if key=='OGLE-IV-late': series=series[series[:,0]>=2459700]
        pars=np.array(results[key]["parameters"])
        ph=((series[:,0]-EPOCH)/PERIOD-pars[3])%1
        baseline=pars[0]+pars[1]*np.cos(4*np.pi*((series[:,0]-EPOCH)/PERIOD))
        edges=np.arange(.315,.731,.014)
        bins=np.digitize(ph,edges)-1
        for j in range(len(edges)-1):
            use=bins==j
            if sum(use)>=3:
                light.scatter(np.median(ph[use]*PERIOD*24),np.median(series[use,1]-baseline[use]),
                              color=clr,s=11,alpha=.55)
        xx=np.linspace(.32,.73,300);dd=wrap(xx+pars[3]-pars[6])
        light.plot(xx*PERIOD*24,pars[5]*np.exp(-.5*(dd/pars[7])**2),lw=2,color=clr,label=label)
    light.set(xlim=(13.5,29.7),ylim=(.39,-.07),xticks=[15,18,21,24,27],
              yticks=[0,.1,.2,.3],ylabel="Dimming in I band (mag)",xlabel="Hours after fitted primary eclipse")
    light.set_title("The secondary eclipse moved later",loc="left",fontweight="bold",pad=12)
    light.legend(loc="lower left",fontsize=8,frameon=False)
    fig.tight_layout(w_pad=2.3)
    fig.savefig(BASE/"results/eclipse-timing.png",dpi=190,bbox_inches="tight",pad_inches=.10)
    plt.close(fig)
if __name__=="__main__":
    figure()
