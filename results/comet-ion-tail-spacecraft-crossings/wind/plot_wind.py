"""Draw the archived population comparison and measured-flow sensitivity."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = Path(__file__).resolve().parent
year = pd.read_csv(p/"wind_2011_three_hour.csv",parse_dates=["start_utc"])
hour = pd.read_csv(p/"wind_hourly.csv",parse_dates=["start_utc"])
ace = pd.read_csv(p/"ace_flow_sensitivity.csv",parse_dates=["start_utc"])
pesa = pd.read_csv(p/"pesa_flow_sensitivity.csv",parse_dates=["start_utc"])
t = pd.Timestamp("2011-09-29 06:00",tz="UTC")
fig,(a,b,c) = plt.subplots(3,1,figsize=(9,8),gridspec_kw={"height_ratios":[1.05,0.8,1]})
fig.subplots_adjust(left=.13,right=.98,bottom=.09,top=.89,hspace=.59)
fig.suptitle("Wind’s 29 September 2011 oxygen-ion pulse",fontsize=16,fontweight="bold",x=.13,ha="left")
valid = year[year.records>=45]
a.scatter(valid.start_utc,100*valid.rate,s=6,color="#4c6b80",alpha=.65,linewidths=0)
a.scatter([t],[100*17/54],color="#c44b13",s=38,zorder=3)
a.annotate("17 of 54 scans",xy=(t,100*17/54),xytext=(-8,-3),
           textcoords="offset points",ha="right",va="center",fontsize=9)
a.xaxis.set_major_locator(dates.MonthLocator(interval=2))
a.xaxis.set_major_formatter(dates.DateFormatter("%b"))
a.set(ylabel="Scans with O⁺ (%)",ylim=(-1,35),title="Three-hour oxygen-ion detection fractions across 2011")
win=hour[(hour.start_utc>=t-pd.Timedelta(hours=6))&(hour.start_utc<t+pd.Timedelta(hours=12))]
for column,label,color in [("o_dc","O⁺ double coincidence","#cc4c16"),("o_tc","O⁺ triple coincidence","#175a80"),("o6_tc","O⁶⁺ triple coincidence","#6c6d70")]:
    b.step(win.start_utc,win[column],where="post",label=label,color=color,lw=1.5)
b.set(ylabel="Positive scans / hour",title="Wind/STICS native moments on 29 September")
b.legend(loc="upper right",frameon=False,fontsize=8)
shift=pd.Timedelta(minutes=30)
for frame,column,label,color,style in [(win,'miss_km','Wind/SWE','#175a80','-'),
        (ace,'miss_at_wind_km','ACE/SWEPAM at Wind position','#74627c',':'),
        (pesa,'miss_as_delivered_km','Wind/3DP as delivered','#27967b','--')]:
    c.plot(frame.start_utc+shift,frame[column]/1e6,label=label,color=color,ls=style,marker='.',ms=4)
c.set(ylabel="Miss distance (million km)",xlabel="29 September 2011 (UTC)",title="Ballistic distance to Elenin’s nominal path")
c.legend(loc="upper right",frameon=False,fontsize=8)
for ax in [b,c]:
    ax.axvspan(t,t+pd.Timedelta(hours=3),color="#f4d8be",alpha=.5,zorder=-1)
    ax.set_xlim(t-pd.Timedelta(hours=6),t+pd.Timedelta(hours=12))
    ax.xaxis.set_major_locator(dates.HourLocator(byhour=list(range(0,25,3))))
    ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
for ax in [a,b,c]:
    label=ax.get_title()
    ax.set_title("")
    ax.set_title(label,loc="left",fontsize=11)
    ax.spines[["right","top"]].set_visible(False)
    ax.grid(axis="y",alpha=.18)
fig.savefig(p/"wind_elenin.png",dpi=135,facecolor="white")
