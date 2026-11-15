"""Render the multipoint event from the bundled measurements."""
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

from replay import DATA, ROOT, calculate


def plot():
    result = calculate()
    fig, axes = plt.subplots(3,1, figsize=(9.2,7.5), sharex=True, layout='constrained')
    colors = {'ace_mag':'#34855a', 'wind_mag':'#245483', 'dscovr':'#af5636'}
    for k,label in [('ace_mag','ACE'),('wind_mag','Wind'),('dscovr','DSCOVR')]:
        axes[0].plot(DATA.index, DATA[k+'.B'].rolling(5,center=True).median(), label=label,color=colors[k],lw=1.6)
    axes[0].set_ylabel('Upstream field (nT)')
    for k,label,c in [('ace_swe','ACE / SWEPAM',colors['ace_mag']),('wind_pp','Wind / 3DP',colors['wind_mag'])]:
        series=DATA[k+'.Np'].dropna();axes[1].plot(series.index,series,'.',ms=2.8,color=c,label=label)
    axes[1].set_ylabel('Protons (cm$^{-3}$)')
    geo=pd.read_csv(ROOT/'goes.csv',index_col=0,parse_dates=True)
    for k,c in [('g16','#5d5c9a'),('g18','#9f4f78')]:
        b=geo[k+'.B'].where(geo[k+'.quality']==0)-result['geostationary'][k]['ambient_nT'];axes[2].plot(b.index,b,label='GOES-'+k[1:],lw=1.6,color=c)
    axes[2].plot(DATA.index,DATA['omni.H']-result['ground']['ambient_median_nT'],label='Ground SYM-H',color='#262b32',lw=1.4)
    for k in ['wind_mag','ace_mag']:
        t=pd.Timestamp(result['crossings'][k]['ballistic_arrival_at_15Re_utc']);axes[2].axvline(t,color=colors[k],lw=1.1,ls='--')
    axes[2].set_ylabel('Change from local\nbaseline (nT)');axes[2].set_xlabel('6 October 2023 (UTC)');axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    axes[-1].set_xlim(pd.Timestamp('2023-10-06T20:25Z'),pd.Timestamp('2023-10-06T23:05Z'));axes[1].set_ylim(0,20)
    for ax in axes:
        ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,ncol=3,fontsize=9,loc='upper left')
    fig.suptitle('A dense solar-wind magnetic cavity and its terrestrial pressure pulse',fontsize=15)
    fig.savefig(ROOT/'event.png',dpi=145)


if __name__=='__main__':plot()
