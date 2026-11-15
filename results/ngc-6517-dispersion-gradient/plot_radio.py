"""Render the two radio panels in dispersion-field.png from the frozen tables."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze import ROOT, offsets, read


def main():
    data = json.loads((ROOT/'radio-analysis.json').read_text())
    fast, gbt = read('fast.csv'), read('green-bank.csv')
    x = offsets(fast)[:,0]
    y = np.array([float(row['dm_yin']) for row in fast])
    g = np.array([float(row['dm_lynch']) for row in gbt])
    g_x = offsets(gbt)[:,0]
    colors = dict(ink='#172c40', fitted='#26659d', gbt='#a34a0f', guide='#718399',
                  grid='#dce3e8')
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
                         'text.color':colors['ink'], 'axes.labelcolor':colors['ink'],
                         'xtick.color':colors['ink'], 'ytick.color':colors['ink'],
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.edgecolor':'#a2afbc','axes.linewidth':.8,
                         'legend.frameon':False,'svg.fonttype':'none'})
    fig, (ax,bx) = plt.subplots(1,2,figsize=(10.6,4.65),gridspec_kw={'width_ratios':[1.24,1]})
    fig.subplots_adjust(left=.082,right=.978,bottom=.18,top=.80,wspace=.29)
    fig.text(.082,.947,'Pulsar dispersion across NGC 6517',fontsize=18,weight='bold')
    for a in (ax,bx):
        a.grid(axis='y',color=colors['grid'],linewidth=.75)
        a.set_axisbelow(True)
        a.tick_params(length=3)
    ax.set_title('Dispersion versus eastward offset',loc='left',fontsize=12,weight='bold',pad=9)
    line_x=np.array([x.min(),x.max()])
    c=data['east_only']['coeff']
    ax.plot(line_x,c[0]+c[1]*line_x,color=colors['fitted'],linewidth=1.15,zorder=1,
            label='FAST east-only fit')
    ax.scatter(x,y,s=31,c=colors['fitted'],edgecolors='white',linewidths=.55,
               label='FAST 2024 · 14 pulsars',zorder=3)
    ax.scatter(g_x,g,s=44,marker='D',facecolors='none',edgecolors=colors['gbt'],
               linewidths=1.25,label='GBT 2011 · 4 of the 14',zorder=4)
    for name, delta in {'D':(-9,-12),'L':(4,6),'G':(4,6),'I':(5,2),'H':(5,2)}.items():
        i=data['data_order'].index(name)
        ax.annotate(name,(x[i],y[i]),xytext=delta,textcoords='offset points',fontsize=9)
    ax.set(xlim=(-21,82),ylim=(173.5,187),xlabel='East of the optical center (arcsec)',
           ylabel='Dispersion measure (pc cm⁻³)')
    ax.legend(loc='lower left',fontsize=8,handlelength=1.7,labelspacing=.5)
    bx.set_title('GBT prediction on 10 other pulsars',loc='left',fontsize=12,weight='bold',pad=9)
    h=data['historical_gbt']
    pred,obs = np.array(h['predicted']), np.array(h['measured'])
    bx.plot([177,186.3],[177,186.3],linestyle='--',linewidth=1.15,color=colors['guide'])
    bx.scatter(pred,obs,s=40,c=colors['fitted'],edgecolors='white',linewidths=.65,zorder=3)
    names='EFGHIKLMNO'
    for name,delta in {'L':(-11,4),'G':(-14,1),'F':(-12,6),'H':(5,0),'I':(6,-4)}.items():
        i=names.index(name)
        bx.annotate(name,(pred[i],obs[i]),xytext=delta,textcoords='offset points',fontsize=9)
    bx.text(.05,.96,'Four GBT sources set the model\nPrediction RMSE  1.05 pc cm⁻³',
            va='top',transform=bx.transAxes,fontsize=8.5,linespacing=1.5)
    bx.set(xlim=(176.9,186.3),ylim=(176.9,186.3),
           xticks=[178,180,182,184,186],yticks=[178,180,182,184,186],
           xlabel='Predicted from GBT (pc cm⁻³)',ylabel='Measured by FAST (pc cm⁻³)')
    fig.text(.082,.045,'GBT = Green Bank Telescope. The dashed line marks exact agreement.',
             color='#44596a',fontsize=8.5)
    fig.savefig(ROOT/'dispersion-field.png',dpi=155,facecolor='white',
                metadata={'Description':'Archival measurements of NGC 6517 pulsars. DM falls to the east. A model based on four old GBT pulsars predicts DM in ten distinct FAST pulsars.'})
    plt.close(fig)


if __name__=='__main__':
    main()
