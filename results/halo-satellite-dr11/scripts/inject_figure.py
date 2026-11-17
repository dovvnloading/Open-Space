"""Render the field-conditional completeness map from its measured table."""
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd
from search import OUT, FIELDS


def run(off_control=False):
    prefix = 'off_control_' if off_control else ''
    table = pd.read_csv(OUT/(prefix+'completeness.csv'))
    fig, axes = plt.subplots(2,3,figsize=(9.6,5.25),sharex=True,sharey=True,layout='constrained')
    title = 'Recovery outside known satellite centres' if off_control else 'Recovery of injected stellar systems in DR11'
    fig.suptitle(title,fontsize=15)
    palette = plt.get_cmap('cividis'); scale = Normalize(0,1)
    for i,field in enumerate(FIELDS):
        for j,rh in enumerate(sorted(table.rhalf_pc.unique())):
            axis=axes[i,j]
            sub = table[(table.field==field)&(table.rhalf_pc==rh)]
            frac = sub.pivot(index='mv_expected',columns='dm',values='fraction').sort_index(ascending=False)
            successes = sub.pivot(index='mv_expected',columns='dm',values='successes').loc[frac.index,frac.columns]
            axis.imshow(frac,cmap=palette,norm=scale,aspect='auto')
            if i==0: axis.set_title(f'Half-light radius {rh:g} pc')
            if j==0: axis.set_ylabel(f'{field.title()} field\nExpected $M_V$')
            if i==1: axis.set_xlabel('Distance modulus (mag)')
            axis.set_xticks(range(len(frac.columns)),[f'{v:g}' for v in frac.columns]); axis.set_yticks(range(len(frac.index)),[f'{v:g}' for v in frac.index])
            for y in range(len(frac.index)):
                for x in range(len(frac.columns)):
                    axis.text(x,y,str(successes.iloc[y,x]),color='white' if frac.iloc[y,x]<.58 else '#141414',ha='center',va='center',fontsize=9)
    fig.colorbar(ScalarMappable(cmap=palette,norm=scale),ax=axes.ravel(),label='Recovery fraction',shrink=.86)
    with matplotlib.rc_context({'svg.fonttype':'none','svg.hashsalt':'dr11-recovery-map-v1'}):
        fig.savefig(OUT/(prefix+'completeness.svg'),metadata={'Date':None})
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--off-control',action='store_true')
    run(parser.parse_args().off_control)
