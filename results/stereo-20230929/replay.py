"""Reproduce the magnetic profile search and the independent rear-transition measurements."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REFERENCE = pd.Timestamp('2023-09-28 20:50Z')
SCALES = (0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
AU_KM = 149597870.7

def load(name):
    return pd.read_csv(ROOT/'data'/name, parse_dates=['time'], index_col='time')

def compare(template, field, scale):
    """Integer-grid starts, common linear time dilation, no rotations or sign flips."""
    size = round((len(template)-1)*scale)+1
    phase = np.linspace(0,1,size)
    model = np.column_stack([np.interp(phase,np.linspace(0,1,len(template)),v) for v in template.T])
    unit = model / np.linalg.norm(model,axis=1)[:,None]
    centered = model-model.mean(axis=0)
    y = field.to_numpy()
    good = np.isfinite(y).all(axis=1)
    y = np.nan_to_num(y)
    yunit = y/np.maximum(np.linalg.norm(y,axis=1)[:,None],1e-12)
    window = np.ones(size)
    numerator = sum(np.correlate(b,a,'valid') for a,b in zip(centered.T,y.T))
    variance = sum(np.convolve(b*b,window,'valid')-np.convolve(b,window,'valid')**2/size for b in y.T)
    shape = numerator/np.maximum(np.sqrt(np.maximum(variance,0))*np.linalg.norm(centered),1e-12)
    direction = sum(np.correlate(b,a,'valid') for a,b in zip(unit.T,yunit.T))/size
    score=(shape+direction)/2
    score[np.convolve(good.astype(float),window,'valid') < size]=np.nan
    return pd.DataFrame(dict(score=score, shape=shape, direction=direction, scale=scale),
                        index=field.index[:len(score)])

def radial_marker(frame, component, start, end, threshold=2):
    # Bin left edge for the first half hour persistently above +2 nT in the radial direction.
    left=frame.loc[start:end,component]
    if component=='ace_bx_gse':
        left=-left  # GSE +X points sunward; -Bx approximates Br near Earth.
    stable=left.gt(threshold).rolling(6,min_periods=6).sum().eq(6)
    return stable[stable].index[0]-pd.Timedelta(minutes=25)

def compute():
    stereo=load('stereo_search_10min.csv.gz')
    wind=load('wind_search_10min.csv.gz')
    event=load('encounter_5min.csv.gz')
    positions=pd.read_csv(ROOT/'data'/'positions.csv')
    template=wind.loc[REFERENCE:pd.Timestamp('2023-09-29 09:40Z')].to_numpy()
    if not np.isfinite(template).all() or len(template)!=78:
        raise ValueError('Incomplete Wind template')
    trials=[compare(template,stereo,s) for s in SCALES]
    near=[x.loc[REFERENCE-pd.Timedelta(hours=12):REFERENCE+pd.Timedelta(hours=12)-pd.Timedelta(minutes=10)] for x in trials]
    winner=max((v.loc[v.score.idxmax()] for v in near),key=lambda r:r.score)
    maxima=pd.concat([r.score for r in trials],axis=1).max(axis=1)
    day=((maxima.index-REFERENCE+pd.Timedelta(hours=12)).total_seconds()//86400).astype(int)
    blocks=maxima.groupby(day).agg(['max','count'])
    controls=blocks[(blocks['count']==144)&(blocks.index!=0)]
    marks={key:radial_marker(event,field,'2023-09-29 '+left,'2023-09-29 '+right) for
           key,field,left,right in [('wind','wind_br','08:00','13:00'),
                                     ('ace','ace_bx_gse','08:00','13:00'),
                                     ('stereo','stereo_br','15:00','18:00')]}
    sensitivity={}
    for threshold in (1, 2, 3, 4):
        wm=radial_marker(event,'wind_br','2023-09-29 08:00','2023-09-29 13:00',threshold)
        sm=radial_marker(event,'stereo_br','2023-09-29 14:00','2023-09-29 18:00',threshold)
        sensitivity[f'{threshold}_nT']=(sm-wm).total_seconds()/3600
    xyz=[]
    for sc in ['wind','stereo']:
        row=positions[(positions.spacecraft==sc)&positions.calendar_tdb.str.contains('12:00:00')].iloc[0]
        xyz.append(row[['x_au','y_au','z_au']].to_numpy(float))
    radii=np.linalg.norm(xyz,axis=1)
    angle=np.rad2deg(np.arccos(np.clip(xyz[0]@xyz[1]/np.prod(radii),-1,1)))
    sta_p=event.loc[winner.name:marks['stereo']].median()
    wind_p=event.loc['2023-09-28 20:55':'2023-09-29 09:45'].median()
    lead=(radii[0]-radii[1])*AU_KM/(wind_p.windplasma_speed+sta_p.stereoplasma_speed)*2/3600
    lag=(marks['stereo']-marks['wind']).total_seconds()/3600
    summary={'template_bins':len(template),'best_stereo_start_utc':str(winner.name),
             'duration_scale':winner.scale,'shape_correlation':winner['shape'],'mean_direction_cosine':winner.direction,
             'combined_score':winner.score,'full_offdate_controls':len(controls),
             'controls_equal_or_higher':int((controls['max']>=winner.score).sum()),
             'highest_control_score':controls['max'].max(),'highest_control_offset_days':int(controls['max'].idxmax()),
             'persistent_radial_markers_utc':{k:str(v) for k,v in marks.items()},
             'lag_sensitivity_hours':sensitivity, 'wind_radius_au':radii[0],'stereo_radius_au':radii[1],'separation_deg':angle,
             'wind_speed_kms':wind_p.windplasma_speed,'stereo_speed_kms':sta_p.stereoplasma_speed,
             'stereo_median_beta':sta_p.stereoplasma_beta,'stereo_median_proton_temperature_k':sta_p.stereoplasma_temperature,
             'radial_lead_expected_hours':lead,'stereo_lag_observed_hours':lag,'timing_difference_hours':lead+lag}
    return summary,blocks,event

def draw(summary,blocks,event,path):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    fig,axes=plt.subplots(3,1,figsize=(9,8),gridspec_kw={'height_ratios':[1,1,.8]},constrained_layout=True)
    colors=['#ad4d28','#18645e','#4863a2']
    for sc,ax in zip(['wind','stereo'],axes[:2]):
        for comp,color in zip(['br','bt','bn'],colors):
            ax.plot(event.index,event[f'{sc}_{comp}'],label='B'+comp[1].upper(),color=color,lw=1)
        marker=pd.Timestamp(summary['persistent_radial_markers_utc'][sc]);ax.axvline(marker,color='#333',linestyle='--',lw=1)
        ax.set_title(f"{'Wind' if sc=='wind' else 'STEREO-A'}   ·   sustained radial recovery {marker:%H:%M} UTC",loc='left')
        ax.set_xlim(pd.Timestamp('2023-09-28 16:00Z'),pd.Timestamp('2023-09-29 19:00Z'))
        ax.set_ylabel('RTN field (nT)');ax.set_ylim(-9,9);ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b\n%H:%M'))
    axes[0].legend(loc='lower left',ncol=3,fontsize=8)
    ctl=blocks[(blocks['count']==144)&(blocks.index!=0)]
    axes[2].hist(ctl['max'],bins=30,color='#56808f')
    axes[2].axvline(summary['combined_score'],color='#ad4d28',lw=2,label='28–29 September')
    axes[2].set_xlabel('Best field-profile score in each complete 24-hour search window')
    axes[2].set_ylabel('Off-date windows');axes[2].legend(loc='upper right')
    fig.savefig(path,dpi=130);plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plot',type=Path)
    opt=parser.parse_args()
    result,blocks,event=compute()
    print(json.dumps(result,indent=2))
    if opt.plot:draw(result,blocks,event,opt.plot)
