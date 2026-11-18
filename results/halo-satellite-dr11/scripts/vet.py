"""Inspect retained peaks and the two field-rank-one control peaks on real data."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from search import weights, FIELDS, ROOT
from vet_fetch import OUT

def separation(ra, dec, xc, yc):
    ra,dec,xc,yc=tuple(np.deg2rad(x) for x in (ra,dec,xc,yc))
    h=np.sin((dec-yc)/2)**2+np.cos(dec)*np.cos(yc)*np.sin((ra-xc)/2)**2
    return 60*np.rad2deg(2*np.arcsin(np.sqrt(np.clip(h,0,1))))

def references():
    full=pd.read_csv(OUT/'lvdb.csv.gz')
    refs={'lvdb':full,'harris':pd.read_csv(OUT/'harris.csv.gz'),'recent_and_discovery':pd.read_csv(OUT/'published_objects.csv')}
    for yr in (2023,2024):
        refs[f'hunt_reffert_{yr}']=pd.read_csv(OUT/f'hunt_reffert_{yr}.csv.gz').rename(columns={'Name':'name','RA_ICRS':'ra','DE_ICRS':'dec'})
    return refs

def masks(radii,aper):
    return radii<=aper, (radii>=2*aper)&(radii<=np.sqrt(5)*aper)

def flags(data, bit, column='maskbits'):
    return (data[column].to_numpy(dtype=np.int64)&(1<<bit)) != 0

def evaluate(field,tag,p,refs,data,gaia,rr,whl):
    xc,yc,dm,aper=p['ra'],p['dec'],p['dm'],p['radius_arcmin']
    selo = data.clean.eq(1)&data.star.eq(1)&data.g.between(17,23.5)&data.r.between(16,23.5)&data.gerr.between(0,.2,inclusive='neither')&data.rerr.between(0,.2,inclusive='neither')
    stars=data[selo].copy()
    stars['weight']=weights(*(stars[c].to_numpy() for c in ['g','r','gerr','rerr']),moduli=np.array([dm]))[0]
    stars['radius']=separation(stars.ra,stars.dec,xc,yc)
    support=stars.weight>=np.exp(-2)
    on,off=masks(stars.radius,aper)
    row = dict(field=field,tag=tag,ra=xc,dec=yc,dm=dm,radius_arcmin=aper,
               stars_on=on.sum(),stars_off=off.sum(),cmd_on=(on&support).sum(),cmd_off=(off&support).sum())
    hb=(stars.g-stars.r).between(-.35,.15)&(stars.g-dm).between(0,1.)
    row.update(blue_hb_box_on=(on&hb).sum(),blue_hb_box_off=(off&hb).sum())
    neighbors=[]
    for source,table in {**refs,'WHL_background':whl}.items():
        if table.empty:
            neighbors.append(dict(field=field,tag=tag,catalog=source,nearest='no entries in padded field',separation_arcmin=np.nan))
            continue
        distances=separation(table.ra,table.dec,xc,yc)
        i=np.argmin(distances)
        neighbors.append(dict(field=field,tag=tag,catalog=source,nearest=table.iloc[i]['name'],separation_arcmin=distances[i]))
    sat=refs['lvdb'].query('table=="dwarf_mw"').copy()
    d=separation(sat.ra,sat.dec,xc,yc)
    idx=np.argmin(d)
    known=sat.iloc[idx]
    matched=d[idx]<=max(aper,known.rhalf) and abs(dm-known.distance_modulus)<=1.
    row['known_system']=known['name'] if matched else ''
    row['known_offset_arcmin']=d[idx] if matched else np.nan
    gon,_=masks(separation(gaia.ra,gaia.dec,xc,yc),aper)
    bright=gaia[gaia.phot_g_mean_mag<13].copy()
    row['nearest_G_lt_13_arcmin']=separation(bright.ra,bright.dec,xc,yc).min() if len(bright) else np.nan
    row['bright_stars_on']=(gon & gaia.phot_g_mean_mag.lt(13)).sum()
    ast=stars[support & stars.ref_cat.eq('G3')].merge(gaia,left_on='ref_id',right_on='source_id',suffixes=('','_gaia'),validate='one_to_one')
    keep=(ast.ruwe<1.4)&(ast.parallax.abs()<3*ast.parallax_error+.1)&ast.pmra_error.between(0,2,inclusive='neither')&ast.pmdec_error.between(0,2,inclusive='neither')&ast.pmra_pmdec_corr.abs().lt(1)
    ast=ast[keep].copy()
    aon,aoff=masks(ast.radius,aper)
    row.update(qualified_pm_on=aon.sum(),qualified_pm_off=aoff.sum())
    if matched:
        x=ast.pmra_gaia-known.pmra
        y=ast.pmdec_gaia-known.pmdec
        v1,v2=ast.pmra_error**2+.05**2,ast.pmdec_error**2+.05**2
        cov=ast.pmra_pmdec_corr*ast.pmra_error*ast.pmdec_error
        good=(v2*x*x+v1*y*y-2*cov*x*y)/(v1*v2-cov**2)<5.991
        row.update(published_pmra=known.pmra,published_pmdec=known.pmdec,pm_compatible_on=(aon&good).sum(),pm_compatible_off=(aoff&good).sum())
    rdata=rr.copy()
    rdata['radius_arcmin']=separation(rdata.ra,rdata.dec,xc,yc)
    ron,roff=masks(rdata.radius_arcmin,aper)
    # Apply published per-star Gaia extinction where present and otherwise assume zero.
    rdata['abs_G_at_peak']=rdata.int_average_g-rdata.g_absorption.fillna(0)-dm
    compat=rdata.abs_G_at_peak.between(0,1.2)
    row.update(rr_on=ron.sum(),rr_off=roff.sum(),rr_distance_on=(ron&compat).sum(),rr_distance_off=(roff&compat).sum())
    rdata=rdata[ron|roff].copy()
    rdata['zone']=np.where(rdata.radius_arcmin<=aper,'on','off')
    rdata['tag']=tag
    rdata['field']=field
    aon_raw,aoff_raw=masks(separation(data.ra,data.dec,xc,yc),aper)
    for name,m in [('on',aon_raw),('off',aoff_raw)]:
        sel=data[m]
        row[f'extended_clean_{name}']=(sel.star.eq(0)&sel.clean.eq(1)).sum()
        row[f'gdepth_median_{name}']=sel.gdepth.median()
        row[f'rdepth_median_{name}']=sel.rdepth.median()
        for bit,label in [(1,'bright'),(11,'medium'),(12,'large_galaxy'),(13,'globular'),(17,'known_dwarf')]:
            row[f'{label}_flag_fraction_{name}']=flags(sel,bit).mean()
        row[f'forced_psf_fraction_{name}']=flags(sel,0,'fitbits').mean()
    a,b,c,e=FIELDS[field]
    row['rectangle_edge_arcmin']=60*min((xc-a)*np.cos(np.deg2rad(yc)),(b-xc)*np.cos(np.deg2rad(yc)),yc-c,e-yc)
    row['WHL_centres_on']=(separation(whl.ra,whl.dec,xc,yc)<=aper).sum()
    catalog_pm=ast[aon|aoff].copy()
    catalog_pm['zone']=np.where(catalog_pm.radius<=aper,'on','off')
    catalog_pm['tag'],catalog_pm['field']=tag,field
    pmcols=['field','tag','zone','source_id','ra','dec','g','r','weight','radius','pmra_gaia','pmdec_gaia','pmra_error','pmdec_error','pmra_pmdec_corr','ruwe','parallax','parallax_error']
    track=pd.read_csv(ROOT/'model/parsec.csv').query('stage==0')
    fig,axes=plt.subplots(ncols=3,figsize=(10.5,4),layout='constrained')
    for ax,sel,title in zip(axes[:2],[on,off],['Peak aperture','Equal-area annulus']):
        ax.scatter((stars.g-stars.r)[sel],stars.g[sel],s=4,c='#363d46',alpha=.50,rasterized=True)
        t=track.sort_values('g')
        ax.plot(t.g-t.r,t.g+dm,color='#c94d29',lw=1.4)
        ax.set(xlabel='Dereddened g − r',ylabel='Dereddened g',title=title,xlim=(-.6,1.5),ylim=(23.6,16.7))
    ax=axes[2]
    for sel,color,lab in [(aoff,'#9fa6a9','Annulus'),(aon,'#1d6592','Aperture')]:
        ax.scatter(ast.loc[sel,'pmra_gaia'],ast.loc[sel,'pmdec_gaia'],s=12,c=color,alpha=.55,label=lab)
    if matched: ax.scatter(known.pmra,known.pmdec,marker='x',s=80,c='#cc452e',label='Published PM')
    ax.set(xlim=(-6,6),ylim=(-6,6),xlabel='μα* (mas/yr)',ylabel='μδ (mas/yr)',title='CMD-selected Gaia DR3')
    ax.legend(fontsize=8)
    title=known['name'].replace('Bootes','Boötes') if matched else field.capitalize()
    fig.suptitle(f'{title} · distance modulus {dm:.1f} · {aper:.0f} arcmin aperture')
    fig.savefig(OUT/f'{tag}_{field}_diagnostic.png',dpi=112)
    plt.close(fig)
    return row,neighbors,rdata,catalog_pm[pmcols]

def run(extra):
    refs=references()
    main=pd.read_csv(ROOT/'measurement/peaks.csv')
    main=main[(main.retained)|(main['rank']==1)].copy()
    main['tag']='primary_'+main['rank'].astype(str)
    if extra:
        table=pd.read_csv(extra)
        table=table[table.retained.astype(str).str.lower().isin(['true','1'])].copy()
        table['tag']='offcontrol_'+table['rank'].astype(str)
        main=pd.concat([main,table],ignore_index=True)
    result,near,rrs,pms=[],[],[],[]
    for field,peaks in main.groupby('field'):
        data=pd.read_csv(ROOT/f'catalog/{field}.csv.gz',dtype={'ref_id':'Int64'})
        gaia=pd.read_csv(OUT/f'gaia_{field}.csv.gz',dtype={'source_id':'Int64'})
        rr=pd.read_csv(OUT/f'rrlyrae_{field}.csv.gz',dtype={'source_id':'Int64'})
        whl=pd.concat([pd.read_csv(OUT/f'whl_{field}_table{t}.csv.gz').rename(columns={col:'name','RAJ2000':'ra','DEJ2000':'dec'}) for t,col in [(3,'WHL'),(4,'WH')]],ignore_index=True)
        for _,p in peaks.iterrows():
            row,neighbors,r,pm=evaluate(field,p.tag,p,refs,data,gaia,rr,whl)
            result.append(row);near.extend(neighbors);rrs.append(r);pms.append(pm)
    for table,name in [(pd.DataFrame(result),'adjudication'),(pd.DataFrame(near),'catalogue_neighbors'),(pd.concat(rrs,ignore_index=True),'rrlyrae_apertures'),(pd.concat(pms,ignore_index=True),'cmd_selected_gaia')]:
        table.to_csv(OUT/f'{name}.csv',index=False,float_format='%.7g')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--extra-peaks',type=Path)
    run(parser.parse_args().extra_peaks)
