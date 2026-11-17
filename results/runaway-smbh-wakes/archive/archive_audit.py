"""Summarize blind HST/JWST searches, with bidirectional geometry matching and catalogued endpoint distances."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree
from matplotlib import pyplot as plt
from wake_filter import compatible
from archive_screen import XDF_WINDOWS

def write(rows,path):
    with path.open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));wr.writeheader();wr.writerows(rows)

def convert(row,w_from,w_to):
    """Map both line endpoints between image frames without assuming they share orientation."""
    pts=np.array([[row[k] for k in ks] for ks in [('x0','y0'),('x1','y1'),('x','y')]])
    a,b,c=w_to.all_world2pix(w_from.all_pix2world(pts,0),0);length=np.linalg.norm(a-b)
    return dict(x0=a[0],y0=a[1],x1=b[0],y1=b[1],x=c[0],y=c[1],length_px=length,
                fwhm_px=row['fwhm_px']*length/row['length_px'],pa_image_deg=np.degrees(np.arctan2(*(b-a)[::-1]))%180)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('images',type=Path);p.add_argument('outputs',type=Path)
    p.add_argument('tables',type=Path);a=p.parse_args();a.tables.mkdir(parents=True,exist_ok=True)
    data=a.images;out=a.outputs
    cat=fits.getdata(next(data.glob('*cat.fits')),1);cat=cat[cat['STAR']==0]
    h606,h814=[fits.open(next(data.glob(f'*f{b}w_v1_sci.fits')))[0] for b in ['606','814']]
    jades=fits.open(data/'jades_f150w_sci.fits')[0]
    whst,wjw=[WCS(x.header).celestial for x in [h606,jades]]
    catxy=whst.all_world2pix(np.column_stack((cat['RA'],cat['DEC'])),0);tree=cKDTree(catxy)
    def nearest(row):
        sep,ids=tree.query([[row['x0'],row['y0']],[row['x1'],row['y1']]])
        k=sep.argmin();source=cat[ids[k]]
        return dict(endpoint=int(k),host_id=int(source['ID']),endpoint_arcsec=float(sep[k]*.06),
              host_ra=float(source['RA']),host_dec=float(source['DEC']),z_bpz=float(source['Z_BPZ']),
              odds_bpz=float(source['ODDS_BPZ']),specz=float(source['SPECZ']),grismz=float(source['GRISMZ']))
    def records(field,band):
        key=f'{field}_{band}' if field!='JADES' else 'jades150'
        d=json.loads((out/(key+'_profile.json')).read_text())
        return d['detections'],d['valid_geometric_deg2']
    measurements=[];pairs=[];summaries=[];stamps=[]
    for field in [*XDF_WINDOWS,'JADES']:
        if field=='JADES':
            curves,area=records(field,None)
            hh=[convert(r,wjw,whst) for r in curves];other,_=records('xdf','814')
        else:
            curves,area=records(field,'606');other,_=records(field,'814');hh=curves
        lookup=[]
        for index,r in enumerate(curves):
            xy=hh[index];pair_idx=[j for j,s in enumerate(other) if compatible(xy,s) and compatible(s,xy)]
            lookup.append(pair_idx)
            for j in pair_idx:
                s=other[j];near=nearest(xy);passed=bool(s['profile_gate'] and r['profile_gate'])
                pairs.append(dict(field=field,index=index,f814_id=j,both_profile=passed,
                   first_profile=r['profile_gate'],f814_profile=s['profile_gate'],**near))
                if passed:
                    stamps.append((field,index,j,xy,h606.data,h814.data))
        summaries.append(dict(field=field,geometric_deg2=area,primary_band='F150W' if field=='JADES' else 'F606W',
               raw=len(curves),f814_raw_comparison=len(other),profile=sum(r['profile_gate'] for r in curves),
               f814_profile_comparison=sum(r['profile_gate'] for r in other),pairs=sum(map(len,lookup)),
               both_profile=sum(p['both_profile'] for p in pairs if p['field']==field)))
        record_set=[('F150W' if field=='JADES' else 'F606W',curves,hh)]
        if field!='JADES':record_set.append(('F814W',other,other))
        for band,raws,coords in record_set:
            for i,(r,c) in enumerate(zip(raws,coords)):
                ra,dec=whst.all_pix2world([[c['x'],c['y']]],0)[0]
                measurements.append(dict(field=field,band=band,index=i,ra=ra,dec=dec,**r,**nearest(c)))
    write(summaries,a.tables/'field_totals.csv')
    write(measurements,a.tables/'archive_triggers.csv')
    write(pairs,a.tables/'archive_pairs.csv')
    def render(items,filename):
        if not items:return
        fig,axes=plt.subplots(int(np.ceil(len(items)/4)),4,figsize=(10,2.5*int(np.ceil(len(items)/4))),squeeze=False)
        for (field,i,j,row,im1,im2),ax in zip(items,axes.ravel()):
            x,y=int(row['x']),int(row['y']);rad=int(.65*row['length_px'])+10
            h,w=im1.shape;x0=max(0,x-rad);x1=min(w,x+rad+1);y0=max(0,y-rad);y1=min(h,y+rad+1)
            imgs=[]
            for image in [im1,im2]:
                z=image[y0:y1,x0:x1]; m=np.nanmedian(z);sigma=np.nanmedian(abs(z-m))*1.4826
                imgs.append(np.arcsinh(np.maximum(z-m,0)/max(2*sigma,1e-20)))
            rgb=np.dstack((imgs[1],.5*(imgs[0]+imgs[1]),imgs[0]))
            ax.imshow(np.clip(rgb/max(2.5,np.nanpercentile(rgb,99.9)),0,1),origin='lower',extent=(x0-.5,x1-.5,y0-.5,y1-.5))
            ax.scatter([row['x0'],row['x1']],[row['y0'],row['y1']],s=22,facecolors='none',edgecolors='yellow',lw=.8)
            ax.set_title(f'{field} {i} / {j}',fontsize=8);ax.set_axis_off()
        for ax in axes.ravel()[len(items):]:ax.set_axis_off()
        fig.tight_layout();fig.savefig(a.tables/filename,dpi=90);plt.close(fig)
    render(stamps,'bilateral_profile.jpg')
    render([('F150W',i,'',r,jades.data,jades.data) for i,r in enumerate(records('JADES',None)[0])], 'jades_raw.jpg')
if __name__=='__main__':main()
