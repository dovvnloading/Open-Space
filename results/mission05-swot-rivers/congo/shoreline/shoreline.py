"""Measure paired wet-season and dry-season shoreline relief on fixed Congo bank segments.

python shoreline.py prepare JAN_SCL JUL_SCL DEM JRC_OCCURRENCE
python shoreline.py analyze
Run with GDAL, NumPy and SciPy. Scene IDs and tile sources are recorded in sources.json.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
NONWATER = (4,5)
REACHES = ('13217000061','13217000091')


def prepare(jan, jul, dem, jrc):
    gdal.UseExceptions()
    ogr.UseExceptions()
    source = ogr.Open(str(ROOT/'reaches.geojson'))
    layer = source.GetLayer()
    ref = osr.SpatialReference(); ref.ImportFromEPSG(32733)
    src = layer.GetSpatialRef(); src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(src,ref)
    janimg,julimg = (gdal.Open(x) for x in (jan,jul))
    gt = janimg.GetGeoTransform()
    assert gt == julimg.GetGeoTransform()
    for feature in layer:
        rid=str(int(feature.GetField('reach_id')))
        line=feature.GetGeometryRef().Clone();line.Transform(transform)
        # Densify to <20 m for along-channel coordinates and fixed segment-side IDs.
        line.Segmentize(10)
        vertices=np.asarray(line.GetPoints())[:,:2]
        chain=np.r_[0,np.cumsum(np.hypot(*np.diff(vertices,axis=0).T))]
        xmin,xmax,ymin,ymax=line.Buffer(1040).GetEnvelope()
        col0=math.floor((xmin-gt[0])/20); row0=math.floor((gt[3]-ymax)/20)
        nx=math.ceil((xmax-gt[0])/20)-col0;ny=math.ceil((gt[3]-ymin)/20)-row0
        west,north=gt[0]+20*col0,gt[3]-20*row0
        x,y=np.meshgrid(west+np.arange(nx)*20+10,north-np.arange(ny)*20-10)
        d,idx=cKDTree(vertices).query(np.c_[x.ravel(),y.ravel()]);d=d.reshape(x.shape);idx=idx.reshape(x.shape)
        direction=vertices[np.minimum(idx+1,len(chain)-1)]-vertices[np.maximum(idx-1,0)]
        side=direction[:,:,0]*(y-vertices[idx,1])-direction[:,:,1]*(x-vertices[idx,0])>0
        seg=(chain[idx]//1000).astype(int)*2+side.astype(int)+1
        seg[(d>1000)|(idx==0)|(idx==len(chain)-1)]=0
        bounds=(west,north-ny*20,west+nx*20,north)
        def resample(path,method):
            ds=gdal.Warp('',path,format='MEM',dstSRS='EPSG:32733',outputBounds=bounds,width=nx,height=ny,resampleAlg=method)
            return ds.ReadAsArray()
        arrays=[janimg.ReadAsArray(col0,row0,nx,ny),julimg.ReadAsArray(col0,row0,nx,ny),resample(dem,'bilinear'),resample(jrc,'near'),seg]
        out=gdal.GetDriverByName('GTiff').Create(str(ROOT/f'grid_{rid}.tif'),nx,ny,5,gdal.GDT_Float32,options=['COMPRESS=DEFLATE','PREDICTOR=3'])
        out.SetProjection(ref.ExportToWkt());out.SetGeoTransform((west,20,0,north,0,-20))
        for i,(name,array) in enumerate(zip(('january_SCL','july_SCL','copdem_m','jrc_occurrence_percent','segment_by_side'),arrays),1):
            band=out.GetRasterBand(i);band.WriteArray(array);band.SetDescription(name)
        out.FlushCache();out=None


def sample(rid, radius, occurrence_cutoff, jitter, bare, alternate_path=None):
    ds=gdal.Open(str(ROOT/f'grid_{rid}.tif'));high,low,dem,occ,segment=ds.ReadAsArray()
    if alternate_path is not None:
        source=gdal.Open(str(alternate_path));full=source.GetGeoTransform();gt=ds.GetGeoTransform()
        col=round((gt[0]-full[0])/20);row=round((full[3]-gt[3])/20)
        high=source.ReadAsArray(col,row,ds.RasterXSize,ds.RasterYSize)
    common=np.isin(high,[4,5,6])&np.isin(low,[4,5,6])
    flooded=(high==6)&(low==5 if bare else np.isin(low,NONWATER))&common
    # Derive distances from confirmed opposite-class pixels, never clouds or corridor edges.
    stable_land=common&np.isin(high,NONWATER)&np.isin(low,NONWATER)
    stable_water=common&(low==6)&(high==6)
    inner=ndimage.distance_transform_edt(~stable_water)<=radius
    outer=ndimage.distance_transform_edt(~stable_land)<=radius
    # Distinct 20 m pixels on opposite margins of the wet-to-dry strip.
    # This intentionally drops narrow bars whose two margins cannot be separated.
    lo=inner&~outer;hi=outer&~inner
    valid=(segment>0)&flooded&(occ>=0)&(occ<occurrence_cutoff)
    ground=np.roll(dem,jitter,axis=(0,1))
    # Require the same connected inundated patch within each fixed bank-side segment.
    local=[]
    for cell in np.unique(segment[segment>0]):
        footprint=valid&(segment==cell)
        objects,n=ndimage.label(footprint,np.ones((3,3),int))
        for idx,box in enumerate(ndimage.find_objects(objects),1):
            region=objects[box]==idx;il=lo[box]&region;ih=hi[box]&region
            if il.sum()<3 or ih.sum()<3:continue
            dz=float(np.median(ground[box][ih])-np.median(ground[box][il]))
            local.append(dict(reach_id=rid,segment_id=int(cell),patch_id=idx,n_low=int(il.sum()),n_high=int(ih.sum()),
                              low_edge_m=round(float(np.median(ground[box][il])),3),high_edge_m=round(float(np.median(ground[box][ih])),3),relief_m=round(dz,3)))
    by_seg={row['segment_id'] for row in local}
    dz=np.array([np.median([r['relief_m'] for r in local if r['segment_id']==seg]) for seg in sorted(by_seg)])
    summ=dict(reach_id=rid,edge_distance_m=radius*20,jrc_occ_under=occurrence_cutoff,dem_shift_y=jitter[0],dem_shift_x=jitter[1],
              july_land='bare' if bare else 'bare_or_vegetation',clear_fraction=round(np.mean(common[segment>0]),4),
              flooded_pixels=int(valid.sum()),n_patches=len(local),n_segments=len(dz),
              median_segment_relief_m=round(float(np.median(dz)),3) if len(dz) else '',
              p25_segment_relief_m=round(float(np.quantile(dz,.25)),3) if len(dz) else '',
              p75_segment_relief_m=round(float(np.quantile(dz,.75)),3) if len(dz) else '')
    return summ,local


def save(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)


def analyze():
    summaries=[];patches=[]
    options=[(1,90,(0,0),True),(2,90,(0,0),True),(1,75,(0,0),True),(1,90,(0,0),False)]
    options.extend((1,90,shift,True) for shift in [(-1,0),(1,0),(0,-1),(0,1)])
    for rid in REACHES:
        for i,args in enumerate(options):
            summary,local=sample(rid,*args);summaries.append(summary)
            if i==0:patches.extend(local)
    save(ROOT/'shoreline_summary.csv',summaries)
    if patches:save(ROOT/'primary_patches.csv',patches)


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare(*sys.argv[2:])
    analyze()
