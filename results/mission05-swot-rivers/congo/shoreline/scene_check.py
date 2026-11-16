"""Scan common clear coverage and historic-permanent-water exclusions for each January SCL scene."""
import sys
from pathlib import Path
import numpy as np
from osgeo import gdal
import shoreline


def main(directory):
    records=[];alternates=[]
    for path in sorted(directory.glob('*202401*.tif')):
        source=gdal.Open(str(path)); refgt=source.GetGeoTransform()
        for rid in shoreline.REACHES:
            ds=gdal.Open(str(shoreline.ROOT/f'grid_{rid}.tif'));gt=ds.GetGeoTransform()
            if refgt[1]!=gt[1] or refgt[5]!=gt[5]:continue
            col=round((gt[0]-refgt[0])/20);row=round((refgt[3]-gt[3])/20)
            high=source.ReadAsArray(col,row,ds.RasterXSize,ds.RasterYSize);low=ds.GetRasterBand(2).ReadAsArray();occ=ds.GetRasterBand(4).ReadAsArray();segment=ds.GetRasterBand(5).ReadAsArray()
            clear=(segment>0)&np.isin(high,[4,5,6])&np.isin(low,[4,5,6]); flooded=clear&(high==6)&(low==5)&(occ<90)
            records.append(dict(scene=path.stem,reach=rid,common_clear_fraction=round(float(clear.sum()/(segment>0).sum()),4),flooded_bare_pixels=int(flooded.sum())))
            if flooded.sum()>25:
                result,_=shoreline.sample(rid,1,90,(0,0),True,alternate_path=path)
                alternates.append(dict(scene=path.stem,**result))
    shoreline.save(shoreline.ROOT/'scene_coverage.csv',records)
    shoreline.save(shoreline.ROOT/'scene_relief.csv',alternates)

if __name__=='__main__':main(Path(sys.argv[1]))
