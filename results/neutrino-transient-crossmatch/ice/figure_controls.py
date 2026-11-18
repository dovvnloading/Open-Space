"""Check the four known source stars on the authors' original raster figures.

Usage: python figure_controls.py /path/to/arxiv-2507.06176v1-source.tar
Requires numpy, scipy, pillow. The pixel distance is useful only within its
own panel. The digitized outlines are historical checks, never discovery masks.
"""
import csv
import hashlib
import io
from pathlib import Path
import sys
import tarfile

import numpy as np
from PIL import Image
from scipy import ndimage as nd

PANELS = {
    # Rectangles encircle the sky-map curves; legends and axis labels are outside.
    "AT2019aalc.png": ("AT2019aalc", "IC191119A", (.17,.26,.90,.68)),
    "AT2019dsg.png": ("AT2019dsg", "IC191001A", (.15,.26,.88,.65)),
    "AT2019fdr_icecat2.png": ("AT2019fdr", "IC200530A", (.37,.22,.81,.69)),
    "txs0506+056.png": ("TXS 0506+056", "IC170922A", (.20,.21,.92,.70)),
}


def analyze(image, box):
    a = np.asarray(image.convert("RGB")).astype(np.int16)
    h,w = a.shape[:2]
    x1,y1,x2,y2 = [round(scale*pos) for scale,pos in zip([w,h,w,h],box)]
    a = a[y1:y2,x1:x2]
    red,green,blue = a.transpose((2,0,1))
    star = (red>180)&(green>180)&(blue<90)
    labels,n = nd.label(star)
    label = np.bincount(labels.ravel())[1:].argmax()+1
    ys,xs = np.where(labels==label)
    # The yellow star has five tips with unequal vertical extent. Use centroid
    # of the colored pixels; retain full symbol radius as the uncertainty band.
    x,y = round(xs.mean()),round(ys.mean())
    star_radius = np.hypot(xs-x,ys-y).max()
    masks = {"IceCat-1":(red>180)&(red-green>100)&(red-blue>80),
             "IceCat-2":(blue>100)&(blue-red>70)&(green-red>35)}
    result = {}
    for version, mask in masks.items():
        distances=[]; areas=[]
        dist_line=nd.distance_transform_edt(~mask)
        for radius in [12,16,20,24,28]:
            # Dilate and close dashed contours, fill their outer envelope and
            # shrink back by the same radius. Nested solid 50% contours cannot
            # change the filled 90% envelope once the dashes have connected.
            filled=nd.binary_fill_holes(dist_line<=radius)
            interior=nd.distance_transform_edt(filled)>radius
            boundary = interior ^ nd.binary_erosion(interior)
            ys2,xs2=np.where(boundary)
            signed=np.min(np.hypot(ys2-y,xs2-x))*(1 if interior[y,x] else -1)
            distances.append(round(float(signed),2));areas.append(int(interior.sum()))
        stable=(max(areas)-min(areas))/max(areas)<.08
        if not stable:
            status="ambiguous_dashed_envelope"
        elif min(distances)>star_radius:
            status="inside"
        elif max(distances)<-star_radius:
            status="outside"
        else:
            status="ambiguous_marker_overlaps_boundary"
        result[version]=(status,distances,areas,round(float(star_radius),2),x+x1,y+y1)
    return result


def main():
    rows=[]
    with tarfile.open(sys.argv[1]) as archive:
        for name,(source,alert,box) in PANELS.items():
            data=archive.extractfile(name).read()
            for version,(status,distances,areas,radius,x,y) in analyze(Image.open(io.BytesIO(data)),box).items():
                rows.append(dict(source=source,alert_id=alert,localization=version,result=status,
                                 marker_x_px=x,marker_y_px=y,marker_radius_px=radius,
                                 signed_edge_distance_px="|".join(map(str,distances)),filled_area_px="|".join(map(str,areas)),
                                 original_figure=name,sha256_png=hashlib.sha256(data).hexdigest()))
    with Path(__file__).with_name("figure_controls.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)


if __name__=="__main__":
    main()
