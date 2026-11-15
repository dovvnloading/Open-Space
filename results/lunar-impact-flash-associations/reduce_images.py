"""Recover the aligned display cutouts from LROC's published Zoomify browse pyramids.

Each archived before/after image is read in its native detector orientation. The
rounded corner metadata supply a starting guess. We then register the surrounding
terrain globally and locally, using median-smoothed displacements on a grid.
"""
import argparse
import concurrent.futures
import csv
from functools import lru_cache
import io
import math
from pathlib import Path
import re
from urllib.request import urlopen

import numpy as np
from PIL import Image
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter
from scipy.optimize import root
from scipy.signal import correlate

ROOT = Path(__file__).resolve().parent
REFERENCE = "M1166565879L"
CENTER = np.array([90.52255, 29.2322])
PIXEL = (277,333)  # y, x in the 520 x 520 reference crop, with (0, 0) upper left
CACHE = Path("nac-cache")

def fetch(url,key):
    path = CACHE/key
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url,timeout=45) as response:
                    data=response.read()
                path.write_bytes(data)
                break
            except OSError:
                if attempt==2:
                    raise
    return path.read_bytes()

@lru_cache(None)
def metadata(product):
    page="https://data.lroc.im-ldi.com/lroc/view_lroc/LRO-L-LROC-3-CDR-V1.0/"+product+"C"
    html=fetch(page,product+".html").decode()
    fields=dict(re.findall(r"<td>(.*?)</th>\s*<td>(.*?)</td>",html))
    base="https://"+re.search(r'"imgPath": "//"\+"([^"]+)"',html)[1]
    xml=fetch(base+"/ImageProperties.xml",product+".xml").decode()
    dims=tuple(int(re.search(k+r'="(\d+)"',xml,re.I)[1]) for k in ("WIDTH","HEIGHT","TILESIZE"))
    return fields,base,dims

def position(product,point):
    fields,_,(width,height,_)=metadata(product)
    def corner(which):
        return np.array([float(fields[which+" longitude"]),float(fields[which+" latitude"])])-point
    ul,ur,lr,ll=map(corner,["Upper left","Upper right","Lower right","Lower left"])
    def residual(fractions):
        x,y=fractions
        return ul*(1-x)*(1-y)+ur*x*(1-y)+lr*x*y+ll*(1-x)*y
    return root(residual,[.5,.5]).x * [width,height]

def native_crop(product,radius):
    _,base,(width,height,tile)=metadata(product)
    x,y=position(product,CENTER).astype(int)-radius
    sizes=[(width,height)]
    while max(sizes[-1])>tile:
        sizes.append(tuple(value//2 for value in sizes[-1]))
    offset=sum(math.ceil(w/tile)*math.ceil(h/tile) for w,h in sizes[1:])
    columns=math.ceil(width/tile)
    level=len(sizes)-1
    out=Image.new("L",(2*radius,2*radius))
    tiles=[(i,j) for i in range(max(0,x//tile),min(columns,(x+2*radius-1)//tile+1))
                 for j in range(max(0,y//tile),min(math.ceil(height/tile),(y+2*radius-1)//tile+1))]
    def one(index):
        i,j=index
        name=f"{level}-{i}-{j}.jpg"
        group=(offset+j*columns+i)//256
        data=fetch(f"{base}/TileGroup{group}/{name}",product+"_"+name)
        return (i*tile-x,j*tile-y),Image.open(io.BytesIO(data)).convert("L")
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        for xy,image in pool.map(one,tiles):
            out.paste(image,xy)
    return np.asarray(out,dtype=float)

def jacobian(product):
    origin=position(product,CENTER)
    return np.stack([(position(product,CENTER+step)-origin)/.001 for step in np.eye(2)*.001],axis=1)

def translation(search,template):
    t=template-template.mean()
    kernel=np.ones_like(t)
    def conv(arr,ker):
        return correlate(arr,ker,mode="valid",method="fft")
    variance=conv(search**2,kernel)-conv(search,kernel)**2/t.size
    score=conv(search,t)/np.sqrt(np.maximum(variance,1e-8)*(t*t).sum())
    return np.array(np.unravel_index(np.argmax(score),score.shape))-(np.array(score.shape)-1)//2

def register(reference,product):
    source=native_crop(product,510)
    jac=jacobian(product) @ np.linalg.inv(jacobian(REFERENCE))
    yy,xx=np.mgrid[-260:260,-260:260]
    sx,sy=np.einsum("ij,jkl->ikl",jac,np.stack([xx,yy]))
    def sample(dy,dx):
        return map_coordinates(source,[sy+510+jac[1,0]*dx+jac[1,1]*dy,sx+510+jac[0,0]*dx+jac[0,1]*dy],order=1)
    approx=sample(0,0)
    hi=reference-gaussian_filter(reference,4)
    dy,dx=translation(approx-gaussian_filter(approx,4),hi[180:-180,180:-180])
    coarse=sample(dy,dx).astype("uint8").astype(float)
    field=coarse-gaussian_filter(coarse,3)
    hi=reference-gaussian_filter(reference,3)
    grid=np.arange(40,reference.shape[0]-39,28)
    shifts=np.empty((len(grid),len(grid),2))
    for iy,y in enumerate(grid):
        for ix,x in enumerate(grid):
            shifts[iy,ix]=translation(field[y-37:y+38,x-37:x+38],hi[y-11:y+12,x-11:x+12])
    yy,xx=np.indices(reference.shape)
    dy,dx=[RegularGridInterpolator((grid,grid),median_filter(shifts[:,:,i],size=3),bounds_error=False,fill_value=None)((yy,xx)) for i in range(2)]
    return map_coordinates(coarse,[yy+dy,xx+dx],order=1)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("--save",type=Path,help="Write rebuilt NPZ to a separate file.")
    args=parser.parse_args()
    global CACHE
    CACHE=args.cache
    CACHE.mkdir(exist_ok=True,parents=True)
    ref=native_crop(REFERENCE,260)
    rebuilt={}
    with (ROOT/"frame_sources.csv").open() as f:
        for row in csv.DictReader(f):
            key=row["product"]
            arr=ref if key==REFERENCE else register(ref,key)
            y,x=PIXEL
            rebuilt[key]=np.rint(arr[y-27:y+28,x-27:x+28]).clip(0,255).astype("uint8")
    frozen=np.load(ROOT/"registered_crops.npz")
    for key,im in rebuilt.items():
        diff=np.abs(im.astype(float)-frozen[key])
        print(key,"mean/max absolute display-level difference",round(diff.mean(),3),int(diff.max()))
        if diff.mean()>1:
            raise ValueError("Registration exceeds the frozen-image replay tolerance")
    if args.save:
        np.savez_compressed(args.save,**rebuilt)

if __name__=="__main__":
    main()
