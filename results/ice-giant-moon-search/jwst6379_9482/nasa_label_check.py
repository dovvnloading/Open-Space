"""Register the independent NASA 6379 labelled image, holding Perdita out of the fit."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from PIL import Image
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

URL=("https://assets.science.nasa.gov/content/dam/science/missions/webb/outreach/migrated/2025/"
     "STScI-01K2286H84SBFR9S9WTEZ276FT.tif/jcr:content/renditions/"
     "Full%20Res%20Annotated%20(For%20Display).png")
SHA256="aa4800ca3a3e60f4bb5bbe8bf9dc35a805eac2efe3dc4e36efa497ed520f3685"
# Approximate moon cores manually located by reading labelled NASA image, at
# 760x760 display scale. We refine within a 15x15 box on the original raster.
REFERENCE={705:(16,312),715:(222,229),709:(375,234),725:(497,210),710:(579,366),
           714:(618,421),708:(444,560),713:(334,595),712:(304,573),711:(206,352),75052:(250,320)}
ANCHORS=[705,715,709,710,714,708,713,712,711]
SUBSETS=[ANCHORS,[705,715,714,711],[705,715,714,711,708],[705,715,714,711,708,713],[705,715,714,711,75052]]


def fit(points, interp, timelimit, bodies, subset, center):
    sel=np.isin(bodies,subset)
    def model(p,parity):
        t,a,b,x,y=p
        return ((interp(t)-center)*[1,parity])@np.array([[a,b],[-b,a]])+[x,y]
    best=None
    for parity in (1,-1):
        for t in np.arange(0,timelimit,.025):
            xy=(interp(t)[sel]-center)*[1,parity]
            u=np.ones(len(xy));v=np.zeros(len(xy))
            design=np.column_stack((np.r_[xy[:,0],xy[:,1]],np.r_[-xy[:,1],xy[:,0]],np.r_[u,v],np.r_[v,u]))
            initial=np.linalg.lstsq(design,np.r_[points[sel,0],points[sel,1]],rcond=None)[0]
            r=least_squares(lambda p:(model(p,parity)[sel]-points[sel]).ravel(),[t,*initial],
                            bounds=([-.006,-10,-10,-5000,-5000],[timelimit+.006,10,10,5000,5000]),
                            loss='soft_l1',f_scale=3.)
            if best is None or r.cost<best[0]:best=(r.cost,r.x,parity,model(r.x,parity))
    return best[1:]


def run(maps,centroids,cache,out):
    if not cache.exists():
        with urlopen(URL,timeout=60) as f:cache.write_bytes(f.read())
    raw=cache.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA256:raise ValueError('NASA image checksum changed')
    img=Image.open(io.BytesIO(raw))
    if img.size!=(1140,1140):raise ValueError('unexpected released-image dimensions')
    grey=np.asarray(img.convert('L'),float)
    points=[]
    for x,y in REFERENCE.values():
        x,y=round(x*1.5),round(y*1.5)
        yy,xx=np.mgrid[y-7:y+8,x-7:x+8];stamp=grey[yy,xx]
        weights=np.maximum(stamp-np.percentile(stamp,25),0)
        points.append([np.sum(xx*weights)/weights.sum(),np.sum(yy*weights)/weights.sum()]
                      if weights.sum()>0 else [x,y])
    points=np.array(points);bodies=np.array(list(REFERENCE));which=list(bodies).index(725)
    with np.load(maps) as f:
        center=f['center_xy'];times=f['path_mjd'][:,f['path_mjd'].shape[1]//2]
        k=f['known_xy'][:,[list(f['known_id']).index(b) for b in bodies],:]
    start=times[0];interp=CubicSpline(times-start,k,axis=0)
    with centroids.open() as f:
        native=np.array([[float(r[c]) for c in ('mjd_utc','x_measured','y_measured')] for r in csv.DictReader(f)])
    native=native[np.isfinite(native).all(axis=1)]
    cases=[]
    for anchors in SUBSETS:
        par,flip,pred=fit(points,interp,times[-1]-start,bodies,anchors,center)
        t,a,b,x,y=par
        source=np.array([np.interp(t+start,native[:,0],native[:,j]) for j in (1,2)])
        observed=(source-center)*[1,flip]@np.array([[a,b],[-b,a]])+[x,y]
        nominal=pred[which];label=points[which]
        cases.append(dict(anchors=anchors,effective_mjd_utc=t+start,scale_display_px_per_native_px=np.hypot(a,b),
                     parity=flip,similarity_parameters=par.tolist(),anchor_rms_display_px=float(np.sqrt(np.mean(np.sum(
                         (pred[np.isin(bodies,anchors)]-points[np.isin(bodies,anchors)])**2,axis=1)))),
                     measured_xy=observed.tolist(),catalog_xy=nominal.tolist(),label_xy=label.tolist(),
                     measured_to_label_px=np.linalg.norm(observed-label),catalog_to_label_px=np.linalg.norm(nominal-label)))
    report=dict(source_url=URL,source_sha256=SHA256,cases=cases,
                centroid_input_sha256=hashlib.sha256(centroids.read_bytes()).hexdigest(),
                registration_input_sha256=hashlib.sha256(maps.read_bytes()).hexdigest())
    out.mkdir(parents=True,exist_ok=True)
    (out/'release_registration.json').write_text(json.dumps(report,indent=2)+'\n')
    primary=cases[0]
    fig,ax=plt.subplots(figsize=(7,4.8),constrained_layout=True)
    ax.imshow(img)
    for key,color,style,label in [('catalog_xy','#e0ac45','x','Nominal URA184 Perdita'),('measured_xy','#b7f5d6','+','Measured February track')]:
        ax.plot(*primary[key],color=color,marker=style,ms=14,ls='none',mew=1.7,label=label)
    ax.set_xlim(688,876);ax.set_ylim(366,267);ax.set_aspect('equal')
    ax.set_xlabel('NASA display x (pixels)');ax.set_ylabel('NASA display y (pixels)')
    ax.set_title("February registration against the discovery team's Perdita label")
    ax.legend(loc='lower right',fontsize=8,facecolor='#ffffff',framealpha=.95)
    fig.text(.02,.005,"Underlying image: NASA, ESA, CSA, M. El Moutamid, M. Hedman; processing J. DePasquale (STScI).",fontsize=7)
    fig.savefig(out/'release_registration.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('maps','centroids','input_png','output'):p.add_argument(key,type=Path)
    a=p.parse_args();run(a.maps,a.centroids,a.input_png,a.output)
