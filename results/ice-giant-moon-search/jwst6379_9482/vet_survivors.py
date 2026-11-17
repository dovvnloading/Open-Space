"""Compute native integration stack contrasts for every surviving positive linked alias."""
import argparse,csv,json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import warnings
import numpy as np
from scipy import ndimage
sys.path.insert(0,str(Path(__file__).resolve().parent))
from native_orbits import search,trials


def evaluate(c,b):
    yy,xx=np.mgrid[-12:13,-12:13]
    ring=(xx**2+yy**2>=6**2)&(xx**2+yy**2<=12**2)
    mid=(xx**2+yy**2<=3**2)
    subset=ring|mid
    offx,offy=xx[subset],yy[subset]
    isring=ring[subset];ismid=mid[subset]
    n=len(c['image']);summary=[]
    for start in range(0,len(b['x']),40):
        k=slice(start,start+40);samples=[]
        for j,(im,valid) in enumerate(zip(c['image'],c['valid'],strict=True)):
            coords=[b['y'][k,j,None]+offy,b['x'][k,j,None]+offx]
            sig=ndimage.map_coordinates(im,coords,order=1,mode='constant')
            ok=ndimage.map_coordinates(valid.astype('f4'),coords,order=1,mode='constant')>=.85
            samples.append(np.where(ok,sig,np.nan))
        arr=np.array(samples);center=np.flatnonzero((offx==0)&(offy==0)).item()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            stacks=[]
            for sub in [arr,arr[:n//2],arr[n//2:]]:
                s=np.nanmedian(sub,axis=0)
                bg=np.nanmedian(s[:,isring],axis=1)
                mad=1.4826*np.nanmedian(abs(s[:,isring]-bg[:,None]),axis=1)
                stacks.append((s-bg[:,None])/mad[:,None])
        for i,z in enumerate(stacks[0]):
            keep=ismid&np.isfinite(z)
            idx=np.argmax(np.where(keep,z,-np.inf))
            summary.append([np.isfinite(arr[:,i,center]).sum(),z[center],z[idx],offx[idx],offy[idx],
                            stacks[1][i,idx],stacks[2][i,idx]])
    return np.asarray(summary)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['map1','map2','links','veto','output']:p.add_argument(k,type=Path)
    a=p.parse_args()
    lr=list(csv.DictReader((a.links/'passing_pairs.csv').open()))
    keep=[int(r['input_row']) for r in csv.DictReader((a.veto/'link_distances.csv').open())
          if r['early_sign']==r['late_sign']=='1' and r['veto']=='False']
    elem=np.zeros((len(keep),6));elem[:,0]=[float(lr[k]['a_km']) for k in keep];elem[:,-1]=[float(lr[k]['phase_first_deg']) for k in keep]
    epoch=json.loads((a.links/'parameters.json').read_text())['epoch1_tdb_jd']
    scores=[]
    with TemporaryDirectory() as td:
        t=Path(td)
        for path in [a.map1,a.map2]:
            b=trials(t/'e.csv',path,t/'b.npz',elem,epoch)
            scores.append(evaluate(search.load_npz(path),b))
    fields=['n','center','peak','dx','dy','early_at_peak','late_at_peak']
    with a.output.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['row','a_km','phase_first_deg',*['feb_'+x for x in fields],*['oct_'+x for x in fields]])
        for i,k in enumerate(keep):w.writerow([k,elem[i,0],elem[i,-1],*scores[0][i],*scores[1][i]])


if __name__=='__main__':main()
