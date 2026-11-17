"""Record exact detector geometry and raw mask words at the best-fit position."""
import io,json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from common import CACHE,HERE,RA,DE
import vis_pixels as v

def main():
    ledger=json.loads((HERE/'coverage.json').read_bytes());ranges=json.loads((HERE/'ranges.json').read_bytes())
    out=[]
    for i,obs in enumerate(ledger):
        exp=obs['product'].split('DET-')[1].split('-0000000')[0]
        all_=[]
        for f in CACHE.glob(exp+'_*SCI_hdr.bin'):
            name=f.name.split('_')[1].removesuffix('SCI')
            with f.open('rb') as fin:h=fits.Header.fromfile(fin)
            x,y=[float(q) for q in WCS(h).all_world2pix(RA[i],DE[i],0)]
            dist=np.hypot(max(-x,0,x-v.NX+1),max(-y,0,y-v.NY+1))
            item=dict(quadrant=name,x=x,y=y,distance_to_quadrant_pix=dist)
            if dist==0:
                rec=next(q for q in obs['quadrants'] if q['quadrant']==name)
                rows=np.memmap(CACHE/(exp+'_'+name+'FLG_rows.bin'),dtype='>i4',mode='r').reshape(-1,v.NX)
                k=int(np.floor(y))-rec['row_start']; l=int(np.floor(x))
                flags=rows[k:k+2,l:l+2];item.update(flag_words=flags.tolist(),invalid_bits=int(np.bitwise_or.reduce((flags&v.BAD_MASK).ravel())))
            all_.append(item)
        nearest=min(all_,key=lambda q:q['distance_to_quadrant_pix'])
        out.append(dict(exposure=exp,ra_deg=float(RA[i]),dec_deg=float(DE[i]),closest=nearest))
    (HERE/'detectors.json').write_text(json.dumps(out,indent=2)+'\n')
    print([(d['exposure'],d['closest']) for d in out])
if __name__=='__main__':main()
