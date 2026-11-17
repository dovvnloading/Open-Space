"""Render paired extended context for the MER line projected near a bright resolved source."""
import argparse,csv,hashlib,io,json
from pathlib import Path
from urllib.request import urlopen
from astropy.io import fits
import numpy as np
from PIL import Image,ImageDraw

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stamps',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--object-id',default='2699220004635136134')
    a=p.parse_args()
    row=json.loads((a.stamps/(a.object_id+'.json')).read_text())
    all_im=[]
    receipts={}
    for band in ('VIS','H'):
        url=row[band+'_input']['url'].replace('size=42arcsec','size=120arcsec')
        with urlopen(url,timeout=90) as s: raw=s.read()
        receipts[band]={'url':url,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
        with fits.open(io.BytesIO(raw)) as hd:
            data=hd[0].data
            v=data[np.isfinite(data)&(data!=0)]
            med,q=np.percentile(v,[50,84])
            z=np.arcsinh((data-med)/(2*max(q-med,1e-12)))+1
            arr=(255*np.clip(np.nan_to_num(z)/4,0,1))[::-1].astype('u1')
            im=Image.fromarray(arr).resize((500,500),Image.Resampling.LANCZOS).convert('RGB')
            d=ImageDraw.Draw(im)
            d.rectangle((163,163,337,337),outline='#408c98',width=1)
            all_im.append(im)
    out=Image.new('RGB',(1000,524),'white')
    d=ImageDraw.Draw(out);d.text((8,5),f'{a.object_id}   120 arcsec per panel  |  VIS left, H right  |  cyan box: original 42 arcsec',fill='black')
    out.paste(all_im[0],(0,24));out.paste(all_im[1],(500,24))
    out.save(a.out/(a.object_id+'_context.jpg'),quality=90)
    (a.out/(a.object_id+'_context.json')).write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__':main()
