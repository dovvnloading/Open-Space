"""Lay out the 14 available slope-epoch image triplets, using original broker PNGs."""
from collections import defaultdict
import csv
from PIL import Image, ImageDraw, ImageFont
from ztf_tap import RAW, ROOT
def main():
    with (ROOT/'derived'/'ztf_public_program.csv').open() as f: rows=list(csv.DictReader(f))
    entries=defaultdict(list)
    for r in rows:
        if r['at_image_epoch']=='True': entries[r['oid']].append(r)
    out=ROOT/'derived'
    font=ImageFont.load_default(size=16)
    for n in range(2):
        ids=sorted(entries)[n*5:(n+1)*5]
        im=Image.new('RGB',(840,len(ids)*193+32),'white'); draw=ImageDraw.Draw(im)
        for ci,title in enumerate(['Science','Reference','Difference','Science','Reference','Difference']):
            draw.text((ci*140+10,5),title,fill='black',font=font)
        for i,oid in enumerate(ids):
            for j,r in enumerate(sorted(entries[oid],key=lambda x:float(x['jd']))):
                x=j*420; y=i*193+30
                draw.text((x+10,y+3),oid+'   MJD '+str(round(float(r['jd'])-2400000.5,5)),fill='black',font=font)
                for a,kind in enumerate(['science','template','difference']):
                    frame=Image.open(RAW/'packets'/(oid+'_'+r['candid']+'_'+kind+'.png')).convert('RGB').resize((140,140),Image.Resampling.LANCZOS)
                    im.paste(frame,(x+140*a,y+28))
                    d=ImageDraw.Draw(im); cx=x+140*a+70; cy=y+28+70
                    d.line((cx-13,cy,cx-7,cy),fill='#e1254d',width=1)
                    d.line((cx,cy+7,cx,cy+13),fill='#e1254d',width=1)
        im.save(out/f'ztf_stamps_{n+1}.jpg',quality=88)
if __name__=='__main__':main()
