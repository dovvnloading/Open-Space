"""Render the clean, extragalactic Q1 MER morphology sample in frozen object-ID order."""
import argparse
import csv
from pathlib import Path
from PIL import Image,ImageDraw

def load(path):
    with path.open() as f: return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selected',type=Path,required=True)
    p.add_argument('--tiles',type=Path,required=True)
    p.add_argument('--stamps',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    deep={r['tile'] for r in load(a.tiles) if r['field']!='LDN1641'}
    targets=sorted((r for r in load(a.selected) if int(r['mer_det_quality_flag'] or -1)==0 and r['tileid'] in deep),key=lambda r:int(r['object_id']))
    a.out.mkdir(exist_ok=True,parents=True)
    index=[]
    cw,ch=366,205
    for i,r in enumerate(targets):
        j=i%15
        if j==0:
            pg=Image.new('RGB',(3*cw,28+5*ch),'white')
            d=ImageDraw.Draw(pg)
            d.text((5,6),'Each pair spans the same 42 arcsec    VIS left    NISP H right    cyan marks bracket the catalogued position',fill='black')
        x,y=(j%3)*cw,28+(j//3)*ch
        oid=r['object_id']
        d.text((x+2,y),f'{i+1:03d}   {oid}',fill='black')
        path=a.stamps/(oid+'.jpg')
        if path.is_file():
            im=Image.open(path).convert('RGB')
            if im.size!=(360,180): raise ValueError('unexpected thumbnail size')
            pg.paste(im,(x,y+16))
            for cx in (x+90,x+270):
                cy=y+106
                # Four outward ticks leave the center and both diagonals clear.
                for k in (-1,1):
                    d.line([(cx+9*k,cy),(cx+13*k,cy)],fill='#43cbd3',width=1)
                    d.line([(cx,cy+9*k),(cx,cy+13*k)],fill='#43cbd3',width=1)
        else:
            d.text((x+20,y+60),'Image pending',fill='black')
        fn=f'catalog_review_{i//15+1:02d}.jpg'
        index.append(dict(index=i+1,object_id=oid,page=fn,panel=j+1,available=int(path.is_file())))
        if j==14 or i==len(targets)-1:
            pg.crop((0,0,pg.width,28+(j//3+1)*ch)).save(a.out/fn,quality=91)
    with (a.out/'sheet_index.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(index[0]),lineterminator='\n')
        w.writeheader();w.writerows(index)
    print(len(index),'targets,',sum(r['available'] for r in index),'images')

if __name__=='__main__': main()
