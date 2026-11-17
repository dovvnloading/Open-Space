"""Run the ancillary search and final per-trigger continuity diagnostic. XDF windows abut without overlap."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import subprocess
import sys

XDF_WINDOWS={'xdf':(2113,2113),'xdf_w':(1089,2113),'xdf_e':(3137,2113),'xdf_s':(2113,1089),'xdf_n':(2113,3137)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('detector',type=Path);p.add_argument('images',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--inject',action='store_true')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    specs=[]
    for field,(x,y) in XDF_WINDOWS.items():
        for band in ['606','814']:
            src=next(a.images.glob(f'*f{band}w_v1_sci.fits'))
            mask=next(a.images.glob(f'*f{band}w_v1_wht.fits'))
            specs.append((f'{field}_{band}',src,['--scale','.06','--mask',mask,'--mask-type','weight',
                    '--window',x,y,1024,1024,'--lengths','4,6,8','--widths','.16']))
    for b,obs in [('606','jety02010'),('814','jety02020')]:
        src=a.images/f'{obs}_drc.fits'
        specs.append((f'control_{b}',src,['--hdu',1,'--mask',src,'--mask-hdu',2,'--mask-type','weight',
            '--scale','.05','--window',2200,1320,384,384,'--lengths','4,6,8','--widths','.12,.20']))
    specs.append(('jades150',a.images/'jades_f150w_sci.fits',['--mask',a.images/'jades_f150w_err.fits',
                  '--mask-type','rms','--scale','.02999474','--lengths','4,6,8','--widths','.16']))
    def command(args):
        subprocess.run([str(v) for v in [sys.executable,*args]],check=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
    def run(spec):
        field,image,opts=spec;out=a.out/f'{field}.json'
        if not out.exists():command([a.detector/'wake_search.py','search',image,*opts,'--out',out])
        quality=a.out/f'{field}_profile.json'
        command([a.detector/'profile_screen.py',out,'--input-dir',a.images,'--out',quality])
    with ThreadPoolExecutor(a.jobs) as ex:list(ex.map(run,specs))
    if a.inject:
        field,image,opts=specs[0]
        command([a.detector/'wake_search.py','inject',image,*opts,'--zeropoint','26.51','--patch','224',
                 '--mu','24,25,26,27,28','--trials','8','--out',a.out/'xdf_inj.json'])
if __name__=='__main__':main()
