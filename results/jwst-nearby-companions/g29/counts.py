"""Summarize the downloaded GOODS-S 24 micron and local Gaia field-control queries."""
import csv
from pathlib import Path
from math import pi
HERE=Path(__file__).resolve().parent
def main():
    fields=list(csv.DictReader((HERE/'catalogue/mips_control.csv').open()))
    area=pi*(60*.06)**2
    records=[]
    for lim in [80,100,120,150,175,200]:
        n=sum(float(row['m1_f_psf'])>=lim for row in fields)
        records.append(('SEIP GOODS-S MIPS 24um',f'>={lim}uJy',n,area))
    n=int(next(csv.DictReader((HERE/'catalogue/gaia_control.csv').open()))['n'])
    records.append(('Gaia DR3 all G in local cone','finite G',n,pi*(60*.1)**2))
    with (HERE/'out/control_counts.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['selection','threshold','N','area_arcmin2','density_arcmin2','expected_full_3_to_9_arcsec'])
        for name,cut,n,area in records:
            den=n/area; w.writerow([name,cut,n,area,den,den*pi*(9**2-3**2)/3600])
if __name__=='__main__': main()
