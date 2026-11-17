"""Select the fixed tile-7119 bright cohort and propagate Gaia astrometry."""
import sys
from screen import load,quality,save_csv
from catalogs import HERE
sys.path.insert(0,str(HERE.parent/'spectra'))
from select_fields import pos,healpix_nest,COLUMNS

def main():
    cohort=[]
    for star in load():
        if quality(star)!='eligible' or float(star['w1mpro'])>=12:
            continue
        star['ref_epoch']='2016.0'
        ra,dec=pos(star)
        tile=healpix_nest(5,ra,dec)
        if tile==7119:
            cohort.append(dict(label=star['source_id'],**{k:star[k] for k in COLUMNS.split(',')},
                               ra_2025_9=ra,dec_2025_9=dec,tile=tile))
    save_csv('dr1_positions.csv',cohort)

if __name__=='__main__':
    main()
