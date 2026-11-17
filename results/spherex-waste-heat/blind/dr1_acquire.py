"""Apply the audited cube photometry to the fixed bright field cohort."""
import sys
import warnings
import numpy as np
from catalogs import HERE,rows
sys.path.insert(0,str(HERE.parent/'spectra'))
from read_cubes import fetch_tile
from photometry import extract

def main():
    cohort=rows('dr1_positions.csv')
    spectra={}
    def reduce(row,stamp):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            strict,loose,coverage,skies,sky_sigma,nsky=extract(stamp)
        spectra[row['label']]=np.vstack([strict,loose,skies]).astype('f4')
        spectra['wavelength_um']=np.asarray(stamp['WAVELENGTH'],dtype='f8')
    fetch_tile((7119,cohort,HERE,HERE),sink=reduce)
    if len(spectra)!=len(cohort)+1: raise ValueError('Incomplete spectral extraction')
    np.savez_compressed(HERE/'dr1_spectra_and_blanks.npz',**spectra)


if __name__=='__main__':
    main()
