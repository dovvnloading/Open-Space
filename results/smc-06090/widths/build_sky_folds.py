"""Prepare alternate TESS aperture data with an independent sky/trend extraction."""
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'joint'))
import tess_template
sys.modules['tess_timing']=tess_template
from sky_reduction_source import reduce

def main():
    if len(sys.argv)!=2:
        raise SystemExit('Usage: python build_sky_folds.py DIRECTORY_WITH_RAW_TESS_FITS')
    records={}
    for sectors in tess_template.GROUPS.values():
        for sector in sectors:
            s,n=reduce(Path(sys.argv[1])/f's{sector:03}.fits')
            records[f's{sector}_sums']=s
            records[f's{sector}_counts']=n
    np.savez_compressed(HERE/'sky_folds.npz',**records)

if __name__=='__main__':
    main()
