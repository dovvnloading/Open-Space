"""Fetch independent image evidence at the only two-band XDF continuity match."""
from pathlib import Path
import sys,json
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
from archive_inputs import jades,other
center=(53.185777249373096,-27.785138788485664)
if __name__=='__main__':
    pixels,receipt=map(Path,sys.argv[1:])
    urls={}
    for b in ['f775w','f850lp']:
        fn=f'hlsp_xdf_hst_acswfc-60mas_hudf_{b}_v1_sci.fits';urls[fn]='https://archive.stsci.edu/hlsps/xdf/'+fn
    for rel in ['ls-dr10','ls-dr11']:
        urls[f'legacy_west_{rel}.fits']='https://www.legacysurvey.org/viewer/fits-cutout?'+urlencode(
            dict(ra=center[0],dec=center[1],layer=rel,pixscale=.262,size=96,bands='grz'))
    with ThreadPoolExecutor(3) as ex:r=list(ex.map(lambda pair:other(pixels,*pair,True),urls.items()))
    r.append(jades(pixels,center,512,'jades_west_f150w'))
    receipt.write_text(json.dumps(r,indent=2)+'\n')
