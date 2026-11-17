"""Read the four DR1 spectra of source 2951241529971223808 and measure local line profiles."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "epochs"))
from catalog_io import RangeStream, byte_range, read_header, table_dtype

SOURCE = "2951241529971223808"
TARGET = 2305843031200339375
URL = "https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/main/backup/209/20969/spectra-main-backup-20969.fits.gz"
SIZE = 781361478


def read_spectra():
    # Fixed-width HDUs are traversed once; only wavelength grids and selected flux/ivar/mask rows survive.
    ranges = RangeStream(URL, SIZE)
    output = {}
    indices = None
    with ranges, gzip.GzipFile(fileobj=io.BufferedReader(ranges)) as f:
        while (h := read_header(f)) is not None:
            dims = [int(h[f'NAXIS{i+1}']) for i in range(int(h['NAXIS']))]
            if not dims:continue
            name = h.get('EXTNAME','')
            unit = abs(int(h['BITPIX']))//8
            length = math.prod(dims)*unit+int(h.get('PCOUNT',0))
            start = f.tell()
            if name == 'FIBERMAP':
                data = np.frombuffer(f.read(length), table_dtype(h), count=dims[1])
                indices = np.flatnonzero(data['TARGETID']==TARGET)
                output['EXPID']=data['EXPID'][indices].tolist()
            elif name.endswith('_WAVELENGTH') and len(dims)==1:
                dtype = '>f8' if unit==8 else '>f4'
                output[name]=np.frombuffer(f.read(length),dtype)
            elif name in {f'{c}_{p}' for c in 'BRZ' for p in ('FLUX','IVAR','MASK')}:
                dtype = '>i4' if name.endswith("MASK") else '>f4'
                assert len(dims)==2 and unit==4 and indices is not None
                rows=[]
                for i in indices:
                    f.seek(start+int(i)*dims[0]*4)
                    rows.append(np.frombuffer(f.read(dims[0]*4),dtype))
                output[name]=np.array(rows)
            f.seek(start + (length+2879)//2880*2880)
    return output,ranges.digest.hexdigest()


def parameters(ledger):
    record=next(p for p in ledger['products'] if p['url'].endswith('rvpix_exp-main-backup.fits'))
    url=record['url']; header,_=byte_range(url,0,65535)
    f=io.BytesIO(header);read_header(f); h=read_header(f);dt=table_dtype(h);offset=f.tell()
    rows=[]
    for v in record['matches']:
        if v['source_id']!=SOURCE:continue
        i=v['row']
        raw,_=byte_range(url,offset+i*dt.itemsize,offset+(i+1)*dt.itemsize-1)
        digest=hashlib.sha256(raw).hexdigest()
        if digest!=v['rvs_sha256']:raise ValueError("RVS row hash mismatch")
        r=np.frombuffer(raw,dt)[0]
        rows.append(dict(row=i,**{k:r[k].item() for k in 'EXPID VRAD VRAD_ERR TEFF TEFF_ERR LOGG LOGG_ERR FEH FEH_ERR VSINI CHISQ_TOT CHISQ_C_TOT SN_B SN_R SN_Z RVS_WARN'.split()},sha256=digest))
    return rows


def main(out):
    root=Path(__file__).resolve().parents[1]
    ledger=json.loads((root/'extension/epochs/desi_depth.json').read_text())
    params=parameters(ledger)
    vel={int(r['expid']):float(r['rv_kms']) for r in csv.DictReader((root/'extension/epochs/desi_visits.csv').open()) if r['source_id']==SOURCE}
    data,digest=read_spectra()
    rows=[]
    # Vacuum wavelengths in angstrom. Central integration and two continuum windows.
    definitions={
        "Hbeta":("B",(4850,4875),(4820,4845),(4880,4905)),
        "Halpha":("R",(6553,6577),(6520,6545),(6585,6610)),
        "CaII8544":("Z",(8533,8556),(8516,8528),(8560,8572)),
    }
    for i,expid in enumerate(data['EXPID']):
        if expid not in vel: continue
        for line,(arm,center,left,right) in definitions.items():
            wave=data[arm+'_WAVELENGTH']/(1+vel[expid]/299792.458)
            flux=data[arm+'_FLUX'][i];ivar=data[arm+'_IVAR'][i];mask=data[arm+'_MASK'][i]
            valid=np.isfinite(flux)&(ivar>0)&(mask==0)
            def win(bounds):return valid&(wave>=bounds[0])&(wave<=bounds[1])
            l,r,j=win(left),win(right),win(center)
            if min(l.sum(),r.sum(),j.sum())<4:continue
            cont=np.interp(wave[j],[np.median(wave[l]),np.median(wave[r])],[np.median(flux[l]),np.median(flux[r])])
            dx=float(np.median(np.diff(wave)))
            ew=dx*np.sum(1-flux[j]/cont)
            formal=dx*np.sqrt(np.sum(1/(ivar[j]*cont**2)))
            rows.append(dict(expid=expid,feature=line,ew_angstrom=float(ew),flux_only_err_angstrom=float(formal),median_normalized_flux=float(np.median(flux[j]/cont)),max_normalized_flux=float(max(flux[j]/cont)),valid_pixels=int(j.sum()),fraction_valid=float(j.sum()/((wave>=center[0])&(wave<=center[1])).sum())))
    assert set(data['EXPID'])>=set(vel) and len(vel)==4
    (out/'desi_lpv_spectra.json').write_text(json.dumps(dict(source_id=SOURCE,url=URL,bytes=SIZE,sha256=digest,rvtab= params,profiles=rows),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);args=p.parse_args();main(args.out)
