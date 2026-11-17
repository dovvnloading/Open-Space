"""Extract DESI DR1 RVSpecFit exposure results by TARGETID, validating upstream SHA-256 manifests."""
import argparse, csv, hashlib, io, json, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from astropy.io import fits
PREFIX='https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/iron/v1.0/rv_output/240521/healpix/'
RVS_COLS='EXPID VRAD VRAD_ERR RVS_WARN SUCCESS RR_Z RR_SPECTYPE SN_B SN_R SN_Z TEFF LOGG FEH'.split()
FIB_COLS='NIGHT MJD FIBERSTATUS'.split()
def epoch_product(group):
    (s,p,h),ids=group
    index=str(int(h)//100)
    filename=f'rvtab_spectra-{s}-{p}-{h}.fits'
    folder=f'{PREFIX}{s}/{p}/{index}/{h}/'
    manifest=f'dr1_vac_dr1_mws_iron_v1.0_rv_output_240521_healpix_{s}_{p}_{index}_{h}.sha256sum'
    for attempt in range(3):
        try:
            with urlopen(folder+filename,timeout=100) as f: data=f.read()
            with urlopen(folder+manifest,timeout=60) as f: checks=f.read().decode()
            actual=hashlib.sha256(data).hexdigest()
            assert any(line.split()[0]==actual and line.split()[-1].lstrip('*').endswith(filename)
                       for line in checks.splitlines()),filename
            break
        except Exception:
            if attempt==2: raise
            time.sleep(2**attempt)
    rows=[]
    with fits.open(io.BytesIO(data)) as f:
        rv=f['RVTAB'].data; fib=f['FIBERMAP'].data
        for a,b in zip(rv,fib,strict=True):
            assert a['TARGETID']==b['TARGETID']
            tid=str(a['TARGETID'])
            if tid not in ids: continue
            row={'targetid':tid,'survey':s,'program':p,'healpix':h}
            row.update({key.lower():a[key] for key in RVS_COLS})
            row.update({key.lower():b[key] for key in FIB_COLS})
            rows.append(row)
    return rows,{'url':folder+filename,'bytes':len(data),'sha256':actual,'n_selected_exposures':len(rows)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('targets',type=Path);ap.add_argument('coadds',type=Path);ap.add_argument('out',type=Path);args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    targetids={row.strip() for row in args.targets.read_text().splitlines() if row.strip()}
    products={}
    with args.coadds.open() as f:
        for r in csv.DictReader(f):
            if r['targetid'] in targetids: products.setdefault((r['survey'],r['program'],r['healpix']),set()).add(r['targetid'])
    with ThreadPoolExecutor(6) as pool: results=list(pool.map(epoch_product,sorted(products.items())))
    rows=[r for batch,_ in results for r in batch]
    with (args.out/'epochs.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['targetid','survey','program','healpix']+[k.lower() for k in RVS_COLS+FIB_COLS]);writer.writeheader();writer.writerows(rows)
    (args.out/'epoch_products.json').write_text(json.dumps([r for _,r in results],indent=2)+'\n')
    print(len(targetids),'IDs requested;',len({r['targetid'] for r in rows}),'IDs returned;',len(rows),'exposure measurements')
if __name__=='__main__':main()
