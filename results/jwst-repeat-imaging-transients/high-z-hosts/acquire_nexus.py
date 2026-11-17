"""Fetch eight public QDR target tables and extract small registered stamps while streaming the mosaics."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import shutil
import time
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

ROOT = 'https://ariel.astro.illinois.edu/nexus/qdr/'
EPOCHS = tuple(range(1, 9))
HERE = Path(__file__).resolve().parent
R = 32
SEED = 510509
DATES = ['2025-06-01','2025-07-30','2025-09-28','2025-11-28','2026-01-28','2026-03-28','2026-05-28','2026-08-07']


class HashReader:
    def __init__(self, response):
        self.response = response
        self.digest = hashlib.sha256()
        self.nbytes = 0

    def read(self, n=-1):
        data = self.response.read(n)
        self.nbytes += len(data)
        self.digest.update(data)
        return data


def make_targets(out):
    manifests, tables = [], []
    frozen = pd.read_csv(HERE/'input_tables.csv').set_index('epoch')
    for ep in EPOCHS:
        name = f'nexus_msa_deep_ep{ep:02d}_targets.fits'
        url = ROOT+'nirspec/'+name
        file = out/name
        if not file.exists():
            with urlopen(url, timeout=60) as inp, file.open('wb') as dest:
                shutil.copyfileobj(inp, dest)
        digest=hashlib.sha256(file.read_bytes()).hexdigest()
        if ep not in frozen.index or frozen.loc[ep,'url'] != url or frozen.loc[ep,'sha256'] != digest:
            raise ValueError(f'Frozen target input mismatch: {file}')
        manifests.append(dict(epoch=ep,url=url,sha256=digest))
        with fits.open(file) as hd:
            a=hd[1].data
            for i in range(len(a)):
                tables.append(dict(id=int(a['MSAID'][i]),ra=a['RA'][i],dec=a['DEC'][i],
                    z=a['z_VI'][i],conf=int(a['z_VI_conf'][i]),priority=a['Priority Class'][i],
                    f444=a['F444W'][i],epoch=ep, f200=a['F200W'][i],
                    e444=a['F444W_e'][i], e200=a['F200W_e'][i]))
    # The deposited receipt is immutable; do not rewrite it from cache contents.
    rows=pd.DataFrame(tables)
    good=rows[rows.conf==2].groupby('id')
    objs=good[['ra','dec','z','priority','f444']].median()
    # A disagreeing secure redshift prevents selection as a known high-z host.
    targets=objs[(good.z.min()>5) & ((good.z.max()-good.z.min())/(1+objs.z)<.015)].copy()
    targets['sample']='high_z'
    con=objs[(good.z.max()<5)&(objs.priority>=2)&(objs.f444>0)].copy()
    # Stratify in catalogue brightness to retain faint hosts and bright registration checks.
    con['bin']=pd.qcut(np.log10(con.f444),8,labels=False,duplicates='drop')
    con=pd.concat([g.sample(min(55,len(g)),random_state=SEED+i) for i,g in con.groupby('bin')])
    con['sample']='control'
    targets=pd.concat([targets,con[targets.columns]]).sort_index().reset_index()
    for name, frame in [('targets.csv',targets),('host_catalog_photometry.csv',rows[rows.id.isin(targets.loc[targets['sample']=='high_z','id'])])]:
        content=frame.to_csv(index=False).encode(); dest=out/name
        if dest.exists() and dest.read_bytes()!=content: raise ValueError(f'Cached target output mismatch: {dest}')
        if not dest.exists(): dest.write_bytes(content)
    return targets


def extract(arg):
    ep,ptype,targets,out,band=arg
    short=band=='f200w'
    token=band if short else 'f444w_60mas'
    name=f'nexus_central_deep_ep{ep:02d}_{token}_i2d_{ptype}.fits.gz'
    url=ROOT+'nircam/'+name
    stem=out/f'ep{ep:02d}_{ptype}'
    frozen = next(x for x in json.loads((HERE/band/'mosaic_inputs.json').read_text()) if x['epoch']==ep and x['kind']==ptype)
    identity=hashlib.sha256(targets.to_csv(index=False).encode()).hexdigest()
    config=dict(targets_sha256=identity, radius=R, seed=SEED, binning=2 if short else 1, source_sha256=frozen['sha256'], source_bytes=frozen['bytes'])
    npz=stem.with_suffix('.npz'); meta_path=stem.with_suffix('.json')
    if npz.exists() != meta_path.exists(): raise ValueError(f'Partial stamp cache: {stem}')
    if npz.exists() and meta_path.exists():
        stored=json.loads(stem.with_suffix('.json').read_text())
        frozen_fields=('epoch','kind','filter','binning','url','bytes','sha256','header')
        if any(stored.get(k)!=frozen.get(k) for k in frozen_fields) or stored.get('cache_identity')!=config or stored.get('stamp_sha256')!=hashlib.sha256(npz.read_bytes()).hexdigest():
            raise ValueError('Use a separate cache directory for each filter')
        return
    for trial in range(3):
        try:
            with urlopen(url, timeout=180) as response:
                reader=HashReader(response)
                with gzip.GzipFile(fileobj=reader) as stream:
                    header=fits.Header.fromfile(stream)
                    assert header['BITPIX']==-32 and header['NAXIS']==2
                    ny,nx=header['NAXIS2'],header['NAXIS1']
                    wcs=WCS(header)
                    px,py=wcs.all_world2pix(targets.ra,targets.dec,0)
                    fac=2 if short else 1
                    if nx%fac or ny%fac: raise ValueError('Mosaic cannot be block averaged')
                    px,py=(px-(fac-1)/2)/fac,(py-(fac-1)/2)/fac
                    cx,cy=np.rint(px).astype(int),np.rint(py).astype(int)
                    stamps=np.full((len(targets),2*R+1,2*R+1),np.nan,np.float32)
                    block=128
                    for y0 in range(0,ny,block):
                        n=min(block,ny-y0)
                        chunk=np.frombuffer(stream.read(n*nx*4),dtype='>f4').reshape(n,nx)
                        yb=y0//fac
                        for j in np.where((cy+R>=yb)&(cy-R<yb+n//fac)&(cx>=R)&(cx<nx//fac-R))[0]:
                            lo,hi=max(cy[j]-R,yb),min(cy[j]+R+1,yb+n//fac)
                            patch=chunk[(lo-yb)*fac:(hi-yb)*fac,(cx[j]-R)*fac:(cx[j]+R+1)*fac]
                            if short:
                                blocks=patch.reshape(hi-lo,2,2*R+1,2)
                                if ptype=='data':
                                    patch=blocks.mean(axis=(1,3))
                                else:
                                    blocks=np.where(blocks>0,blocks,np.nan)
                                    patch=np.sqrt((blocks**2).sum(axis=(1,3)))/4
                            stamps[j,lo-cy[j]+R:hi-cy[j]+R]=patch
                    stream.read()  # Validate gzip CRC and consume the FITS pad.
                assert reader.nbytes==int(response.headers['Content-Length'])
            if reader.nbytes!=frozen['bytes'] or reader.digest.hexdigest()!=frozen['sha256'] or url!=frozen['url']:
                raise ValueError('Mosaic differs from frozen receipt')
            np.savez_compressed(stem.with_suffix('.npz'),stamps=stamps,px=px,py=py)
            fields=['NAXIS1','NAXIS2','CDELT1','CDELT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','PIXAR_SR','MAG_ZP','BUNIT']
            meta=dict(epoch=ep,kind=ptype,filter=band,binning=fac,url=url,bytes=reader.nbytes,sha256=reader.digest.hexdigest(),
                      header={k:header.get(k) for k in fields}, cache_identity=config, stamp_sha256=hashlib.sha256(npz.read_bytes()).hexdigest())
            stem.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
            print(ep,ptype,reader.nbytes,flush=True)
            return
        except Exception:
            if trial==2: raise
            time.sleep(2**trial)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path('.'))
    p.add_argument('--band',choices=['f200w','f444w'],default='f444w')
    p.add_argument('--workers',type=int,default=3)
    p.add_argument('--kind',choices=['data','error','both'],default='both')
    args=p.parse_args();args.cache.mkdir(exist_ok=True,parents=True)
    targets=make_targets(args.cache)
    pts=['data','error'] if args.kind=='both' else [args.kind]
    tasks=[(ep,t,targets,args.cache,args.band) for t in pts for ep in reversed(EPOCHS)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(extract,tasks))
