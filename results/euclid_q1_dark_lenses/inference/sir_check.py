"""Extract and test NISP red-grism line features at fixed redshift hypotheses for the two M/L leads."""
import csv
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from infer import HERE, write


def run(cache):
    cache.mkdir(parents=True,exist_ok=True)
    spectra=csv.DictReader(io.StringIO(json.loads((HERE/'lead_spectra_query.json').read_text())['raw_csv']))
    scenario = {
        '2659753895653871052': [('SPE',1.4359829426),('PHZ_PP',2.457281828),('PHZ',.17),('Legacy_DR9',.32003)],
        '2655306065654633189': [('SPE',2.540083647),('PHZ',.02),('Legacy_DR9',.720801)]}
    features=[('Halpha_complex',6564.61),('Hbeta',4862.68),('OIII5008',5008.24),('Pa_beta',12821.6)]
    records=[]
    for s in spectra:
        obj,fname=s['source_id'],s['file_name']
        tile=fname.split('_')[3]
        url=f"https://irsa.ipac.caltech.edu/api/spectrumdm/convert/euclid/q1/SIR/{tile}/{fname}?dataset_id=euclid_combspec&hdu={s['hdu_index']}"
        f=cache/(obj+'.vot')
        if not f.exists():
            f.write_bytes(urllib.request.urlopen(url,timeout=100).read())
        raw=f.read_bytes()
        t=ET.fromstring(raw).find('.//{*}TABLE')
        assert t.attrib['ID']==f'Euclid_{obj}_Combined'
        fields=[e.attrib['name'] for e in t.findall('./{*}FIELD')]
        assert fields==['WAVELENGTH','SIGNAL','MASK','QUALITY','VAR','NDITH','UNCERTAINTY']
        raw_values=[[float(c.text or 'nan') for c in tr] for tr in t.findall('.//{*}TR')]
        wave,flux,mask,quality,var,ndith,err=np.array(raw_values).T
        keep=(mask.astype(int)%2==0)&(mask<64)&(ndith>=2)&(err>0)&np.isfinite(flux+err)
        for label,z in scenario[obj]:
            for name,rest in features:
                line=rest*(1+z)
                if not 12500<line<18500:
                    continue
                sel=keep&(abs(wave-line)<=260)
                rel=wave[sel]-line
                if len(rel)<10 or not np.any(abs(rel)<25):
                    continue
                for fwhm in (50.,85.):
                    sigma=fwhm/2.35482
                    template=np.exp(-.5*(rel/sigma)**2)/(sigma*np.sqrt(2*np.pi))
                    if name=='OIII5008':
                        other=4960.295*(1+z)
                        template+=np.exp(-.5*((wave[sel]-other)/sigma)**2)/(2.98*sigma*np.sqrt(2*np.pi))
                    design=np.stack([np.ones_like(rel),rel/260,template],axis=1)/err[sel,None]
                    cov=np.linalg.inv(design.T@design)
                    coeff=cov@design.T@(flux[sel]/err[sel])
                    flux_unc=np.sqrt(cov[-1,-1])
                    records.append(dict(object_id=obj,hypothesis=label,z=z,line=name,lambda_obs_A=line,FWHM_A=fwhm,
                        n_bins=len(rel),flux_in_1e16_units=coeff[-1],error_in_1e16_units=flux_unc,
                        SNR=coeff[-1]/flux_unc,sha256=hashlib.sha256(raw).hexdigest(),url=url))
    write(HERE/'sir_line_checks.csv',records)

if __name__=='__main__':
    run(Path(sys.argv[1]))
