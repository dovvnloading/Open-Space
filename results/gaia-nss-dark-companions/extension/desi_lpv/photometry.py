"""Independent passband check of the Gaia LPV classification at the NSS period and its harmonic."""
import csv
import hashlib
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
SID='2951241529971223808'
PERIOD=276.4694701770533
AIP='https://gaia.aip.de/tap/sync?'
ESA='https://gea.esac.esa.int/tap-server/tap/sync?'


def main():
    tables=[('epoch_photometry.xml',AIP,'gaiadr3.epoch_photometry'),
            ('dr3_lpv.csv',ESA,'gaiadr3.vari_long_period_variable'),
            ('fpr_lpv.csv',ESA,'gaiafpr.vari_long_period_variable'),
            ('fpr_rv.csv',ESA,'gaiafpr.vari_epoch_radial_velocity'),('vari_summary.csv',ESA,'gaiadr3.vari_summary')]
    receipt=[]
    for fn,endpoint,table in tables:
        url=endpoint+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=f'SELECT * FROM {table} WHERE source_id={SID}'))
        path=HERE/fn
        if not path.exists():path.write_bytes(urllib.request.urlopen(url,timeout=65).read())
        receipt.append(dict(file=fn,url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    doc=ET.parse(HERE/'epoch_photometry.xml')
    names=[c.attrib['name'] for c in doc.findall('.//{*}FIELD')]
    data=dict(zip(names,[c.text for c in doc.findall('.//{*}TD')]))
    def arr(key):
        return np.array([np.nan if v in ('None','null','NaN') else float(v) for v in data[key].split()])
    period=1/float(next(csv.DictReader((HERE/'dr3_lpv.csv').open()))['frequency'])
    rec={}
    for band in ('g','bp','rp'):
        pre='g_transit' if band=='g' else band
        t=arr('g_transit_time' if band=='g' else band+'_obs_time')
        mag=arr(pre+'_mag')
        ok=np.isfinite(t)&np.isfinite(mag)&(np.array(data['variability_flag_'+band+'_reject'].split())=='False')
        if band!='g':ok&=(np.array(data['photometry_flag_'+band+'_reject'].split())=='False')
        if band=='g':ok&=arr('g_transit_n_obs')>0
        t,mag=t[ok],mag[ok]
        # Give the same weight to each observing night to limit close-pair leverage.
        night=np.floor(t)
        t=np.array([t[night==d].mean() for d in np.unique(night)])
        mag=np.array([mag[night==d].mean() for d in np.unique(night)])
        out={}
        for label,p in [('orbital',PERIOD),('half_orbital',PERIOD/2)]:
            angle=2*np.pi*t/p
            X=np.column_stack((np.ones(len(t)),np.sin(angle),np.cos(angle)))
            co=np.linalg.lstsq(X,mag,rcond=None)[0]
            out[label]=dict(semi_amplitude_mag=float(np.hypot(*co[1:])),rms_mag=float(np.sqrt(np.mean((mag-X@co)**2))))
        rec[band]=dict(retained_transits=int(ok.sum()),nights=len(t),span_d=float(np.ptp(t)),model=out)
    (HERE/'photometry.json').write_text(json.dumps(dict(inputs=receipt,orbital_period=PERIOD,lpv_period=period,
        twice_photometric_over_orbital=2*period/PERIOD,bands=rec),indent=2)+'\n')


if __name__=='__main__':main()
