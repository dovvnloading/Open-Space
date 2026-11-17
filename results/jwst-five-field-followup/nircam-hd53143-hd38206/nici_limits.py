"""Evaluate both old sky tracks at the published 2008 NICI H-band completeness radii for HD 38206."""
import argparse,csv,hashlib,json,subprocess,tempfile
from datetime import datetime,timezone
from urllib.request import urlopen
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import get_body_barycentric

def main(out,cache):
    url='https://arxiv.org/pdf/1307.0818';pdf=cache/'wahhaj2013.pdf'
    if not pdf.exists():
        with urlopen(url,timeout=90) as r:pdf.write_bytes(r.read())
    # Extract only target-specific tabulated values used below.
    body=subprocess.check_output(['pdftotext','-layout',str(pdf),'-']).decode()
    selected=[s.strip() for s in body.splitlines() if s.strip().startswith('HD 38206')]
    tab=next(s for s in selected if 'H-band' in s).split('H-band')[1].split()
    delta=[float(v) for v in tab[3:9]];radius=np.array([1,1.5,2,3,4,5])
    assert any('2008-12-18' in s and 'ADI' in s for s in selected)
    assert len(delta)==6 and delta==[12.6,14.6,15.4,16.1,16.2,16.2]
    epoch=Time('2008-12-18T06:00:00',scale='utc');rows=[]
    host=next(r for r in csv.DictReader((out/'HD38206-gaia.csv').open()) if r['source_id']=='2967287596507432704')
    ra,de=np.deg2rad([float(host['ra']),float(host['dec'])]);mu=np.array([float(host['pmra']),float(host['pmdec'])]);plx=float(host['parallax'])
    ax=np.array([[-np.sin(ra),np.cos(ra),0],[-np.sin(de)*np.cos(ra),-np.sin(de)*np.sin(ra),np.cos(de)]])
    for v in ['hd38206_o141','hd38206_o150']:
        meta=json.loads((out/(v+'_measurement.json')).read_text());h=fits.getheader(cache/meta['science_filename'],1);t=Time(meta['expstart_mjd'],format='mjd')
        delta_xyz=get_body_barycentric('earth',epoch).xyz.to_value('au')-np.array([h['JWST_'+k] for k in 'XYZ'])/149597870.7
        bg=(-mu*(epoch.jyear-t.jyear)+plx*ax@delta_xyz)/1000
        for trigger in csv.DictReader((out/(v+'_peaks.csv')).open()):
            origin=np.array([float(trigger['east_arcsec']),float(trigger['north_arcsec'])])
            for name,rel in [('co-moving',origin),('stationary',origin+bg)]:
                rho=np.linalg.norm(rel)
                # Report nearest tabulated bracket; 5-7 arcsec coverage fraction varies, hence no position-specific constraint.
                rad=[*radius,7];dm=[*delta,float(tab[9])]
                hi=np.searchsorted(rad,rho);lo=hi-1
                rows.append(dict(visit=v,x=trigger['x0'],y=trigger['y0'],hypothesis=name,east_arcsec=rel[0],north_arcsec=rel[1],rho_arcsec=rho,lower_r=rad[lo],upper_r=rad[hi],delta_H_lo=dm[lo],delta_H_hi=dm[hi],apparent_Hmag_lo=5.84+dm[lo],apparent_Hmag_hi=5.84+dm[hi],coverage='annular limits; up to 7 arcsec only 86% coverage at the 7 arcsec radius' if rho>5 else 'tabulated 95% completeness at bracketing radii; exact PA completeness unavailable'))
    with (out/'NICI_HD38206_loci.csv').open('w') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'NICI_HD38206_receipt.json').write_text(json.dumps(dict(url=url,doi='10.1088/0004-637X/773/2/179',sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),bytes=pdf.stat().st_size,observed_date='2008-12-18',assumed_hour='06:00 UTC for parallax',host_Hmag=5.84,contrast_source='Table 4 95% completeness, 1 to 7 arcsec; Table 6/7 has no HD 38206 candidate row; no source position recovered',checked_utc=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);p.add_argument('cache',type=Path);a=p.parse_args();main(a.out,a.cache)
