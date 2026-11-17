"""Run the original F444W RDI, detection and recovery on one independent exposure per epoch/roll.
No changes to the merged reducer. Each process gets only its visit's matching PSF references.
"""
import argparse, csv, json, math, os, sys, tempfile
from pathlib import Path

VISITS={
 'hd53143_o001':('03973','001','003',5350,'HD53143'),
 'hd53143_o008':('03973','008','003',5350,'HD53143'),
 'hd38206_o141':('06012','141','142',9600,'HD38206'),
 'hd38206_o150':('06012','150','151',9600,'HD38206')}

def run(key,parent,cache,out):
    from astropy.io import fits
    sys.path.insert(0,str(parent.resolve()/'nircam'))
    import nircam_reduce as rdi
    from nircam_fetch import digest
    program,obs,ref,temp,host=VISITS[key]
    # Nominal photosphere interpolated from 2MASS Ks. Publish recovery in absolute
    # microJy as primary; contrasting with a different stellar SED scales ratios.
    with (out/(host+'-2mass.csv')).open() as f:
        rows=list(csv.DictReader(f));assert len(rows)==1;mag=float(rows[0]['k_m'])
    lk,lf=2.159,4.421
    ratio=(lk/lf)**3*math.expm1(14387.77/(temp*lk))/math.expm1(14387.77/(temp*lf))
    host_jy=666.7*10**(-0.4*mag)*ratio
    with (out/'inputs.csv').open() as f: rows=list(csv.DictReader(f))
    sel=[x for x in rows if x['expected_sha256'] or
         (x['program']==str(int(program)) and x['obs'] in [str(int(obs)),str(int(ref))] and 'nrcalong_' in x['filename'])]
    targets=[x for x in sel if x['role']=='science']
    refs=[x for x in sel if x['role']=='reference']
    assert len(targets)==1 and len(refs)==5
    with tempfile.TemporaryDirectory() as stage:
        for x in sel:
            path=(cache/x['filename']).resolve()
            if digest(path)!=x['sha256']:raise IOError(path)
            (Path(stage)/path.name).symlink_to(path)
        target=targets[0]['filename'];rdi.FIELDS[key]=('jw'+program,target,len(refs))
        sys.argv=['nircam_reduce.py','--field',key,'--inputs',stage,'--out',str(out),'--host-jy',str(host_jy)]
        rdi.main()
    meta=out/(key+'_measurement.json')
    data=json.loads(meta.read_text())
    data.update(host_normalization_method='Planck monochromatic 2MASS Ks to F444W nominal wavelength; approximate',
                host_ks_2mass_vega=mag,host_temp_assumed_K=temp,host_ks_wavelength_um=lk,
                host_f444w_wavelength_um=lf,reducer='nircam/nircam_reduce.py',science_filename=target,
                reference_filenames=[x['filename'] for x in refs])
    meta.write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('visit',choices=list(VISITS));p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.visit,a.parent,a.cache,a.out)
