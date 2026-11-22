"""Optional one-filter follow-up. Run with already-published primary/wide evidence.

The two-filter catalogue and its outcome are inputs; this program writes a
separate census and packets for individual inspection. No automatic promotion.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits

import deflector
import imaging
from study import flag, load_components


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def number(value):
    n=deflector.numeric(value)
    return n if n is not None else None


def desi(centre, folder):
    columns=['targetid','survey','program','healpix','target_ra','target_dec',
             'mean_fiber_ra','mean_fiber_dec','spectype','z','zerr','zwarn','deltachi2',
             'zcat_primary','objtype','coadd_fiberstatus']
    query=(f"SELECT {','.join(columns)} FROM desi_dr1.zpix WHERE q3c_radial_query("
           f"mean_fiber_ra,mean_fiber_dec,{centre[0]:.10f},{centre[1]:.10f},7.0/3600.0)")
    adql=(f"SELECT {','.join(columns)} FROM desi_dr1.zpix WHERE 1=CONTAINS("
          f"POINT('ICRS',mean_fiber_ra,mean_fiber_dec),CIRCLE('ICRS',{centre[0]:.10f},{centre[1]:.10f},"
          "0.0019444444444444444))")
    routes=[(deflector.SQL,{'sql':query,'ofmt':'csv','out':'','drop':'False','async':'False'}),
            (deflector.TAP,{'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':adql,'MAXREC':1000})]
    receipt=dict(errors=[])
    for url,params in routes:
        try:
            resp=requests.get(url,params=params,headers=deflector.DATALAB_PUBLIC,timeout=(12,55))
            resp.raise_for_status()
            reader=csv.DictReader(io.StringIO(resp.content.decode('utf-8-sig')))
            if not reader.fieldnames or set(columns)-set(reader.fieldnames):
                raise ValueError('DESI result schema mismatch: '+resp.text[:200])
            rows=[]
            for rec in reader:
                primary=rec['zcat_primary'].lower() in ('true','t','1')
                robust=(primary and rec['objtype'].strip()=='TGT' and number(rec['zwarn'])==0 and
                        number(rec['coadd_fiberstatus'])==0 and (number(rec['deltachi2']) or 0)>=25)
                rows.append(dict(catalog='DESI DR1',identifier=rec['targetid'],classification=rec['spectype'].strip(),
                            ra=number(rec['mean_fiber_ra']),dec=number(rec['mean_fiber_dec']),
                            target_ra=number(rec['target_ra']),target_dec=number(rec['target_dec']),
                            z=number(rec['z']),z_err=number(rec['zerr']),quality_pass=robust,
                            fiber_radius_arcsec=.75,quality=rec))
            if len(rows)>=1000:
                raise ValueError('Unexpected or truncated DESI cone')
            (folder/'desi.csv.gz').write_bytes(gzip.compress(resp.content,mtime=0))
            receipt.update(url=url,parameters=params,sha256=hashlib.sha256(resp.content).hexdigest(),rows=len(rows),status='complete')
            return rows,receipt
        except Exception as exc:
            receipt['errors'].append(dict(url=url,message=str(exc)[:500]))
    return [],dict(receipt,status='unavailable')


def sdss_centres(centres, evidence):
    """One reread of DR19 for all leads, at measured fiber coordinates."""
    if not centres:
        return {},dict(status='no_eligible_positions')
    receipt_path=evidence/'catalog_input_checksums.csv'
    checks=pd.read_csv(receipt_path,dtype=str,keep_default_na=False)
    old=checks[checks.catalog.eq('SDSS DR19') & checks.url.str.contains('spAll-lite-v6_1_3.fits.gz',regex=False)]
    urls=list(dict.fromkeys(list(old.url)+[
        host+'/sas/dr19/spectro/boss/redux/v6_1_3/spAll-lite-v6_1_3.fits.gz'
        for host in ('https://dr19.sdss.org','https://data.sdss.org','https://dev-mirror.sdss.org')]))
    expected=set(old.sha256)-{''}
    rec=dict(status='unavailable',errors=[])
    ids=list(centres)
    sky=SkyCoord([centres[x][0] for x in ids],[centres[x][1] for x in ids],unit='deg')
    results={ident:[] for ident in ids}
    with tempfile.TemporaryDirectory() as tmp:
        archive=Path(tmp)/'sp.fits.gz'
        for url in urls:
            try:
                digest=hashlib.sha256()
                with requests.get(url,stream=True,timeout=(12,75)) as resp, archive.open('wb') as dest:
                    resp.raise_for_status()
                    for block in resp.iter_content(1024*1024):
                        digest.update(block)
                        dest.write(block)
                if expected and digest.hexdigest() not in expected:
                    raise ValueError('Archive differs from primary DR19 release checksum')
                expanded=Path(tmp)/'spectra.fits'
                with gzip.open(archive) as src, expanded.open('wb') as target:
                    shutil.copyfileobj(src,target)
                rec.update(url=url,sha256=digest.hexdigest(),bytes=archive.stat().st_size)
                break
            except Exception as exc:
                rec['errors'].append(str(exc)[:350])
        else:
            return results,rec
        try:
            with fits.open(expanded,memmap=True) as ff:
                data=ff[1].data
                need={'FIBER_RA','FIBER_DEC','RACAT','DECCAT','CATALOGID','CLASS','SPECPRIMARY',
                      'Z','Z_ERR','ZWARNING','SN_MEDIAN_ALL','MJD','FIELD'}
                if not need.issubset(data.names):
                    raise ValueError(f'DR19 columns absent: {need-set(data.names)}')
                for start in range(0,len(data),40000):
                    group=data[start:start+40000]
                    ra,dec=group['FIBER_RA'],group['FIBER_DEC']
                    valid=np.flatnonzero(np.isfinite(ra)&np.isfinite(dec)&(ra>=0)&(ra<360)&(np.abs(dec)<=90))
                    if not len(valid):
                        continue
                    ii,jj,_,_=search_around_sky(sky,SkyCoord(ra[valid]*u.deg,dec[valid]*u.deg),7*u.arcsec)
                    for i,j in zip(ii,valid[jj]):
                        r=group[j]
                        quality=dict(specprimary=int(r['SPECPRIMARY']),zwarning=int(r['ZWARNING']),
                                     sn_median_all=number(r['SN_MEDIAN_ALL']),field=int(r['FIELD']),mjd=int(r['MJD']))
                        robust=(quality['specprimary']==1 and quality['zwarning']==0 and quality['field']!=16181
                                and (quality['sn_median_all'] or 0)>=5)
                        cls=r['CLASS']
                        cls=cls.decode().strip() if isinstance(cls,bytes) else str(cls).strip()
                        results[ids[i]].append(dict(catalog='SDSS DR19',identifier=str(r['CATALOGID']),
                            classification=cls,ra=number(ra[j]),dec=number(dec[j]),
                            target_ra=number(r['RACAT']),target_dec=number(r['DECCAT']),
                            z=number(r['Z']),z_err=number(r['Z_ERR']),quality_pass=robust,
                            fiber_radius_arcsec=1.0,quality=quality))
                rec.update(status='complete',archive_rows=len(data),matched=sum(map(len,results.values())))
        except Exception as exc:
            rec['errors'].append(str(exc)[:450])
            return {i:[] for i in ids},rec
    return results,rec


def assess_spectra(records,centre,positions):
    """All associations use measured fiber coordinates plus nominal targets.

    Clearance is fiber radius + 0.65 arcsec to the competing object. Seeing
    wings can still contaminate a nominally isolated fiber.
    """
    galaxies,qsos=[],[]
    for row in records:
        fiber=[row['ra'],row['dec']]
        pointing=[row['target_ra'],row['target_dec']]
        if None in fiber or None in pointing:
            continue
        centre_d=deflector.separation(fiber,centre)
        d=np.array([deflector.separation(fiber,xy) for xy in positions])
        nearest=int(d.argmin())
        floor=row['fiber_radius_arcsec']+.65
        z,e=row['z'],row['z_err']
        reliable=(row['quality_pass'] and z is not None and z>0 and e is not None and 0<e<.01*(1+z))
        on_gal=(centre_d<=.50 and np.min(d)>floor and
                deflector.separation(pointing,centre)<=.65 and centre_d+.1<min(d))
        on_image=(d[nearest]<=.50 and centre_d>floor and np.delete(d,nearest).min()>floor and
                  deflector.separation(pointing,positions[nearest])<=.65)
        row.update(distance_centre_arcsec=centre_d,image_distances_arcsec=d.tolist(),nearest_image_index=nearest,
                   reliable_z=bool(reliable),isolated_at_centre=bool(on_gal),isolated_at_image=bool(on_image),
                   clearance_arcsec=floor)
        if reliable and on_gal and row['classification']=='GALAXY':
            galaxies.append(row)
        if reliable and on_image and row['classification']=='QSO':
            qsos.append(row)
    pairs=[]
    for gal in galaxies:
        for qso in qsos:
            disjoint=(gal['catalog'],gal['identifier'])!=(qso['catalog'],qso['identifier'])
            gap=qso['z']-gal['z']
            limit=max(.005,5*math.hypot(qso['z_err'],gal['z_err']))
            if disjoint and gap>limit:
                pairs.append(dict(galaxy=[gal['catalog'],gal['identifier']],quasar=[qso['catalog'],qso['identifier']],
                                  z_lens=gal['z'],z_source=qso['z'],minimum_z_gap=limit))
    conflicts=any(r.get('isolated_at_image') and r['quality_pass'] and r['classification']=='STAR' for r in records)
    if len(qsos)>1:
        zs=[r['z'] for r in qsos]
        conflicts |= max(zs)-min(zs)>.01*(1+min(zs))
    if len(galaxies)>1:
        zs=[r['z'] for r in galaxies]
        conflicts |= max(zs)-min(zs)>.01*(1+min(zs))
    return pairs,bool(conflicts)


def run(evidence_dir,output_dir):
    evidence,out=Path(evidence_dir),Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    fields=pd.read_csv(evidence/'system_dispositions.csv',dtype=str,keep_default_na=False)
    comps=load_components(evidence)
    imgs={r['id']:r for r in (json.loads(s) for s in (evidence/'imaging'/'imaging.jsonl').read_text().splitlines())}
    known={}
    for path in [evidence/'known_recovery.csv',evidence/'wide'/'known_recovery.csv']:
        if path.exists():
            for r in pd.read_csv(path,dtype=str,keep_default_na=False).itertuples():
                for ident in r.fields.split():
                    known.setdefault(ident,set()).add(r.name)
    summaries=[]
    leads={}
    for r in fields.itertuples():
        image=imgs[r.representative]
        light={key:b for key,b in image.get('bands',{}).items() if b.get('extended_light')
               and b.get('resolved_from_every_single_host')}
        colors=[t for t in image.get('two_band_photometry',[]) if t.get('within_055mag')]
        eligible=[key for key in light if any(key in t['bands'] for t in colors)]
        row=dict(system_id=r.system_id,representative=r.representative,cohort=r.cohort,
                 geometry='wide' if r.system_id[1:2]=='W' else 'primary',known_names='|'.join(sorted(known.get(r.system_id,[]))),
                 prior=flag(r.prior),original_full_morphology=flag(r.full_morphology),
                 bands_with_light_and_host_test=len(light),two_filter_colour=bool(colors),chosen_band='',
                 optical_eligibility=bool(eligible and not image.get('matched_bands')),qso_anchor=flag(r.qso_anchor),
                 prior_conflict=flag(r.stellar_or_redshift_conflict),targeted=False)
        if row['optical_eligibility'] and row['qso_anchor'] and not row['prior_conflict']:
            # Register the colour-compatible filter with strongest worst-case
            # host-vs-free-centre comparison before looking at spectroscopy.
            band=min(eligible,key=lambda k:(-min(light[k]['bic_single_host_minus_free_galaxy']),k))
            row.update(targeted=True,chosen_band=band)
            if r.representative not in leads:
                comp=comps.loc[comps.system_id==r.system_id].sort_values('source_id').to_dict('records')
                leads[r.representative]=dict(id=r.representative,components=comp,
                    center_radec_deg=light[band]['galaxy_radec_deg'],chosen_band=band,band_evidence=light[band],
                    compatible_colors=[t for t in colors if band in t['bands']])
        summaries.append(row)
    sdss,sdss_receipt=sdss_centres({k:s['center_radec_deg'] for k,s in leads.items()},evidence)
    save(out/'sdss_receipt.json',sdss_receipt)
    verdicts={}
    for ident,obj in leads.items():
        folder=out/ident
        folder.mkdir(exist_ok=True)
        centre=obj['center_radec_deg']
        positions=[[float(c['ra']),float(c['dec'])] for c in obj['components']]
        coords=SkyCoord([p[0] for p in positions],[p[1] for p in positions],unit='deg')
        origin=SkyCoord(*centre,unit='deg')
        try:
            model=imaging.fit_lens(imaging.sky_offsets(coords,origin),origin,np.zeros(2),
                                   [float(c['g']) for c in obj['components']])
            model['center_prior_basis']='measured one-filter extended-light centroid'
            model['mass_to_measured_light_arcsec']=deflector.separation(centre,model['center_radec_deg'])
        except Exception as exc:
            model=dict(physical_geometry=False,error=str(exc)[:500])
        deblends=deflector.assess([obj],out/'catalogue_deblends')
        fg_catalog=any(x['status']=='central_extended_object_and_distinct_catalogued_images' for x in deblends)
        spectra,desi_receipt=desi(centre,folder)
        spectra.extend(sdss.get(ident,[]))
        matches,conflicts=assess_spectra(spectra,centre,positions)
        lens_ok=bool(model['physical_geometry'] and model['mass_to_measured_light_arcsec']<=.35)
        result=dict(centroid_model_pass=lens_ok,catalogued_central_deblend=fg_catalog,
            secure_foreground_and_source_spectra=bool(matches),new_spectral_conflict=conflicts,
            warrants_individual_review=bool(lens_ok and fg_catalog and matches and not conflicts),
            desi_lookup=desi_receipt['status'],sdss_lookup=sdss_receipt['status'])
        verdicts[ident]=result
        save(folder/'secondary_evidence.json',dict(chosen_band=obj['chosen_band'],measured_centre_radec_deg=centre,
            image_source_ids=[str(c['source_id']) for c in obj['components']],band_evidence=obj['band_evidence'],
            color_evidence=obj['compatible_colors'],lens_model=model,independent_catalogues=deblends,
            spectra=spectra,foreground_source_pairs=matches,desi_receipt=desi_receipt,result=result))
    for row in summaries:
        if row['targeted']:
            row.update(verdicts[row['representative']])
    cols=['system_id','representative','cohort','geometry','known_names','prior','original_full_morphology',
          'bands_with_light_and_host_test','two_filter_colour','chosen_band','optical_eligibility','qso_anchor',
          'prior_conflict','targeted','centroid_model_pass','catalogued_central_deblend',
          'secure_foreground_and_source_spectra','new_spectral_conflict','warrants_individual_review',
          'desi_lookup','sdss_lookup']
    pd.DataFrame(summaries,columns=cols).to_csv(out/'field_adjudication.csv',index=False)
    known_aliases={s['representative'] for s in summaries if s['prior'] or s['known_names']}
    choose={}
    for s in summaries:
        ident=s['representative']
        if s['cohort']=='science' and s.get('warrants_individual_review') and ident not in known_aliases:
            centre=leads[ident]['center_radec_deg']
            choose.setdefault(ident,dict(system_id=s['system_id'],ra=centre[0],dec=centre[1]))
    pd.DataFrame(list(choose.values()),columns=['system_id','ra','dec']).to_csv(out/'bibliography_targets.csv',index=False)
    report=dict(git_commit=os.getenv('GITHUB_SHA'),run_id=os.getenv('GITHUB_RUN_ID'),strata={},
                independent_image_sets_targeted=len(leads),automatic_promotions=0)
    for arm in ('primary','wide'):
        for group in ('science','control','known_benchmark'):
            sample=[s for s in summaries if s['geometry']==arm and (bool(s['known_names'])
                    if group=='known_benchmark' else s['cohort']==group)]
            report['strata'][arm+'_'+group]=dict(gaia_selected_fields=len(sample),
                optical_eligible=sum(s['optical_eligibility'] for s in sample),
                targeted=sum(s['targeted'] for s in sample),
                targeted_unique_image_sets=len({s['representative'] for s in sample if s['targeted']}),
                review_fields=sum(bool(s.get('warrants_individual_review')) for s in sample))
    paths=['system_dispositions.csv','imaging/imaging.jsonl','gaia_components.csv','wide/gaia_components.csv',
           'known_recovery.csv','wide/known_recovery.csv','catalog_input_checksums.csv']
    report['input_sha256']={p:hashlib.sha256((evidence/p).read_bytes()).hexdigest() for p in paths if (evidence/p).is_file()}
    save(out/'secondary_summary.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('evidence')
    parser.add_argument('output')
    args=parser.parse_args()
    run(args.evidence,args.output)
