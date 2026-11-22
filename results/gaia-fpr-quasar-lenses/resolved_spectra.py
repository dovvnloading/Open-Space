"""Targeted examination of DESI DR1 and SDSS DR19 coadded spectral pixels."""
import gzip
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sys
import tempfile
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
from bs4 import BeautifulSoup
from scipy.optimize import least_squares, minimize_scalar
from scipy.stats import chi2

from deflector import desi_spectra, separation
from secondary import sdss_centres
from study import flag, load_components

# Rest-frame vacuum wavelengths in angstroms. Profiles are measured in observed
# pixels, without instrumental deconvolution. Close blends get separate labels.
LINES={'LyAlpha':1215.67,'CIV_complex':1549.48,'CIII_complex':1908.73,'MgII_complex':2798.75,
       'OII_complex':3728.48,'Hgamma':4341.68,'Hbeta':4862.68,'OIII4960':4960.30,
       'OIII5008':5008.24,'NII6550':6549.86,'Halpha':6564.61,'NII6585':6585.27,
       'SII6718':6718.29,'SII6733':6732.67}
PROFILE={'CIV_complex','CIII_complex','MgII_complex','Hbeta','Hgamma'}
C=299792.458


def save(path,payload):
    path.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')


def targets(evidence,followup):
    fields=pd.read_csv(evidence/'system_dispositions.csv',dtype=str,keep_default_na=False)
    sec=pd.read_csv(followup/'field_adjudication.csv',dtype=str,keep_default_na=False)
    optics=set(sec.loc[sec.optical_eligibility.map(flag),'representative'])
    prior=set(fields.loc[fields.prior.map(flag),'representative'])
    prior.update(sec.loc[sec.known_names.ne('') | sec.prior.map(flag),'representative'])
    eligible=fields.loc[fields.cohort.eq('science') &
                        (fields.full_morphology.map(flag) | fields.representative.isin(optics))]
    imgs={r['id']:r for r in map(json.loads,(evidence/'imaging'/'imaging.jsonl').read_text().splitlines())}
    comp=load_components(evidence)
    result={}
    for ident,group in eligible.groupby('representative'):
        if ident in prior:
            continue
        sysid=group.system_id.iloc[0]
        rows=comp[comp.system_id.eq(sysid)].sort_values('source_id')
        img=imgs[ident]
        centre=img.get('lens_model',{}).get('center_radec_deg')
        basis='two-filter SIE+shear centre'
        secpath=followup/ident/'secondary_evidence.json'
        if not any(group.full_morphology.map(flag)):
            if secpath.is_file():
                centre=json.loads(secpath.read_text())['measured_centre_radec_deg']
                basis='single-filter extended-light centroid'
            else:
                colors={band for m in img.get('two_band_photometry',[]) if m.get('within_055mag') for band in m['bands']}
                bands=[(k,v) for k,v in img.get('bands',{}).items() if k in colors and
                       v.get('extended_light') and v.get('resolved_from_every_single_host') and v.get('galaxy_radec_deg')]
                if bands:
                    _,band=min(bands,key=lambda item:(-min(item[1]['bic_single_host_minus_free_galaxy']),item[0]))
                    centre=band['galaxy_radec_deg']
                    basis='single-filter extended-light centroid'
        if centre is None or not np.all(np.isfinite(centre)):
            raise ValueError(f'No valid centre for {ident}')
        result[ident]=dict(aliases=sorted(group.system_id),center=centre,center_basis=basis,
                           images=rows[['ra','dec']].to_numpy(float).tolist(),ids=rows.source_id.astype(str).tolist())
    return result


def spath(cache,url,getlog):
    if url in cache:
        if cache[url] is None:
            raise ValueError('previous retrieval failed')
        return cache[url]
    path=cache['dir']/('s-'+hashlib.sha256(url.encode()).hexdigest()[:24]+'.fits')
    digest=hashlib.sha256()
    try:
        with requests.get(url,stream=True,timeout=(10,70)) as response, path.open('wb') as dest:
            response.raise_for_status()
            for chunk in response.iter_content(1024*1024):
                digest.update(chunk)
                dest.write(chunk)
        with fits.open(path,memmap=True) as h:
            if len(h)<2:
                raise ValueError('FITS has no data extension')
        getlog.append(dict(url=url,bytes=path.stat().st_size,sha256=digest.hexdigest()))
    except Exception as exc:
        cache[url]=None
        getlog.append(dict(url=url,error=str(exc)[:250]))
        if path.exists():
            path.unlink()
        raise
    cache[url]=path
    return path


def as_float(value):
    try:
        x=float(value)
        return x if np.isfinite(x) else None
    except (TypeError,ValueError):
        return None


def legacy_spectra(systems,cache,log):
    """DR19's archival SDSS I-IV redshift catalogue plus plug-coordinate provenance."""
    result={k:[] for k in systems}
    if not systems:
        return result,dict(status='no_inputs')
    rec=dict(status='unavailable',errors=[])
    urls=[f'https://{host}/sas/dr17/sdss/spectro/redux/specObj-dr17.fits'
          for host in ('dr17.sdss.org','data.sdss.org','dev-mirror.sdss.org')]
    keys=[k for k,s in systems.items() for _ in [s['center'],*s['images']]]
    coords=[x for s in systems.values() for x in [s['center'],*s['images']]]
    anchors=SkyCoord(coords,unit='deg')
    for url in urls:
        try:
            file=spath(cache,url,log)
            with fits.open(file,memmap=True) as h:
                tab=h[1].data
                col={s.upper():s for s in tab.names}
                def resolve(*names):
                    return next((col[n] for n in names if n in col),None)
                ra_c,dec_c=resolve('PLUG_RA','RA'),resolve('PLUG_DEC','DEC')
                fiber_c,plate_c=resolve('FIBERID','FIBER'),resolve('PLATE','PLATEID')
                if not all([ra_c,dec_c,fiber_c,plate_c]) or not {'MJD','Z','Z_ERR','CLASS','ZWARNING','RUN2D'}<=set(col):
                    raise ValueError('Legacy spectrum catalogue schema missing: '+str(list(col)))
                if len(tab)<4_000_000:
                    raise ValueError('Unexpected DR17 specObj record count')
                matches={k:{} for k in systems}
                for off in range(0,len(tab),45000):
                    group=tab[off:off+45000]
                    ra,dec=group[ra_c],group[dec_c]
                    valid=np.flatnonzero(np.isfinite(ra)&np.isfinite(dec)&(ra>=0)&(ra<360)&(np.abs(dec)<=90))
                    if not len(valid):
                        continue
                    ii,jj,_,_=search_around_sky(anchors,SkyCoord(ra[valid],dec[valid],unit='deg'),8*u.arcsec)
                    for i,j in zip(ii,valid[jj]):
                        r=group[j]
                        plate,fiber,mjd=int(r[plate_c]),int(r[fiber_c]),int(r[col['MJD']])
                        run=str(r[col['RUN2D']]).strip()
                        ident=f'{plate:04d}-{mjd}-{fiber:04d}'
                        primary_col=resolve('SPECPRIMARY','SCIENCEPRIMARY')
                        sn_col=resolve('SN_MEDIAN_ALL','SNMEDIANALL')
                        wa=int(r[col['ZWARNING']])
                        source=str(r[col['CLASS']]).strip()
                        sn=as_float(r[sn_col]) if sn_col else None
                        # When summary SNR is unavailable, remeasure it in the coadd before any resolved inference.
                        robust=wa==0 and (sn is None or sn>=5)
                        quality=dict(plate=plate,mjd=mjd,fiber=fiber,run2d=run,zwarning=wa,
                                     sn_median_all=sn,specprimary=int(r[primary_col]) if primary_col else None,
                                     plug_coordinate_columns=[ra_c,dec_c])
                        survey='SDSS-I/II' if run.isdigit() else 'BOSS/eBOSS'
                        matches[keys[i]][ident]=dict(catalog='SDSS I-IV DR17',identifier=ident,
                           classification=source,ra=float(ra[j]),dec=float(dec[j]),
                           target_ra=float(ra[j]),target_dec=float(dec[j]),z=as_float(r[col['Z']]),
                           z_err=as_float(r[col['Z_ERR']]),quality_pass=robust,quality=quality,
                           fiber_radius_arcsec=1.5 if survey=='SDSS-I/II' else 1.,
                           position_basis='Legacy catalogue drilled fiber position; nominal observing centre')
                result={k:list(rows.values()) for k,rows in matches.items()}
                rec.update(status='complete',url=url,rows=len(tab),matched=sum(map(len,result.values())),
                           aperture='8 arcsec from any image or adopted centre',coordinate_columns=[ra_c,dec_c])
            # Free the large catalogue before storing per-target coadds.
            file.unlink()
            return result,rec
        except Exception as exc:
            rec['errors'].append(dict(url=url,error=str(exc)[:350]))
            cached=cache.get(url)
            if cached is not None and cached.is_file():
                cached.unlink()
    return result,rec


def read_boss_coadd(h):
    tabs=[ext for ext in h[1:] if isinstance(ext,fits.BinTableHDU) and
          {'loglam','flux','ivar','or_mask'}<={k.lower() for k in ext.columns.names}]
    if len(tabs)!=1:
        raise ValueError('Expected one unambiguous BOSS coadded flux table')
    co=tabs[0].data
    keys={k.lower():k for k in co.names}
    return {'BOSS':(10**np.asarray(co[keys['loglam']],dtype=float),
             np.asarray(co[keys['flux']],dtype=float).copy(),np.asarray(co[keys['ivar']],dtype=float).copy(),
             np.asarray(co[keys['or_mask']],dtype=np.uint32).copy())}


def legacy_pixels(row,cache,log):
    q=row['quality']
    run,plate,mjd,fiber=q['run2d'],q['plate'],q['mjd'],q['fiber']
    name=f'spec-{plate:04d}-{mjd}-{fiber:04d}.fits'
    tree='sdss' if run.isdigit() else 'eboss'
    versions=list(dict.fromkeys([run,'26' if tree=='sdss' else 'v5_13_2']))
    for host in ['dr17.sdss.org','data.sdss.org','dev-mirror.sdss.org']:
        for vers in versions:
            for subset in ['lite/','']:
                url=(f'https://{host}/sas/dr17/{tree}/spectro/redux/{vers}/spectra/'
                     f'{subset}{plate:04d}/{name}')
                try:
                    with fits.open(spath(cache,url,log),memmap=True) as h:
                        hdr=h[0].header
                        if any(int(hdr.get(k,expected))!=expected for k,expected in
                               [('PLATEID',plate),('MJD',mjd),('FIBERID',fiber)]):
                            raise ValueError('Legacy coadd identifier does not match')
                        pra,pdec=as_float(hdr.get('PLUG_RA')),as_float(hdr.get('PLUG_DEC'))
                        if pra is not None and pdec is not None:
                            if separation([pra,pdec],[row['ra'],row['dec']])>.6:
                                raise ValueError('Legacy plug coordinates inconsistent with catalogue')
                            row.update(ra=pra,dec=pdec,position_basis='Drilled fibre coordinates in raw coadd header')
                        bands=read_boss_coadd(h)
                        _,f,iv,m=bands['BOSS']
                        good=(m==0)&(iv>0)&np.isfinite(f)&np.isfinite(iv)
                        measured=float(np.median(f[good]*np.sqrt(iv[good]))) if good.any() else 0
                        q['coadd_median_pixel_snr']=measured
                        if q['sn_median_all'] is None:
                            row['quality_pass']=row['quality_pass'] and measured>=5
                        return bands,url
                except Exception as exc:
                    error=str(exc)[:350]
    raise ValueError('Legacy coadd retrieval failed: '+error)


def pixels(row,cache,log):
    """Validate fibre identity before recording calibrated pixels."""
    if row['catalog']=='SDSS I-IV DR17':
        return legacy_pixels(row,cache,log)
    target=row['identifier']
    q=row['quality']
    if row['catalog']=='DESI DR1':
        survey=q.get('survey','')
        program=q.get('program','')
        hp=as_float(q.get('healpix'))
        if not survey or not program or hp is None or not target.isdigit():
            raise ValueError('Missing native DESI survey/program/HEALPix/target identifier')
        hp=int(hp)
        rel=f'/spectro/redux/iron/healpix/{survey}/{program}/{hp//100}/{hp}/coadd-{survey}-{program}-{hp}.fits'
        urls=[host+rel for host in ('https://data.desi.lbl.gov/public/dr1','https://webdav-hdfs.pic.es/data/public/DESI/DR1')]
        for url in urls:
            try:
                with fits.open(spath(cache,url,log),memmap=False) as h:
                    fm=h['FIBERMAP'].data
                    ix=np.flatnonzero(fm['TARGETID']==int(target))
                    if len(ix)!=1:
                        raise ValueError('DESI coadd TARGETID missing or duplicated')
                    j=ix[0]
                    if separation([float(fm['TARGET_RA'][j]),float(fm['TARGET_DEC'][j])],
                                  [row['target_ra'],row['target_dec']])>.25:
                        raise ValueError('DESI coadd target coordinates disagree with catalogue')
                    # DESI stores unscaled wavelength, flux, inverse variance and mask per camera.
                    result={}
                    for band in 'BRZ':
                        w=np.array(h[band+'_WAVELENGTH'].data,dtype=float)
                        f=np.array(h[band+'_FLUX'].data[j],dtype=float)
                        iv=np.array(h[band+'_IVAR'].data[j],dtype=float)
                        mask=np.asarray(h[band+'_MASK'].data[j],dtype=np.uint32)
                        result[band]=(w,f,iv,mask)
                    return result,url
            except Exception as exc:
                last=exc
        raise ValueError(str(last))
    field,mjd=int(q['field']),int(q['mjd'])
    base='https://dr19.sdss.org/sas/dr19/spectro/boss/redux/v6_1_3/spectra/'
    # Directory layout differs between plate and FPS coadds in DR19.
    for family in ('lite','full'):
        for d in dict.fromkeys([f'{field:06d}',str(field),f'{field:06d}p',str(field)+'p']):
            listing=urljoin(base,f'{family}/{d}/')
            if listing not in cache:
                try:
                    r=requests.get(listing,timeout=(8,18))
                    r.raise_for_status()
                    html=r.content
                    cache[listing]=[a.get('href','') for a in BeautifulSoup(html,'html.parser').find_all('a')]
                    log.append(dict(url=listing,bytes=len(html),sha256=hashlib.sha256(html).hexdigest()))
                except requests.RequestException:
                    cache[listing]=[]
            hrefs=cache[listing]
            for folder in ([listing,urljoin(listing,f'{mjd}/')] if f'{mjd}/' in hrefs else [listing]):
                if folder!=listing and folder not in cache:
                    try:
                        r=requests.get(folder,timeout=(8,18));r.raise_for_status()
                        cache[folder]=[a.get('href','') for a in BeautifulSoup(r.content,'html.parser').find_all('a')]
                        log.append(dict(url=folder,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest()))
                    except requests.RequestException:
                        cache[folder]=[]
                for name in cache[folder]:
                    if not name.startswith('spec-') or not name.endswith('.fits'):
                        continue
                    tokens=Path(name).stem.split('-')
                    if len(tokens)<4 or tokens[-2]!=str(mjd) or not tokens[-1].isdigit() or int(tokens[-1])!=int(target):
                        continue
                    url=urljoin(folder,name)
                    with fits.open(spath(cache,url,log),memmap=True) as h:
                        return read_boss_coadd(h),url
    raise ValueError('No DR19 daily coadd found for exact field, MJD and CATALOGID')


def measure_line(bands,z,label,rest):
    center=rest*(1+z)
    trials=[]
    for band,(wave,flux,ivar,mask) in bands.items():
        v=C*(wave/center-1)
        use=(mask==0)&np.isfinite(flux)&np.isfinite(ivar)&(ivar>0)&(abs(v)<8200)
        # Suppress other expected transitions when estimating continuum/shape.
        for other,lam in LINES.items():
            if other!=label:
                use &= abs(C*(wave/(lam*(1+z))-1))>540
        sides=use&(abs(v)>4500)
        if min(np.sum(sides&(v<0)),np.sum(sides&(v>0)))<5:
            continue
        scale=5000.
        continuum=np.polyfit(v[sides]/scale,flux[sides],1,w=np.sqrt(ivar[sides]))
        residual=flux-np.polyval(continuum,v/scale)
        inner=use&(abs(v)<2600)
        if inner.sum()<12 or min(np.sum(inner&(v<0)),np.sum(inner&(v>0)))<4:
            continue
        vin,fin,weight=v[inner],residual[inner],np.sqrt(ivar[inner])
        peak=float(np.max(fin))
        if peak<=0:
            continue
        fit=least_squares(lambda p:(p[0]*np.exp(-.5*((vin-p[1])/p[2])**2)-fin)*weight,
                          [peak,0,350],bounds=([0,-1000,45],[np.inf,1000,3900]),max_nfev=250)
        amp,offset,sigma=fit.x
        template=np.exp(-.5*((vin-offset)/sigma)**2)
        # Fixed-template amplitude precision, conditional on fitted width and centroid.
        sn=float(amp*np.sqrt(np.sum((template*weight)**2)))
        # Reject a spike or unconstrained boundary fit as line-shape evidence.
        detected=bool(sn>=5 and abs(offset)<950 and sigma<3800 and
                      np.sum((abs(vin-offset)<max(2*sigma,150))&(fin*weight>1.5))>=3)
        profile=[]
        for left in range(-2400,2400,300):
            segment=inner&(v>=left)&(v<left+300)
            if segment.sum()<2:
                profile.append(None)
            else:
                w=ivar[segment]
                profile.append([float(np.average(residual[segment],weights=w)),float(1/np.sqrt(w.sum()))])
        trials.append(dict(line=label,rest_wave=rest,camera=band,detected=detected,
                           conditional_amplitude_snr=sn,velocity_offset_kms=float(offset),
                           observed_fwhm_kms=float(2.35482*sigma),profile=profile,
                           continuum_bins=int(sides.sum()),line_bins=int(inner.sum())))
    return max(trials,key=lambda v:(v['detected'],v['conditional_amplitude_snr'])) if trials else None


def compare_profiles(a,b):
    """Independent, isolated spectra only; one free positive flux scale per line."""
    result=[]
    for label in PROFILE:
        if label not in a or label not in b or not a[label]['detected'] or not b[label]['detected']:
            continue
        pairs=[(x,y) for x,y in zip(a[label]['profile'],b[label]['profile']) if x is not None and y is not None]
        if len(pairs)<9:
            continue
        f,e,g,s=np.asarray([[*x,*y] for x,y in pairs]).T
        def distance(logratio):
            scale=np.exp(logratio)
            return float(np.sum((g-scale*f)**2/(s*s+scale*scale*e*e)))
        fit=minimize_scalar(distance,bounds=(-3.9,3.9),method='bounded')
        dof=len(f)-1
        result.append(dict(line=label,relative_flux_scale=float(np.exp(fit.x)),chi2=float(fit.fun),
                           approximate_dof=dof,descriptive_p=float(chi2.sf(fit.fun,dof))))
    return sorted(result,key=lambda x:x['line'])


def assess_row(row,system):
    target=[row['target_ra'],row['target_dec']]
    fiber=[row['ra'],row['dec']]
    if None in fiber or None in target:
        return dict(association='missing_fiber_position',reliable=False)
    d=[separation(fiber,p) for p in system['images']]
    nearest=int(np.argmin(d))
    other=min(v for i,v in enumerate(d) if i!=nearest)
    center=separation(fiber,system['center'])
    robust=row['quality_pass'] and row['z'] is not None and row['z_err'] is not None
    robust=bool(robust and (row['classification']=='STAR' or 0<row['z_err']<.01*(1+row['z'])))
    if row['catalog']=='SDSS I-IV DR17' and row['quality']['sn_median_all'] is None:
        robust = robust and row['quality'].get('coadd_median_pixel_snr',0)>=5
    on_image=d[nearest]<=.50 and separation(target,system['images'][nearest])<=.65 and other>=2.5
    on_center=center<=.50 and separation(target,system['center'])<=.65 and min(d)>row['fiber_radius_arcsec']+.65
    if on_image and center>row['fiber_radius_arcsec']+.65:
        loc='isolated_image'
    elif on_center:
        loc='isolated_centre'
    elif on_image:
        loc='image_with_central_light_in_aperture'
    else:
        loc='blended_or_offset'
    return dict(association=loc,reliable=robust,nearest_image_id=system['ids'][nearest],
                nearest_image_arcsec=d[nearest],other_image_arcsec=other,center_arcsec=center)


def conclusions(records):
    good=[r for r in records if r['association']['reliable']]
    image=[r for r in good if r['association']['association']=='isolated_image']
    center=[r for r in good if r['association']['association']=='isolated_centre']
    outcome=[]
    pairs=[]
    for a,b in combinations(image,2):
        if a['association']['nearest_image_id']==b['association']['nearest_image_id']:
            continue
        delta=abs(a['z']-b['z'])
        err=np.hypot(a['z_err'] or 0,b['z_err'] or 0)
        mean=(a['z']+b['z'])/2
        dz_reject=max(.01*(1+mean),5*err)
        close=bool(delta<=max(.002*(1+mean),3*err))
        cats=[a['classification'],b['classification']]
        row=dict(image_ids=[a['association']['nearest_image_id'],b['association']['nearest_image_id']],
                 spectroscopic_ids=[a['identifier'],b['identifier']],delta_z=delta,classes=cats,
                 redshift_conflict_threshold=dz_reject,redshifts_agree=close)
        if 'STAR' in cats:
            outcome.append('stellar_contaminant_on_isolated_image')
        elif delta>dz_reject:
            outcome.append('different_redshifts_on_two_isolated_images')
        elif close and cats==['QSO','QSO']:
            row['profile_comparisons']=compare_profiles(a['lines'],b['lines'])
            outcome.append('two_separate_quasar_spectra_at_consistent_redshift_physical_pair_remains_possible')
        elif close and all(c=='GALAXY' for c in cats) and mean<.25:
            outcome.append('two_low_redshift_galaxy_spectra_consistent_with_nuclei_or_emitting_knots')
        elif close:
            outcome.append('distinct_components_at_similar_redshift_with_inconclusive_spectral_classes')
        else:
            outcome.append('distinct_components_with_marginal_redshift_difference')
        pairs.append(row)
    for a in center:
        for b in image:
            if a['classification']=='GALAXY' and b['classification']=='QSO' and b['z']-a['z']>max(.01,5*np.hypot(a['z_err'],b['z_err'])):
                outcome.append('lower_redshift_central_galaxy_and_background_quasar_lensing_unresolved')
    if any(r['classification']=='STAR' for r in image):
        outcome.append('stellar_contaminant_on_isolated_image')
    if not pairs and not outcome:
        low=[r for r in good if r['classification'] in ('GALAXY','QSO') and r['z'] is not None and 0<r['z']<.25]
        if low:
            outcome.append('low_redshift_spectral_emission_present_second_quasar_image_unestablished')
        elif good:
            outcome.append('catalogue_spectrum_present_no_independent_comparison_of_images')
        elif records:
            outcome.append('only_low_quality_or_unassigned_spectra')
        else:
            outcome.append('no_spectral_observation_in_the_search_apertures')
    return dict(outcomes=sorted(set(outcome)),distinct_isolated_image_ids=sorted({r['association']['nearest_image_id'] for r in image}),
                independent_image_comparisons=pairs,raw_spectra_read=sum(r['pixel_status']=='complete' for r in records),
                line_measurement_caution='Observed widths use a single Gaussian and adjacent-line masks. Complexes may blend. Agreement cannot distinguish a physical quasar pair; microlensing and variability can change profiles.')


def run(evidence,followup):
    evidence,followup=Path(evidence),Path(followup)
    out=followup/'resolved_spectra'
    out.mkdir(parents=True,exist_ok=True)
    systems=targets(evidence,followup)
    sdss,sdss_receipt=sdss_centres({k:s['center'] for k,s in systems.items()},evidence)
    log=[]
    summaries=[]
    with tempfile.TemporaryDirectory() as tmp:
        cache={'dir':Path(tmp)}
        legacy,legacy_receipt=legacy_spectra(systems,cache,log)
        for ident,system in systems.items():
            folder=out/ident
            folder.mkdir(exist_ok=True)
            desi,receipt=desi_spectra(system['center'],system['images'],folder,evidence)
            records=[]
            for i,row in enumerate(desi+sdss.get(ident,[])+legacy.get(ident,[])):
                assoc=assess_row(row,system)
                item=dict(row,association=assoc,lines={},pixel_status='quality_or_aperture_not_selected')
                # Inspect even a blended high-quality spectrum if it can illuminate a nucleus.
                if row['quality_pass'] and None not in [row['ra'],row['dec']] and (
                    min(assoc.get('nearest_image_arcsec',100),assoc.get('center_arcsec',100))<=2.0):
                    try:
                        bands,url=pixels(row,cache,log)
                        item.update(row,association=assess_row(row,system))
                        traces=[]
                        for name,(w,f,iv,m) in bands.items():
                            traces.append(pd.DataFrame(dict(camera=name,wavelength_vacuum_A=w,flux=f,ivar=iv,mask=m)))
                        path=folder/f'coadd-{i}.csv.gz'
                        pd.concat(traces).to_csv(path,index=False,float_format='%.8g')
                        if row['z'] is not None and row['z']>0:
                            for label,rest in LINES.items():
                                line=measure_line(bands,row['z'],label,rest)
                                if line:
                                    item['lines'][label]=line
                        item.update(pixel_status='complete',coadd_url=url,pixels=path.name,pixel_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                    except Exception as exc:
                        item.update(pixel_status='unavailable',pixel_error=str(exc)[:350])
                records.append(item)
            result=conclusions(records)
            if not records and (not receipt.get('status','').startswith('complete') or sdss_receipt.get('status')!='complete' or legacy_receipt.get('status')!='complete'):
                result['outcomes']=['no_returned_spectrum_incomplete_archive_coverage']
            save(folder/'spectral_evidence.json',dict(system_id=ident,**system,desi_receipt=receipt,sdss_status=sdss_receipt['status'],legacy_status=legacy_receipt['status'],
                                                       spectra=records,result=result))
            summaries.append(dict(system_id=ident,aliases=system['aliases'],**result))
    save(out/'summary.json',dict(examined_systems=len(systems),systems=summaries,sdss_receipt=sdss_receipt,legacy_receipt=legacy_receipt))
    save(out/'archive_receipts.json',log)
    print('Dedicated spectral fields',len(summaries),'coadds read',sum(r['raw_spectra_read'] for r in summaries),flush=True)


if __name__=='__main__':
    run(*sys.argv[1:])
