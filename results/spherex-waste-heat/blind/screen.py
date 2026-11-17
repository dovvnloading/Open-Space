"""Photometric triage of a fixed Gaia/AllWISE field; rerun from local inputs."""
import csv
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

from catalogs import HERE, rows


def number(row, key):
    return float(row.get(key) or np.nan)


def save_csv(name, records):
    with (HERE/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)


def verify():
    receipts = json.loads((HERE/'catalog_receipts.json').read_text())
    for name, rec in receipts.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != rec['sha256']:
            raise ValueError(f'Input digest mismatch: {name}')


def load():
    stars = rows('stars.csv')
    matches = {r['source_id']: r for p in sorted(HERE.glob('xm_[0-9]*.csv')) for r in rows(p)}
    wise = {r['designation']: r for p in sorted(HERE.glob('wise_[0-9]*.csv')) for r in rows(p)}
    if len(wise) != len({r['original_ext_source_id'] for r in matches.values()}):
        raise ValueError('Incomplete AllWISE photometry')
    for s in stars:
        cross = matches.get(s['source_id'], {})
        s.update(cross)
        phot = wise.get(cross.get('original_ext_source_id'), {})
        # Keep Gaia astrometry; WISE has its own two position columns.
        s.update({key:v for key,v in phot.items() if key not in ['ra','dec']})
    return stars


def quality(s):
    color = number(s,'bp_rp')
    mg = number(s,'phot_g_mean_mag') + 5*np.log10(number(s,'parallax'))-10
    ridge = np.interp(color, np.arange(.5,3.51,.5), [2.5,4.6,6.4,8.2,9.5,10.8,12.3])
    if not (abs(mg-ridge)<=1.4 and number(s,'ruwe')<1.4 and s['duplicated_source']=='false'):
        return 'Gaia astrometry or CMD'
    if not s.get('original_ext_source_id'):
        return 'no AllWISE match'
    if not (number(s,'number_of_neighbours')==1 and number(s,'number_of_mates')==0):
        return 'nonunique match'
    if not s.get('designation'):
        raise ValueError('Missing requested AllWISE row')
    if not (number(s,'ext_flg')==0 and number(s,'nb')==1 and s['cc_flags'][0]=='0'
            and number(s,'w1snr')>=10 and number(s,'w1mpro')>=8 and number(s,'w1sat')==0
            and 0<number(s,'w1rchi2')<=2 and number(s,'w1sigmpro')>0):
        return 'W1 or morphology quality'
    return 'eligible'


def band_info(s, band, zp):
    prefix = f'w{band}'
    flux, err = (number(s, prefix+field)*10**(-0.4*zp) for field in ['flux','sigflux'])
    valid = (s.get('cc_flags','0000')[band-1:band]=='0' and number(s,prefix+'sat')==0
             and np.isfinite(flux) and err>0 and np.isfinite(err))
    chi = number(s,prefix+'rchi2')
    # An undetected source has no profile chi square. Its signed measurement still
    # constrains the photospheric locus and can receive a flux-level injection.
    if np.isfinite(chi) and not 0<chi<=2:
        valid = False
    if flux>=2*err and not np.isfinite(chi):
        valid = False
    return flux,err,valid


class Screen:
    def __init__(self, stars, zp):
        self.stars=stars
        self.ids=np.array([s['source_id'] for s in stars])
        self.bp=np.array([number(s,'bp_rp') for s in stars])
        self.mag=np.array([number(s,'w1mpro') for s in stars])
        self.f1=10**(-.4*self.mag)
        self.fe1=np.array([number(s,'w1sigmpro')*np.log(10)/2.5 for s in stars])
        self.band={}
        for b in (3,4):
            triplets=[band_info(s,b,zp[b]) for s in stars]
            flux,err,valid=map(np.array,zip(*triplets))
            self.band[b]=(flux,err,valid)

    def model(self, i, b, w1):
        flux,err,valid=self.band[b]
        nearby=(np.abs(self.bp-self.bp[i])<=.4) & (np.abs(self.mag-w1)<=1.5) & valid
        nearby[i]=False
        ix=np.flatnonzero(nearby)
        scale=((self.bp[ix]-self.bp[i])/.4)**2+((self.mag[ix]-w1)/1.5)**2
        ix=ix[np.argsort(scale,kind='stable')[:25]]
        if len(ix)<5: return None
        q=flux[ix]/self.f1[ix]
        center=np.median(q)
        if center<=0: return None
        scat=1.4826*np.median(np.abs(q-center))
        err0=np.median(np.hypot(err[ix]/self.f1[ix],self.fe1[ix]*center))
        intrinsic=np.sqrt(max((.03*center)**2,scat**2-err0**2))
        loc_err=1.253*np.hypot(intrinsic,err0)/np.sqrt(len(ix))
        return center,intrinsic,loc_err,ix

    def score(self,i,b,flux1,fluxb):
        actual,err,valid=self.band[b]
        w1=-2.5*np.log10(flux1)
        if w1<8.-1e-8 or not valid[i]: return None
        model=self.model(i,b,w1)
        if model is None: return None
        center,intrinsic,loc_err,ix=model
        expected=flux1*center
        sigma=np.sqrt(err[i]**2+flux1**2*(intrinsic**2+loc_err**2+(center*self.fe1[i])**2))
        z=(fluxb-expected)/sigma
        ratio=fluxb/expected
        snr=fluxb/err[i]
        longmag=-2.5*np.log10(fluxb) if fluxb>0 else np.inf
        detect=(snr>=5 and z>=5 and ratio>=10**.1 and longmag>={3:3.8,4:-.4}[b])
        return dict(z=z,flux_over_model=ratio,snr=snr,excess=bool(detect),n=len(ix),
                    expected=expected,ix=ix)


def zero_points(stars):
    zps,diag={},{}
    for b in (3,4):
        p=f'w{b}'
        z=[number(s,p+'mpro')+2.5*np.log10(number(s,p+'flux')) for s in stars
           if number(s,p+'flux')>5*number(s,p+'sigflux')>0 and number(s,p+'sigmpro')>0
           and number(s,p+'sat')==0]
        med=float(np.median(z)); dev=float(np.max(np.abs(np.array(z)-med)))
        if len(z)<5 or dev>.03: raise ValueError(f'Unstable instrumental zero point W{b}: n={len(z)} range={dev}')
        zps[b]=med
        diag[str(b)]=dict(zp=med,max_abs_deviation_mag=dev,n=len(z))
    return zps,diag


@lru_cache(maxsize=2048)
def inband(b, temp):
    wave,r,_=np.loadtxt(HERE/f'RSR-W{b}.txt').T
    return np.trapezoid(r / wave**4 / np.expm1(14387.76877/(wave*temp)),wave)


@lru_cache(maxsize=3)
def zero_mag_counts(b):
    # Wright et al. (2010), F_lambda at the WISE isophotal wavelengths.
    # A constant F_lambda with this value gives the Vega-system zero count.
    wave, response, _ = np.loadtxt(HERE/f'RSR-W{b}.txt').T
    f_lambda_zero = {1: 8.180e-15, 3: 6.515e-17, 4: 5.090e-18}
    return f_lambda_zero[b]*np.trapezoid(wave*response,wave)


def main():
    verify()
    stars=load()
    first={s['source_id']:quality(s) for s in stars}
    eligible=[s for s in stars if first[s['source_id']]=='eligible']
    zp,diag=zero_points(stars)
    phot=Screen(eligible,zp)
    scanned=[]; models={}; controls=[]; candidates=[]
    out={s['source_id']:dict(source_id=s['source_id'],ra=s['ra'],dec=s['dec'],
         allwise=s.get('designation',''),status=first[s['source_id']],z3='',z4='',
         ratio3='',ratio4='',snr3='',snr4='',selected3='',selected4='') for s in stars}
    ctrl_rows=[]
    for i,s in enumerate(eligible):
        caught=False; scores=[]
        for b,(flux,err,valid) in phot.band.items():
            test=phot.score(i,b,phot.f1[i],flux[i])
            if test is None: continue
            models[i,b]=test
            for col,key in [('z','z'),('ratio','flux_over_model'),('snr','snr')]:
                out[s['source_id']][col+str(b)]=round(test[key],5)
            out[s['source_id']]['selected'+str(b)]=int(test['excess'])
            ctrl_rows.append(dict(source_id=s['source_id'],band=b,
                                  comparison_ids=' '.join(phot.ids[test['ix']])))
            caught |= test['excess']
            scores.append(test['z'])
        if scores: scanned.append(i)
        if caught:
            candidates.append(s['source_id'])
            out[s['source_id']]['status']='photometric excess'
        elif scores:
            out[s['source_id']]['status']='screened'
            if max(abs(x) for x in scores)<3: controls.append(i)
        else:
            out[s['source_id']]['status']='no usable long-band controls/measurement'
    save_csv('screen.csv',list(out.values()))
    save_csv('comparison_ids.csv',ctrl_rows)
    temps={}
    for i in controls:
        t=number(eligible[i],'teff_gspphot')
        if np.isfinite(t) and 2500<=t<=8000:
            temps[i]=round(t/10)*10
    injection=[]
    for dust in [100,200,300,400,500,600]:
        factors={(i,b):(t/dust)**4 * inband(b,dust)/inband(1,t)
                               * zero_mag_counts(1)/zero_mag_counts(b)
                 for i,t in temps.items() for b in (1,3,4)}
        for gamma in [0.,.0001,.0003,.001,.003,.01,.03,.1]:
            detects=Counter()
            for i in temps:
                flux1=phot.f1[i]*(1+gamma*factors[i,1]); hits=[]
                for b,(flux,err,valid) in phot.band.items():
                    if (i,b) not in models: continue
                    fluxb=flux[i]+gamma*factors[i,b]*phot.f1[i]
                    test=phot.score(i,b,flux1,fluxb)
                    if test and test['excess']: hits.append(b)
                for b in hits: detects[b]+=1
                detects['either']+=int(bool(hits))
            injection.append(dict(dust_K=dust,L_excess_over_L_star=gamma,controls=len(temps),
                                  recovered_w3=detects[3],recovered_w4=detects[4],
                                  recovered_either=detects['either'],
                                  recovery_fraction=round(detects['either']/len(temps),5)))
    save_csv('injections.csv',injection)
    save_csv('injection_controls.csv',[dict(source_id=phot.ids[i],teff_K=t) for i,t in temps.items()])
    counts=dict(field=len(stars),first_failure_counts=dict(Counter(first.values())),
                catalog_zp_audit=diag,screened_either=len(scanned),
                screened_w3=sum(b==3 for _,b in models),screened_w4=sum(b==4 for _,b in models),
                photometric_candidates=candidates,within_3sigma_controls=len(controls),
                injection_controls=len(temps))
    (HERE/'counts.json').write_text(json.dumps(counts,indent=2)+'\n')
    print(json.dumps(counts,indent=2))


if __name__=='__main__':
    main()
