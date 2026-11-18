"""Derive conditional infrared angular sizes at two KMT event positions."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

from sightlines import HERE,CFA,read,catalog,angle,urlencode
from source_size import clump

TARGETS = ('KMT-2024-BLG-0824','KMT-2024-BLG-3278')


def tsv(b):
    lines=[l for l in b.decode().splitlines() if l.strip() and not l.startswith('#')]
    if not lines or not lines[0].startswith('_r\t'):
        raise ValueError('invalid VizieR field')
    return [dict(zip(lines[0].split('\t'),[x.strip() for x in l.split('\t')])) for l in lines[3:]]


def field(event, cache, ra, dec):
    if event.endswith('0824'):
        name='virac-clump';k='Ksmag';j='Jmag'
        params={'-source':'II/387/virac2','-c':f'{ra:.9f} {dec:.9f}','-c.rs':120,'Ksmag':'11..15',
                '-out.add':'_r','-sort':'_r','-out.max':30000,'-out':'srcid,Jmag,Ksmag,o_Ksmag,MedKsmag'}
    else:
        name='sbc-clump';k='Ksmag3';j='Jmag3'
        params={'-source':'II/348/vvv2','-c':f'{ra:.9f} {dec:+.9f}','-c.rs':120,'-out.add':'_r','-sort':'_r',
                '-out.max':30000,'-out':'srcid,raJ2000,deJ2000,Jmag3,Hmag3,Ksmag3,Jperrbits,Ksperrbits'}
    buf,receipt=read(cache,name,CFA+urlencode(params));rows=tsv(buf)
    if event.endswith('0824'):
        rows=[r for r in rows if r[j] and r[k] and int(r['o_Ksmag'])>=8 and float(r['MedKsmag'])<.15]
    else:
        rows=[r for r in rows if r[j] and r[k] and int(r['Jperrbits'])<256 and int(r['Ksperrbits'])<256]
    a=np.array([[float(r['_r']),float(r[j]),float(r[k])] for r in rows])
    return dict(receipt=receipt,n=len(rows),clumps=clump(a))


def radius(j,k,ej,ek,reddening):
    # Kervella et al. 2004 Table 4, C0=K, C1=J. The calibration is for
    # dereddened Johnson magnitudes and dwarfs/subgiants; retain only
    # plausible intrinsic model colors in the approximation below.
    vals=[]
    color=np.array([1.,-1.,-1.,1.,0.])
    bounds=[(j-ej,j+ej),(k-ek,k+ek),(min(reddening),max(reddening)),(-.10,.10),(-.05,.05)]
    for coeff in (.428,.528):
        c=.3192*color-.2*np.array([0.,1.,-coeff,0.,1.])
        for sign,scatter in [(1.,1-.0587),(-1.,1+.0587)]:
            fit=linprog(sign*c,A_ub=[color,-color],b_ub=[.85,.02],bounds=bounds)
            if fit.success:
                vals.append(500*10**(.5089+c@fit.x)*scatter)
    return [round(min(vals),4),round(max(vals),4)] if vals else None


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args();all_rows=json.loads((HERE/'sightline-receipts.json').read_text())
    results=[]
    for event in TARGETS:
        r=next(x for x in all_rows if x['event']==event)
        ra,dec=r['ra_deg'],r['dec_deg']
        c=args.cache/event
        fld=field(event,c,ra,dec)
        pars={'-source':'II/387/virac2','-c':f'{ra:.9f} {dec:.9f}','-c.rs':3,'-out.add':'_r','-sort':'_r',
              '-out.max':50,'-out':'srcid,Jmag,e_Jmag,o_Jmag,Ksmag,e_Ksmag,o_Ksmag'}
        url=CFA+urlencode(pars)
        b,rec=read(c,'errors',url)
        if b'#INFO\tError=' in b:
            (c/'errors.receipt').unlink()
            b,rec=read(c,'errors',url)
        errs={d['srcid']:d for d in tsv(b)}
        excess=[d['J_K_centroid']-.68 for d in fld['clumps']]
        e0=float(np.median(excess));e=[max(0,min(excess)-.05),max(excess)+.05]
        src=[]
        for x in r['catalogs'][1]['rows']:
            if float(x['_r'])>1.0:continue
            j,k=float(x['Jmag']),float(x['Ksmag'])
            er=errs[x['srcid']]
            # Use the scatter among individual epochs, not the much smaller standard error
            # from treating all repeated exposures as independent. Numbers are sensitivity ranges.
            ej,ek=max(.03,float(er['e_Jmag'])),max(.03,float(er['e_Ksmag']))
            color=j-k-e0
            nominal=500*10**(.3192*color+.5089-.2*(k-.478*e0))
            runit=radius(j,k,ej,ek,e)
            front=radius(j,k,ej,ek,[0,e[0],e[1]])
            rho_values=[2.0611,2.5973] if event.endswith('0824') else [2.5707,2.5820]
            src.append(dict(srcid=x['srcid'],sep_arcsec=float(x['_r']),J=j,Ks=k,J_K_0_full_screen=round(color,3),
                J_epoch_sd=ej,K_epoch_sd=ek,nJ=er['o_Jmag'],nK=er['o_Ksmag'],theta_star_uas_nominal=round(nominal,4),
                full_clump_screen_radius_uas=runit,foreground_to_clump_radius_uas=front,
                tested_rho_geometry=rho_values,thetaE_at_full_screen_nominal_uas=[round(nominal/v,4) for v in rho_values],
                thetaE_full_screen_sensitivity_uas=[round(runit[0]/max(rho_values),4),round(runit[1]/min(rho_values),4)],
                note='Conditional on this isolated dwarf/subgiant supplying all of the catalog J and Ks source light.'))
        results.append(dict(event=event,field=fld,error_receipt=rec,intrinsic_RC_JK=.68,observed_clump_EJK=e0,
                            full_clump_reddening_spread=e,source_options=src))
    (HERE/'angular.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
