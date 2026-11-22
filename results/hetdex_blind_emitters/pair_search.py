"""Blind Ly-alpha/He II pair search, including displaced-pair and injection controls."""
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path

import numpy as np
from scipy.ndimage import convolve1d, maximum_filter1d, uniform_filter1d

import spectral as sp

RATIO = sp.UV["heii"] / sp.UV["lya"]
OFFSETS = (0, -76., 76.)
WIDTHS = (1.8, 2.4, 3.5, 5.0)
METALS = ("civ1","civ2","oiii_uv1","oiii_uv2","ciii1","ciii2")
VEL = (-300.,600.)
SEED = 482319


def jsonline(obj):
    from vetting import atom
    return json.dumps(atom(obj),allow_nan=False,separators=(",",":"))


def conv_score(y,e,dx):
    finite=np.isfinite(y)&np.isfinite(e)&(e>0)
    weights=np.divide(1.,e**2,out=np.zeros_like(e),where=finite)
    yw=np.where(finite,y,0.)*weights
    den=uniform_filter1d(weights,55,axis=1,mode="constant")
    base=np.divide(uniform_filter1d(yw,55,axis=1,mode="constant"),den,
                   out=np.zeros_like(yw),where=den>0)
    residual=yw-base*weights
    score=np.zeros_like(yw)
    for sigma in (2.2,4.2):
        offsets=np.arange(-math.ceil(3.6*sigma/dx),math.ceil(3.6*sigma/dx)+1)*dx
        kernel=np.exp(-.5*(offsets/sigma)**2)
        num=convolve1d(residual,kernel,axis=1,mode="constant")
        norm=np.sqrt(convolve1d(weights,kernel**2,axis=1,mode="constant"))
        sn=np.divide(num,norm,out=np.zeros_like(yw),where=norm>0)
        np.maximum(score,sn,out=score)
    return np.where(finite,score,0.)


def centroid(w,score,j):
    if j==0 or j+1==len(w):
        return float(w[j])
    left,mid,right=score[j-1:j+2]
    curvature=left-2*mid+right
    sub=np.clip((left-right)/(2*curvature),-.5,.5) if curvature<-1e-8 else 0.
    return float(w[j]+(w[1]-w[0])*sub)


def local_fit(w,y,e,center):
    side=sp._local_side(w,y,e,center,5.0)
    if side is None:
        return None
    pairs=[]
    for width in WIDTHS:
        measurement=sp._profile(w,y,e,center,width,side,w[1]-w[0])
        if measurement is not None:
            area,error=measurement
            pairs.append((area/error,area,error,width))
    if not pairs:
        return None
    sn,area,err,sigma=max(pairs)
    sl=sp._slice(w,center-2.5*sigma,center+2.5*sigma)
    positive=np.maximum(y[sl]-side[0],0)*np.isfinite(e[sl])
    single=positive.sum()<=0 or positive.max()>.68*positive.sum()
    sky=np.min(abs(sp.SKY_AIR-center))<=max(6.,2*sigma)
    badcal=any(a-2*sigma<=center<=b+2*sigma for a,b in sp.BAD_AIR)
    return dict(wave=center,wave_vac=float(sp.air_to_vac(center)),flux=area*1e-17,
                flux_err=err*1e-17,snr=sn,sigma_A=sigma,sky_near=bool(sky),
                badcal=bool(badcal),single_pixel=bool(single),
                continuum=side[0]*1e-17,continuum_err=side[1]*1e-17)


def predicted(anchor_vac,offset,lo,hi):
    # Use a cyclic displacement only for controls, keeping the true redshift
    # hypothesis unshifted. Require the same anchor domain in all streams.
    expected=anchor_vac*RATIO
    if offset:
        return lo+np.mod(expected-lo+offset,hi-lo)
    return expected


def noise_scale(y,e):
    weight=np.divide(1.,e**2,out=np.zeros_like(e),where=np.isfinite(e))
    den=uniform_filter1d(weight,55,mode="constant")
    mean=np.divide(uniform_filter1d(y*weight,55,mode="constant"),den,
                   out=np.zeros_like(y),where=den>0)
    signed=((y-mean)/e)[np.isfinite(e)]
    return max(1.,float(np.median(abs(signed-np.median(signed)))*1.4826)) if signed.size else 1.


def evaluate(w,y,e,scale,ly,he,offset,target,stats):
    if ly is None or he is None:
        return None,"no_profile"
    if ly["snr"]<7:
        return None,"lya_below_7"
    if he["snr"]<4.5 or he["snr"]>=7:
        return None,"heii_outside_4p5_to_7"
    z=he["wave_vac"]/sp.UV["heii"]-1
    velocity=(target/he["wave_vac"]-1)*sp.C
    if not 1.9<=z<=2.4 or not VEL[0]<=velocity<=VEL[1]:
        return None,"redshift_or_velocity"
    for line in (he,ly):
        if line["single_pixel"] or line["sky_near"] or line["badcal"]:
            return None,"single_pixel_sky_or_badcal"
    if any(2.355*l["sigma_A"]*sp.C/l["wave"]>=900 for l in (he,ly)):
        return None,"width_above_900_kms"
    ratio=he["flux"]/ly["flux"]
    if ratio<.10:
        return None,"heii_to_lya_below_0p1"
    cont=he["continuum"]
    cerr=he["continuum_err"]
    ew=he["flux"]/(cont*(1+z)) if cont>0 else None
    ewl=max(0.,he["flux"]-1.96*he["flux_err"])/((max(cont,0.)+1.96*cerr)*(1+z))
    if cont>=3*cerr and ew<3.:
        return None,"heii_EW_below_3A_measured_continuum"
    metals={k:sp._force(w,y,e,w[1]-w[0],sp.UV[k]*(1+z),he["sigma_A"]) for k in METALS}
    if any(m.get("snr",-np.inf)>=3 for m in metals.values()):
        return None,"metal_at_or_above_3sigma"
    return dict(he,channel="heii",stream="weak_heii_pair",noise_scale=scale,
                z_heii=z,z_lya=he["wave_vac"]/sp.UV["lya"]-1,
                pair_lya_z=ly["wave_vac"]/sp.UV["lya"]-1,pair_lya=ly,
                pair_velocity_kms=velocity,pair_heii_to_lya=ratio,
                pair_offset_A=offset,pair_metals=metals,
                heii_rest_ew_A=ew,heii_rest_ew_lower_A=ewl,
                ew0=None,ew0_lower=None), "passed"


def scan_block(w,flux,error,start,stats,offsets=OFFSETS):
    dx=w[1]-w[0]
    vac=sp.air_to_vac(w)
    lo=2.9*sp.UV["heii"]
    hi=min(3.4*sp.UV["heii"],vac[-1])
    anchor_min,anchor_max=lo*(1+VEL[0]/sp.C),hi*(1+VEL[1]/sp.C)
    ly_domain=(vac*RATIO>=anchor_min)&(vac*RATIO<=anchor_max)
    blue=np.flatnonzero((vac*RATIO>=anchor_min-90)&(vac*RATIO<=anchor_max+90))
    red=np.flatnonzero((vac>=lo-85)&(vac<=hi+85))
    def take(indices,rows):
        # Use contiguous slices to preserve C layout for direct convolution.
        s=slice(int(indices[0]),int(indices[-1])+1)
        f=np.ascontiguousarray(flux[rows,s])
        e=np.array(error[rows,s],copy=True,order="C")
        for a,b in sp.BAD_AIR:
            e[:,(w[s]>=a)&(w[s]<=b)]=np.inf
        return conv_score(f,e,dx)
    bs=take(blue,slice(None))
    bl=(bs>=5.15)&(maximum_filter1d(bs,5,axis=1)==bs)&ly_domain[blue]
    zi,bi=np.where(bl)
    stats["total"]["observation_rows"]+=len(flux)
    stats["total"]["ly_pre_peaks"]+=len(bi)
    stats["total"]["ly_pre_rows"]+=len(np.unique(zi))
    if len(bi)==0:
        return []
    rows,inv=np.unique(zi,return_inverse=True)
    rs=take(red,rows)
    rm=(rs>=3.45)&(maximum_filter1d(rs,5,axis=1)==rs)
    # Catch boundary peaks spectrally narrower than the integer-grid match.
    rm&=(vac[red]>=lo-2)&(vac[red]<=hi+2)
    anchor=vac[blue[bi]]
    # Advanced indexing into 13 fixed nearby red pixels avoids per-anchor
    # profile fitting on spectra without a possible secondary line.
    steps=np.arange(-8,6)
    measurements={}
    matched={}
    passed=[]
    for offset in offsets:
        counts=stats[str(offset)]
        target=predicted(anchor,offset,anchor_min,anchor_max)
        mid=np.searchsorted(vac[red],target)
        idx=mid[:,None]+steps
        valid=(idx>=0)&(idx<len(red))
        idx=np.clip(idx,0,len(red)-1)
        v=(target[:,None]/vac[red[idx]]-1)*sp.C
        possible=valid & (v>=VEL[0]-185)&(v<=VEL[1]+185)
        usable=possible & np.isfinite(error[zi[:,None],red[idx]]) & (error[zi[:,None],red[idx]]>0)
        for a,b in sp.BAD_AIR:
            usable&=(w[red[idx]]<a)|(w[red[idx]]>b)
        counts["anchor_opportunities"]+=int(np.count_nonzero(usable.any(axis=1)))
        ai,di=np.where(usable & rm[inv[:,None],idx])
        counts["pairs_prescreened"]+=len(ai)
        for a,d in zip(ai,di):
            row=int(zi[a])
            hj=int(idx[a,d]); bj=int(bi[a])
            if row not in measurements:
                yy,ee,_=sp._inputs(w,flux[row],error[row],1e-17)
                scale=noise_scale(yy[0],ee[0])
                measurements[row]=(yy[0],ee[0]*scale,scale)
            yy,ee,scale=measurements[row]
            lkey=(row,"ly",bj); hkey=(row,"he",hj)
            if lkey not in matched:
                center=centroid(w[blue],bs[row],bj)
                matched[lkey]=local_fit(w,yy,ee,center)
            if hkey not in matched:
                center=centroid(w[red],rs[inv[a]],hj)
                matched[hkey]=local_fit(w,yy,ee,center)
            ly,he=matched[lkey],matched[hkey]
            refined=predicted(ly["wave_vac"],offset,anchor_min,anchor_max) if ly else np.nan
            rec,reason=evaluate(w,yy,ee,scale,ly,he,offset,refined,counts)
            counts[reason]+=1
            if rec:
                rec["row_index"]=start+row
                passed.append(rec)
    return passed


def inject(w,f,e,rng):
    """One random Gaussian pair per spectrum; keep raw errors and underlying continuum."""
    z=rng.uniform(1.9,2.4)
    v=rng.uniform(*VEL)
    he=sp.vac_to_air(sp.UV["heii"]*(1+z))
    ly=sp.vac_to_air(sp.UV["lya"]*(1+z)*(1+v/sp.C))
    width=float(rng.choice((2.,3.,4.)))
    requested_he=float(rng.choice((4.5,5.5,6.5,8.)))
    requested_ly=float(rng.choice((9.,14.,22.)))
    y,err,_=sp._inputs(w,f,e,1e-17)
    scale=noise_scale(y[0],err[0])
    err=err[0]*scale
    data=f.copy()
    truth=dict(z=z,velocity_kms=v,wave_heii=he,wave_lya=ly,sigma_A=width,
               nominal_sn_heii=requested_he,nominal_sn_lya=requested_ly)
    for key,center,snr in (("heii",he,requested_he),("lya",ly,requested_ly)):
        if center<w[0] or center>w[-1]:
            return None,dict(truth,status=key+"_outside_band")
        side=sp._local_side(w,y[0],err,center,5.0)
        profile=sp._profile(w,y[0],err,center,width,side,w[1]-w[0])
        if profile is None:
            return None,dict(truth,status=key+"_insufficient_pixels")
        area=snr*profile[1]
        truth["flux_"+key]=float(area*1e-17)
        data+=area*np.exp(-.5*((w-center)/width)**2)/(np.sqrt(2*np.pi)*width)
    truth["flux_ratio"]=truth["flux_heii"]/truth["flux_lya"]
    return data,truth


def known_eight(cat,observed):
    from vetting import PRIOR
    w=cat.wave
    # Position or original detection id includes repeats in either data release.
    ra=np.asarray(cat.info["ra"],dtype=float)
    dec=np.asarray(cat.info["dec"],dtype=float)
    ids=np.asarray(cat.info["detectid"])
    byrow=defaultdict(list)
    for p in observed:
        byrow[p["row_index"]].append(p)
    catalog=[]
    for ident,cra,cdec,z in PRIOR:
        dra=((ra-cra+180)%360-180)*np.cos(np.deg2rad(cdec))
        sep=np.hypot(dra,dec-cdec)*3600.
        rows=np.flatnonzero((sep<=3.)|(ids==ident))
        entry=dict(prior_detectid=ident,ra=cra,dec=cdec,z_reference=z,observations=[])
        for row in rows:
            f,e=cat.read(int(row),int(row)+1)
            y,err,_=sp._inputs(w,f,e,1e-17)
            y,err=y[0],err[0]
            scale=noise_scale(y,err)
            err*=scale
            anchors={}
            for label in ("heii","lya"):
                wave=sp.vac_to_air(sp.UV[label]*(1+z))
                # Reference redshifts tabulated to three decimal places. Allow
                # Ly-alpha escape offsets; these maxima are diagnostic only.
                deltav=(-200,-100,0,100,200) if label=="heii" else (-300,-150,0,150,300,450,600)
                trials=[local_fit(w,y,err,wave*(1+v/sp.C)) for v in deltav]
                anchors[label]=max((x for x in trials if x),key=lambda d:d["snr"],default=None)
            he,ly=anchors["heii"],anchors["lya"]
            _,reason=evaluate(w,y,err,scale,ly,he,0,ly["wave_vac"]*RATIO if ly else np.nan,{})
            recovered=any(abs(p["z_heii"]-z)<.006 for p in byrow[int(row)])
            entry["observations"].append(dict(row_index=int(row),shotid=int(cat.shot_ids[row]),
                separation_arcsec=float(sep[row]),heii=he,lya=ly,forced_stage=reason,
                weak_pair_blind_recovery=recovered))
        catalog.append(entry)
    return catalog


def pair_search(cat,out):
    """Return all observed screening survivors for the same independent vetting as ≥7σ peaks."""
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    w=cat.wave
    stats=defaultdict(Counter)
    observed=[]
    shuffle_stats=defaultdict(Counter)
    trials=[]
    # Deterministic stratified injection positions. Flux and statistical errors
    # come from the actual source, retaining masks and pre-existing emission.
    rng=np.random.default_rng(SEED)
    inrows={int(rng.integers(a,min(a+1536,cat.size))) for a in range(0,cat.size,1536)}
    with gzip.open(out/"weak-heii-pairs.jsonl.gz","wt") as ledger:
        for start,f,e in cat.blocks(512):
            recs=scan_block(w,f,e,start,stats)
            for rec in recs:
                row=rec["row_index"]
                meta=cat.metadata(row)
                rec.update({k:meta.get(k) for k in ("ra","dec","ra_det","dec_det","detectid","source_id","shotid")})
                ledger.write(jsonline(rec)+"\n")
                if rec["pair_offset_A"]==0:
                    observed.append(rec)
            # Apply the identical physical pair search to randomly selected
            # blocks of wavelength-shuffled flux/error fibres, preserving
            # wavelength-specific noise, calibration residuals and input masks.
            block=start//512
            if len(f)>1 and block%16==int(np.random.default_rng(817993+block//16).integers(16)):
                fy,fe=sp.shuffle_fibres(f,e,722731+start)
                controls=scan_block(w,fy,fe,start,shuffle_stats,offsets=(0,))
                for rec in controls:
                    rec.update(stream="weak_heii_wavelength_shuffled",row_is_synthetic=True)
                    ledger.write(jsonline(rec)+"\n")
            iloc=[r-start for r in range(start,start+len(f)) if r in inrows]
            for i in iloc:
                d,truth=inject(w,f[i],e[i],rng)
                truth["row_index"]=i+start
                if d is not None:
                    fixed=scan_block(w,d[None,:],e[i:i+1],i+start,defaultdict(Counter),offsets=(0,))
                    found=any(abs(p["wave"]-truth["wave_heii"])<7 and
                              abs(p["pair_lya"]["wave"]-truth["wave_lya"])<7 for p in fixed)
                    preexisting=any(r["row_index"]==i+start and r["pair_offset_A"]==0 and
                                    abs(r["wave"]-truth["wave_heii"])<7 and
                                    abs(r["pair_lya"]["wave"]-truth["wave_lya"])<7 for r in recs)
                    truth.update(recovered=found,preexisting=preexisting,
                                 status="preexisting" if preexisting else "injected")
                trials.append(truth)
    known=known_eight(cat,observed)
    (out/"weak-heii-known-eight.json").write_text(json.dumps(json.loads(jsonline(known)),indent=2)+"\n")
    with gzip.open(out/"weak-heii-injections.jsonl.gz","wt") as ih:
        for r in trials:
            ih.write(jsonline(r)+"\n")
    bins={}
    for sn in (4.5,5.5,6.5,8.):
        subset=[r for r in trials if r["nominal_sn_heii"]==sn]
        unit=[r for r in subset if r["status"]=="injected"]
        bins[str(sn)]=dict(attempts=len(subset),in_band_with_valid_profiles=len(unit),
                           recovered=sum(r["recovered"] for r in unit),
                           preexisting=sum(r["status"]=="preexisting" for r in subset))
    controls={str(off):dict(stats[str(off)]) for off in OFFSETS}
    rates={str(off):stats[str(off)]["passed"]/max(stats[str(off)]["anchor_opportunities"],1) for off in OFFSETS}
    summary=dict(total=dict(stats["total"]),stages=controls,rate_per_anchor=rates,
        expected_displaced_pairs_at_true_exposure=stats["0"]["anchor_opportunities"]*
        (stats["-76.0"]["passed"]+stats["76.0"]["passed"])/
        max(stats["-76.0"]["anchor_opportunities"]+stats["76.0"]["anchor_opportunities"],1),
        shuffled_fibre_control=dict(total=dict(shuffle_stats["total"]),stages=dict(shuffle_stats["0"]),
            stratification="One seeded random 512-row block per 16 blocks; independently permute "
            "spectrum indices in each wavelength bin while moving the flux and error together."),
        observed_clean_pairs=len(observed),unique_sources=len({r["source_id"] for r in observed}),
        known_eight_with_rows=sum(bool(k["observations"]) for k in known),
        known_eight_blind_recovered=sum(any(x["weak_pair_blind_recovery"] for x in k["observations"]) for k in known),
        injection_seed=SEED,injection_bins=bins,
        domain="All catalog observations. z_HeII 1.9–2.4 truncated at spectral coverage. "
        "Lyα ≥7σ, HeII 4.5≤S/N<7 with observed flux ratio ≥0.1. "
        "Both <900 km/s FWHM, clean multi-pixel profiles. "
        "HeII EW0 ≥3 Å when continuum ≥3σ. Any ≥3σ C IV, C III], O III] component vetoes. "
        "Undetected or out-of-band metals still require deep external spectroscopy.",
        control="±76 vacuum-Angstrom cyclic displacement of the HeII prediction over identical "
        "anchor range. Same anchor, masks, fit, velocity −300 to +600 km/s, profile, ratio and metal cuts. "
        "Rates are spectral preselection chance-alignment diagnostics; offsets are not independently vetted.",
        injections="Random z uniform 1.9–2.4, source spectra including real noise, calibration masks and lines; "
        "Gaussian widths 2, 3, 4 Å, uniform velocity −300..600 km/s. Matched nominal HeII S/N 4.5, 5.5, "
        "6.5, 8 and Lyα S/N 9, 14, 22. Recovery requires both centroids within 7 Å, full spectral "
        "preselection, and absence of a pre-existing pair. Additional blank-sky fibres are separately "
        "measured by the principal search. No injection success implies astrophysical completeness.")
    (out/"weak-heii-summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return observed
