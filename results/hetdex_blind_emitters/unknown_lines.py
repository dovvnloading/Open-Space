"""Audit a blind peak at independently measured redshift and alternative redshifts."""
from base64 import b64decode
import gzip
from hashlib import sha256
from io import BytesIO
import math
import re
import shlex
from threading import Lock

from astropy.table import Table
import numpy as np
import requests

import spectral as spec


LOCK = Lock()
ATLAS = None
# Vacuum Angstroms. Include close doublets as correlated measurements rather
# than treating two generally blended components as independent evidence.
FAMILIES = {
    "lya": (1215.6701, ("nv1","civ1","heii","oiii_uv2","ciii2")),
    "heii": (1640.42, ("lya","nv1","civ1","oiii_uv2","ciii2")),
    "civ1548": (1548.204, ("lya","civ2","heii","ciii2")),
    "ciii1909": (1908.734, ("lya","civ1","heii","oiii_uv2","ciii1")),
    "mgii2796": (2796.354, ("mgii2","nev3426","oii1","neiii")),
    "mgii2803": (2803.531, ("mgii1","nev3426","oii1","neiii")),
    "oii_blend": (3728.48, ("neiii","hdelta","hgamma","hbeta","oiii5007","halpha")),
    "nev3426": (3426.864, ("nev3346","oii1","neiii","hbeta","oiii5007")),
    "neiii3869": (3869.86, ("oii1","hdelta","hgamma","hbeta","oiii5007")),
    "hdelta": (4102.892, ("oii1","hgamma","hbeta","oiii5007","halpha")),
    "hgamma": (4341.684, ("oii1","hdelta","hbeta","oiii5007","halpha")),
    "hbeta": (4862.683, ("oii1","hgamma","oiii4959","oiii5007","halpha")),
    "oiii4959": (4960.295, ("oii1","hbeta","oiii5007","halpha")),
    "oiii5007": (5008.240, ("oii1","hbeta","oiii4959","halpha")),
    "halpha": (6564.61, ("hgamma","hbeta","oiii5007")),
}
REST = {**spec.UV,**spec.OPT,"mgii1":2796.354,"mgii2":2803.531,
        "nev3346":3346.783,"nev3426":3426.864}


def _fetch_tree(repo,path):
    url=f"https://api.github.com/repos/{repo}/contents/{path}"
    r=requests.get(url,timeout=(6,15),headers={"User-Agent":"Open-Space-astro"})
    r.raise_for_status()
    return {t["name"]:t for t in r.json() if t["type"]=="file"}


def _blob(item):
    if item.get("size",0)>16_000_000:
        raise ValueError("reference blob exceeds 16 MB")
    r=requests.get(item["git_url"],timeout=(6,22),headers={"User-Agent":"Open-Space-astro"})
    r.raise_for_status()
    b=b64decode(r.json()["content"])
    if len(b)!=item["size"]:
        raise ValueError("truncated reference atlas")
    return b,dict(file=item["name"],git_blob=item["sha"],url=item["git_url"],
                  bytes=len(b),sha256=sha256(b).hexdigest())


def _manga(raw):
    lines=[]
    for text in raw.decode().splitlines():
        if not re.match(r"^\s*DAPEML\s+[0-9]",text):
            continue
        row=shlex.split(text)
        if len(row)<5:
            continue
        name,reference=row[2],row[4].lower()
        try:
            w=float(row[3])
        except ValueError:
            continue
        if reference.startswith("air"):
            w=float(spec.air_to_vac(w))
        elif not reference.startswith("vac"):
            continue
        if 900<w<12000:
            lines.append((name,w,"SDSS MaNGA DAP emission"))
    return lines


def _linetools(raw,name):
    if name.endswith(".gz"):
        raw=gzip.decompress(raw)
    if ".fits" in name:
        tab=Table.read(BytesIO(raw),format="fits")
    else:
        txt=raw.decode("utf-8")
        if not txt.lstrip().startswith("# %ECSV"):
            # The distributed ASCII linelist can also be an astropy fixed-width table.
            tab=Table.read(txt.splitlines(),format="ascii")
        else:
            tab=Table.read(txt.splitlines(),format="ascii.ecsv")
    keys={k.lower():k for k in tab.colnames}
    if not {"wrest","name"}<=keys.keys():
        return []
    result=[]
    for label, rest in zip(tab[keys["name"]],tab[keys["wrest"]]):
        if isinstance(label,bytes):
            label=label.decode("ascii","replace")
        rest=float(rest)
        if 900<rest<12000 and np.isfinite(rest):
            result.append((str(label).strip(),rest,"linetools atomic/ISM; emission not inferred"))
    return result


def atlas():
    global ATLAS
    with LOCK:
        if ATLAS is not None:
            return ATLAS
        data=[]
        sources=[]
        for repo,path in (("sdss/mangadap","mangadap/data/emission_lines"),
                          ("linetools/linetools","linetools/data/lines")):
            try:
                files=_fetch_tree(repo,path)
                if repo.startswith("sdss"):
                    names=sorted((n for n in files if n.startswith("elp") and n.endswith(".par")),
                                 key=lambda n:(n!="elpmpl11.par",n))[:4]
                else:
                    names=sorted((n for n in files if "linelist" in n.lower() or
                                  n.lower().startswith(("gal_em","galaxy_","ism_"))),
                                 key=lambda n:(not n.startswith("linelist."),n))[:4]
                if not names:
                    raise ValueError("no reference files at listed path")
                for name in names:
                    try:
                        raw,meta=_blob(files[name])
                        rows=_manga(raw) if repo.startswith("sdss") else _linetools(raw,name)
                        sources.append(dict(meta,repo=repo,parsed_entries=len(rows)))
                        data.extend(rows)
                        if len(rows)>=20:
                            break
                    except Exception as ex:
                        sources.append(dict(repo=repo,file=name,error=str(ex)[:220]))
            except Exception as ex:
                sources.append(dict(repo=repo,error=str(ex)[:220]))
        for k,w in REST.items():
            data.append((k,w,"internal reference transitions"))
        ATLAS=(sorted(set(data),key=lambda t:t[1]),sources)
        return ATLAS


def unknown_spectrum(wave,flux,err,peak,*,density_unit=1.0):
    w,dx=spec._grid(wave)
    f,e,_=spec._inputs(w,flux,err,density_unit)
    f,e=f[0],e[0]*peak["noise_scale"]
    lam,vac,sig=peak["wave"],peak["wave_vac"],peak["sigma_A"]
    alternatives=[]
    for name,(rest,keys) in FAMILIES.items():
        z=vac/rest-1
        if z<-.001:
            continue
        measured={key:spec._force(w,f,e,dx,REST[key]*(1+z),sig) for key in keys}
        valid={k:m for k,m in measured.items() if m["status"]=="measured"}
        supporting=[k for k,m in valid.items() if m["snr"]>=5 and
                    abs(m["air_wave"]-lam)>max(2.355*sig,6)]
        contrary=[]
        if name in ("oiii4959","oiii5007"):
            key,ratio=("oiii5007",2.98) if name=="oiii4959" else ("oiii4959",1/2.98)
            if key in valid:
                expected=ratio*peak["flux"]
                variance=math.hypot(ratio*peak["flux_err"],valid[key]["flux_err"])
                difference=(valid[key]["flux"]-expected)/variance
                if abs(difference)>=3:
                    contrary.append(dict(line=key,ratio_residual_sigma=difference))
        alternatives.append(dict(identity=name,z=z,companions=measured,
                                 independent_wavelength_support=supporting,ratio_tensions=contrary))
    blue=spec._slice(w,lam-145,lam-max(20,4*sig))
    red=spec._slice(w,lam+max(20,4*sig),lam+145)
    cb=spec._weighted_cont(f[blue],e[blue],minimum=20)
    cr=spec._weighted_cont(f[red],e[red],minimum=20)
    profile={"blue_continuum":None if cb is None else [cb[0]*1e-17,cb[1]*1e-17],
             "red_continuum":None if cr is None else [cr[0]*1e-17,cr[1]*1e-17]}
    if cb and cr:
        profile["red_minus_blue_continuum_snr"]=(cr[0]-cb[0])/math.hypot(cb[1],cr[1])
    # Use equal-width mirrored bins around the fitted line. Interpolation
    # splits the centroid bin by its geometric overlap with each wing.
    side=spec._local_side(w,f,e,lam,sig)
    off=w-lam
    radius=2.5*sig
    red_weights=np.clip(np.minimum(off+dx/2,radius)-np.maximum(off-dx/2,0),0,dx)
    blue_weights=np.clip(np.minimum(off+dx/2,0)-np.maximum(off-dx/2,-radius),0,dx)
    eligible=(red_weights+blue_weights)>0
    if side and np.all(np.isfinite(e[eligible])):
        weight=(red_weights-blue_weights)[eligible]
        diff=float(weight@(f[eligible]-side[0]))
        se=math.hypot(float(np.linalg.norm(weight*e[eligible])),side[1]*weight.sum())
        if se>0:
            profile.update(red_minus_blue_wing_flux=diff*1e-17,wing_err=se*1e-17,wing_snr=diff/se)
    return dict(alternatives=alternatives,lya_profile=profile,
                profile_scope="Red-side asymmetry and a continuum break can favor Ly-alpha; neither is obligatory. "
                              "Wing uncertainty conditions on a centroid fitted to the same data.")


def reassess_unknown(rec,audit,cache=None):
    """A transition coincidence or an unsupported alternate is not an identification."""
    table,manifest=atlas()
    anchors=[]
    desi=audit.get("DESI") or {}
    if (desi.get("dex_flag")==1 and desi.get("target_sep_arcsec") is not None and
        float(desi["target_sep_arcsec"])<=1.25 and desi.get("vi_quality") is not None and
        float(desi["vi_quality"])>=3 and desi.get("vi_z") is not None):
        anchors.append(("DESI HETDEX visual VAC",float(desi["vi_z"])))
    for m in (audit.get("SPARCL") or {}).get("matches",[]):
        try:
            if float(m["redshift_warning"])==0:
                anchors.append(("SPARCL "+str(m.get("sparcl_id",m.get("specid",""))),float(m["redshift"])))
        except (KeyError,ValueError,TypeError):
            continue
    anchors=[(a,z) for a,z in anchors if np.isfinite(z) and z>=-.001]
    fits=[]
    vac=rec.get("wave_vac",float(spec.air_to_vac(rec["wave"])))
    for source,z in anchors:
        rest=vac/(1+z)
        tol=max(400.,spec.C*1.5*float(rec.get("sigma_A",2.4))/vac)
        hits=[{"name":name,"rest_vac_A":rw,"offset_km_s":spec.C*(rest/rw-1),"catalog":catalog}
              for name,rw,catalog in table if abs(spec.C*(rest/rw-1))<=tol]
        hits.sort(key=lambda r:abs(r["offset_km_s"]))
        fits.append(dict(source=source,z=z,rest_vac_A=rest,velocity_halfwidth_km_s=tol,
                         n_matches=len(hits),closest=hits[:15]))
    measurements=rec.get("unknown_spectrum") or {}
    support=[a["identity"] for a in measurements.get("alternatives",[]) if a["independent_wavelength_support"]
             and not a["ratio_tensions"]]
    external_size=sum(m.get("parsed_entries",0) for m in manifest)
    status=("no independent spectroscopic anchor" if not fits else
            "transition coincidence at independent redshift" if any(a["n_matches"] for a in fits) else
            "common alternative has a supporting companion" if support else
            "incomplete reference atlas" if external_size<20 else "unmatched at independent redshift")
    return dict(status=status,anchor_tests=fits,supported_alternatives=support,
                unique_reference_entries=len(table),atlas_receipts=manifest,
                empirical=measurements,
                scope="Coincident catalog wavelengths include weak/absorption transitions. "
                      "A projected single-line emitter cannot be eliminated by a redshift measured for its neighbor.")
