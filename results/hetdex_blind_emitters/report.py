"""Produce the branch catalog entry and claim-first report from committed evidence."""
from collections import Counter
import csv
import gzip
import json
import math
from pathlib import Path
import re

from calibration import summarize, write_table

DIR=Path(__file__).resolve().parent
ROOT=DIR.parents[1]
EVIDENCE=DIR/'evidence'
RELEASES=('sc1','sc2')
CHANNELS=('heii','ew','blob','unknown')
NAMES={'heii':'He II','ew':'Extreme Lyα EW','blob':'Extended Lyα','unknown':'Unassigned-line screen'}
FIELDS=['release','source_id','detectid','shotid','ra_deg','dec_deg','channel','stream','wave_air_A',
        'trial_redshift','observed_EW_A','flux_cgs','flux_err_cgs','snr','rest_EW_A','rest_EW_lower_A','catalog_class',
        'catalog_redshift','disposition','independent_shot_max_snr','checks_passed','reasons']


def read(path):
    return json.loads(path.read_text())


def number(x):
    return f'{x:,}'


def separation(a,b):
    x,y=map(math.radians,(a['dec_deg'],b['dec_deg']))
    d=math.radians(a['ra_deg']-b['ra_deg'])
    return math.degrees(2*math.asin(min(1,math.sqrt(math.sin((x-y)/2)**2+
                                            math.cos(x)*math.cos(y)*math.sin(d/2)**2))))*3600


def compact(release,p,v):
    channel=p['channel']
    z=p.get('z_heii') if channel=='heii' else p.get('z_lya') if channel!='unknown' else None
    ew_obs=(p['ew0']*(1+p['z_lya']) if p.get('ew0') is not None else
            p.get('spectral_EW0_A')*(1+p['z_lya']) if p.get('spectral_EW0_A') is not None else None)
    low_obs=(p['ew0_lower']*(1+p['z_lya']) if p.get('ew0_lower') is not None else
             p.get('spectral_EW0_lower_A')*(1+p['z_lya']) if p.get('spectral_EW0_lower_A') is not None else None)
    if channel=='heii' and p.get('heii_rest_ew_lower_A') is not None and z is not None:
        low_obs=p['heii_rest_ew_lower_A']*(1+z)
        ew_obs=p['heii_rest_ew_A']*(1+z) if p.get('heii_rest_ew_A') is not None else None
    return dict(release=release,source_id=p.get('source_id'),detectid=p.get('detectid'),shotid=p.get('shotid'),
                ra_deg=p.get('ra'),dec_deg=p.get('dec'),channel=channel,stream=p.get('stream','spatial'),
                wave_air_A=p.get('wave'),trial_redshift=z,observed_EW_A=ew_obs,flux_cgs=p.get('flux'),flux_err_cgs=p.get('flux_err'),
                snr=p.get('snr'),rest_EW_A=ew_obs/(1+z) if z is not None and ew_obs is not None else None,
                rest_EW_lower_A=low_obs/(1+z) if z is not None and low_obs is not None else None,
                catalog_class=p.get('source_type'),catalog_redshift=p.get('z_hetdex'),
                disposition=v['disposition'],independent_shot_max_snr=max(
                    (r['snr'] for r in v.get('repeats',[]) if r.get('snr') is not None),default=None),
                checks_passed=' | '.join(v.get('passed',[])),reasons=' | '.join(v.get('reasons',[])))


def root_catalog(line,label):
    p=ROOT/'README.md'
    lines=[s for s in p.read_text().splitlines() if f'](results/{DIR.name}/README.md)' not in s]
    at=lines.index('### Galaxies and the distant universe')+2
    lines.insert(at,line)
    counts=Counter()
    for s in lines:
        if s.lstrip().startswith('- '):
            m=re.findall(r'(?<!\*)\*(Measured|Candidate|Prediction|Correction|Null result)\*(?!\*)',s)
            if m:
                counts[m[-1]]+=1
    assert label in counts
    old=lines[8]
    for pat,val in [('(?<=has )\\d+(?= entries)',sum(counts.values())),
                    (r'\d+(?= new measurements)',counts['Measured']),
                    (r'\d+(?= candidate discoveries)',counts['Candidate']),
                    (r'\d+(?= forecasts of specific events)',counts['Prediction']),
                    (r'\d+(?= corrections to published catalogs)',counts['Correction']),
                    (r'\d+(?= carefully recorded null)',counts['Null result'])]:
        old,n=re.subn(pat,str(val),old)
        assert n==1,(pat,lines[8])
    lines[8]=old
    p.write_text('\n'.join(lines)+'\n')


def main():
    allstats={}
    new=[]
    followed=[]
    counts=Counter()
    calib={}
    for tag in RELEASES:
        d=EVIDENCE/tag
        scan=read(d/'scan-summary.json')
        vet=read(d/'vet-summary.json')
        blanks=read(d/'offsource-summary.json')
        repeat=read(d/'repeat-summary.json')
        assert scan['observed']['spectra']==scan['input_rows'] and scan['shuffled']['spectra']>0
        assert scan['injections'].get('unblended',0)>0 and sum(vet['shortlist_counts'].values())>0
        calib[tag]=summarize(d)
        allstats[tag]=(scan,vet,blanks,repeat)
        tallied=Counter()
        with gzip.open(d/'vetted.jsonl.gz','rt') as f:
            for l in f:
                entry=json.loads(l)
                p,v=entry['peak'],entry['vet']
                tallied[(p['channel'],v['disposition'])]+=1
                if 'NED' in v:
                    followed.append(compact(tag,p,v))
                if v['disposition']=='candidate':
                    new.append(compact(tag,p,v))
        for ch,t in vet['status_counts'].items():
            for status,n in t.items():
                assert tallied[ch,status]==n
        counts.update(tallied)
    write_table(EVIDENCE/'candidates.csv',new,FIELDS)
    write_table(EVIDENCE/'followup.csv',followed,FIELDS)
    groups=[]
    for row in new:
        if not any(row['channel']==other['channel'] and abs(row['wave_air_A']-other['wave_air_A'])<5
                   and separation(row,other)<=1.25 for other in groups):
            groups.append(row)
    nobs=sum(v[0]['observed']['spectra'] for v in allstats.values())
    npeaks=sum(v[0]['observed'].get('peaks',0) for v in allstats.values())
    ninj=sum(v[0]['injections']['unblended'] for v in allstats.values())
    nrec=sum(v[0]['injections']['unblended_recovered'] for v in allstats.values())
    nblank=sum(v[2].get('spectra',0) for v in allstats.values())
    nunres=sum(n for (k,c),n in counts.items() if c=='unresolved')
    label='Candidate' if new else 'Null result'
    claim=(f'**{len(groups)} candidate source/channel associations passed the stated tests**'
           if new else '**No new emitter passed all confirmation checks**')
    sections=['# Rare emission lines in the two public HETDEX source catalogs','',
      f'{claim} in {number(nobs)} source observations. The blind scan found {number(npeaks)} fitted peaks at signal-to-noise of at least 7. '
      f'The follow-up ledger records {len(followed):,} ranked line or spatial features with external checks and {nunres:,} unresolved entries across all stages. '
      'The releases share targets and exposures. Counts describe observations or feature hypotheses and can include repeated measurements of a galaxy.',
      '', '## Search coverage','',
      '| Release | Source observations scanned | Distinct release source IDs | Fitted peaks | Shuffled spectra | Off-source apertures | Independent-position stacks |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for tag,(s,v,b,r) in allstats.items():
        sections.append(f"| {tag.upper()} | {s['input_rows']:,} | {s['distinct_release_source_ids']:,} | {s['observed'].get('peaks',0):,} | {s['shuffled']['spectra']:,} | {b.get('spectra',0):,} | {r['counts'].get('physical_groups',0):,} |")
    if new:
        sections+=['','### Sources passing the checks','',
                   '| Channel | Release / detection ID | RA, Dec (deg) | Trial z | Flux ± 1σ (cgs) | Rest EW (Å) |',
                   '|---|---|---|---:|---:|---:|']
        for q in new:
            ee=q['rest_EW_A']
            ew=f"{ee:.1f}" if ee is not None else 'continuum unconstrained'
            z=f"{q['trial_redshift']:.5f}" if q['trial_redshift'] is not None else 'unassigned'
            sections.append(f"| {NAMES[q['channel']]} | {q['release'].upper()} / {q['detectid']} | {q['ra_deg']:.6f}, {q['dec_deg']:.6f} | {z} | {q['flux_cgs']:.3g} ± {q['flux_err_cgs']:.2g} | {ew} |")
    sections+=['','The input is the published PSF-weighted extraction at every released source observation. A diffuse source with no catalog detection can enter only through the sampled cube measurements. The public spectra combine fibres and dithers. A same-shot map checks spatial morphology but provides no independent exposure.',
              '', '## Leads and dispositions','',
              '| Channel | Single-visit screen, SC1 / SC2 | Excluded | Already published | Unresolved | Passed | External follow-up, SC1 / SC2 |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for ch in CHANNELS:
        ss=' / '.join(str(allstats[t][0]['observed'].get(ch,0)) if ch!='blob' else 'spatial' for t in RELEASES)
        rem=' / '.join(str(allstats[t][1]['remote_assessed'].get(ch,0)) for t in RELEASES)
        sections.append(f"| {NAMES[ch]} | {ss} | {counts[ch,'excluded']:,} | {counts[ch,'published']:,} | {counts[ch,'unresolved']:,} | {counts[ch,'candidate']:,} | {rem} |")
    sections+=['','Spatial proposals, subthreshold helium–hydrogen line pairs, imaging-led EW measurements and independent-visit stacks enter after the single-visit screen. A hypothesis can appear in more than one channel. Ranked follow-up gives priority to independent visits or DESI counterparts. Rows outside external follow-up retain an unresolved disposition.',
      '', 'An isolated line anywhere on this wavelength grid admits a Lyα redshift. The unassigned-line channel retains strong features displaced from common transitions at the catalog redshift. The search supplies no evidence for a new atomic transition. The 5540 Å air red edge cuts off He II at about z = 2.378. Both C III] components fall within VIRUS only below about z = 1.903 in the requested interval. At higher redshift the metal test needs a redder spectrum. The evidence distinguishes an out-of-band component from a measured upper limit. A candidate galaxy spectrum cannot on its own determine a Population III or direct-collapse origin.',
      '', '[All external follow-ups](evidence/followup.csv) · [Rows passing every gate](evidence/candidates.csv). Each release’s compressed vetting ledger contains IDs, comparisons, rejected alternatives, and pending checks. Fluxes are in erg s⁻¹ cm⁻², spectral densities per Å, wavelengths in observed air Å and EWs in rest-frame Å. Each trial redshift depends on the named identification.',
      '', '## Controls and recovery','',
      f'{nrec:,} of {ninj:,} successfully injected, unblended primary lines were recovered ({100*nrec/ninj:.1f}%) over the specified mix of amplitudes and widths. This is a spectral recovery fraction conditional on an extracted source. It does not include the probability of entering the public catalog, sky coverage, spatial completeness, or success in external follow-up.',
      '', '| Stream | Spectra | Fitted peaks | He II leads | Extreme EW leads | Off-list leads |','|---|---:|---:|---:|---:|---:|']
    for tag,(s,v,b,r) in allstats.items():
        for key,a in [('observed',s['observed']),('shuffled',s['shuffled']),('off-source',b)]:
            sections.append('| '+tag.upper()+' '+key+' | '+' | '.join(f'{a.get(k,0):,}' for k in ['spectra','peaks','heii','ew','unknown'])+' |')
    controls=[]
    for tag in RELEASES:
        d=EVIDENCE/tag
        known=read(d/'positive-controls.json')
        blank=read(d/'halo-injections-summary.json')
        null=read(d/'shuffle-vet'/'vet-summary.json')
        mc=read(d/'matched-shuffle-vet'/'vet-summary.json')
        weak=read(d/'weak-heii-summary.json')
        k=sum(c.get('candidate',0) for c in null['status_counts'].values())
        n=sum(null['shortlist_counts'].values())
        b=blank['usable_position_wavelength_pairs']
        ctr=allstats[tag][3]['counts']
        controls.append(f"{tag.upper()} recovered {known['target_heii_screen_recoveries']} of {known['position_or_id_matches']} matched positions in the eight-object published He II sample at the single-line threshold. The weaker paired-line search found {weak['observed_clean_pairs']:,} spectral leads against {weak['expected_displaced_pairs_at_true_exposure']:.1f} expected from displaced matches, before source vetting. The matched wrong-position comparison produced {ctr.get('control_groups',0):,} independent-visit control stacks. In the bounded raw-shuffle follow-up, {k} passed among {n:,} feature hypotheses. The continuum-preserving control yielded {sum(mc['shortlist_counts'].values()):,} feature hypotheses before source vetting. Halo injection tests had {b} usable blank position/wavelength pairs.")
    sections+=['',' '.join(controls),'',
               '| Release | Nominal injected S/N | Clean trials | Primary line recovered |',
               '|---|---:|---:|---:|']
    for t in RELEASES:
        for sn in (5,7,10,15,25):
            bins=[a for a in calib[t] if a['dimension']=='nominal_snr' and float(a['level'])==sn and a['selection']=='peak']
            n=sum(a['trials'] for a in bins)
            yes=sum(a['recovered'] for a in bins)
            if n:
                sections.append(f"| {t.upper()} | {sn} | {n} | {yes} ({100*yes/n:.1f}%) |")
    sections+=['','The raw wavelength-by-wavelength shuffle mixes bright and faint continua and destroys adjacent-pixel covariance. Its very large positive counts make it a stress test. A separate same-shot block shuffle preserves source continua and supplies the matched control. Off-source apertures preserve cube interpolation and spatial noise. Empty sky has no catalog redshift, so its off-list test is undefined. The per-release tables break recovery down by input flux, S/N, wavelength, profile width and input rest EW. Wilson intervals describe binomial trial uncertainty conditional on the injected backgrounds. The paired-visit controls substitute unrelated sources observed in the same shots.',
               '', '## Inputs, evidence and reproduction','']
    sections+=['| Release input | Bytes | SHA-256 |','|---|---:|---|']
    for t in RELEASES:
        for item in read(EVIDENCE/t/f'inputs-{t}.json'):
            sections.append(f"| {item['file']} | {item['bytes']:,} | `{item['sha256']}` |")
    for t in RELEASES:
        a=f'evidence/{t}'
        sections+=['',f'{t.upper()} · [Input checksums and mirrors]({a}/inputs-{t}.json) · [Ancillary and cube receipts]({a}/vet-inputs.json) · [Vetting ledger]({a}/vetted.jsonl.gz) · [Injections]({a}/injections.jsonl.gz) · [Completeness]({a}/completeness.csv) · [He II positive controls]({a}/positive-controls.json) · [Weaker paired lines]({a}/weak-heii-summary.json) · [Imaging-led EW measurements]({a}/ew-rescue-summary.json) · [Raw-shuffle follow-up]({a}/shuffle-vet/vet-summary.json) · [Continuum-preserving control]({a}/nuisance-summary.json) · [Matched control vetting]({a}/matched-shuffle-vet/vet-summary.json) · [Halo recovery]({a}/halo-injections-summary.json) · [Repeat controls]({a}/repeat-summary.json) · [Off-source spectra selection]({a}/offsource.jsonl.gz).']
    sections+=['','The workflow pins Python dependencies, fetches both releases, checks the publisher MD5 values and hashes the received bytes. Run its **full** dispatch on the research branch, or use the same commands on a runner with network access and sufficient disk. Run from the repository root.',
               '', '```sh','D=results/hetdex_blind_emitters','python -m pip install -r "$D/requirements.txt"','python "$D/check_pipeline.py"',
               'for r in sc1 sc2; do','  python "$D/scan.py" "$r" scan --output "$D/evidence/$r" --cache "/tmp/hetdex-$r" --workers 3',
               '  python "$D/scan.py" "$r" vet --output "$D/evidence/$r" --cache "/tmp/hetdex-$r"','done',
               'python "$D/report.py"','python -m open_space.verify','python -m unittest discover -s tests','```',
               '', '[Spectral method](SPECTRAL_METHODS.md) · [Spatial, imaging and literature checks](VETTING.md) · [HETDEX source catalog 1 and public data](https://doi.org/10.5281/zenodo.7448504) · [HETDEX source catalog 2 and public data](https://doi.org/10.5281/zenodo.19581262) · [HETDEX data products and associated catalogs](https://hetdex.org/data-results/) · [Prior HETDEX He II search](https://arxiv.org/abs/2609.16134).',
               '', 'These measurements use public HETDEX observations obtained with VIRUS on the Hobby–Eberly Telescope. The release records list the HETDEX collaboration, telescope partners, funding and requested acknowledgements.']
    (DIR/'README.md').write_text('\n'.join(sections)+'\n')
    if new:
        entry=f'**Rare ultraviolet emitters in HETDEX.** {len(groups)} source/channel associations passed the search criteria in {nobs:,} source observations.'
    else:
        entry=f'**A search for rare emitters in HETDEX.** Across {nobs:,} source observations, no new source passed all checks for unusually strong helium or hydrogen emission, extended hydrogen clouds or an unidentified line.'
    root_catalog(f'- {entry} *{label}* · [Details](results/{DIR.name}/README.md)', label)


if __name__=='__main__':
    main()
