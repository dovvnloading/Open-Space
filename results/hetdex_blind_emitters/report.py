"""Produce the branch catalog entry and claim-first report from committed evidence."""
from collections import Counter
import csv
import gzip
import json
import math
import os
from pathlib import Path
import re

from calibration import summarize, write_table
from evidence_summary import summarize_evidence, adjudications, additional_ew

DIR=Path(__file__).resolve().parent
ROOT=DIR.parents[1]
EVIDENCE=DIR/'evidence'
RELEASES=('sc1','sc2')
CHANNELS=('heii','ew','blob','unknown')
NAMES={'heii':'He II','ew':'Extreme Lyα EW','blob':'Extended Lyα','unknown':'Unassigned-line screen'}
FIELDS=['release','source_id','detectid','shotid','ra_deg','dec_deg','channel','stream','wave_air_A',
        'secondary_ra_deg','secondary_dec_deg','secondary_wave_air_A','morphology',
        'trial_redshift','observed_EW_A','flux_cgs','flux_err_cgs','snr','rest_EW_A','spectral_rest_EW_lower_A',
        'imaging_rest_EW_lower_A','catalog_class','catalog_redshift','disposition',
        'independent_shot_max_snr','checks_passed','reasons','evidence_ledger','jsonl_line']


def read(path):
    return json.loads(path.read_text())


def number(x):
    return f'{x:,}'


def separation(a,b):
    x,y=map(math.radians,(a[1],b[1]))
    d=math.radians(a[0]-b[0])
    return math.degrees(2*math.asin(min(1,math.sqrt(math.sin((x-y)/2)**2+
                                            math.cos(x)*math.cos(y)*math.sin(d/2)**2))))*3600


def same_source(a,b):
    if a['channel']!=b['channel']:
        return False
    def endpoints(row):
        points=[(row['ra_deg'],row['dec_deg'],row['wave_air_A'])]
        if row['channel']=='blob' and row['secondary_ra_deg'] is not None:
            points.append((row['secondary_ra_deg'],row['secondary_dec_deg'],row['secondary_wave_air_A']))
        return points
    pa,pb=endpoints(a),endpoints(b)
    if len(pa)!=len(pb):
        return False
    def same(x,y):
        return (None not in x and None not in y and abs(x[2]-y[2])<5 and separation(x,y)<=1.25)
    return (all(same(x,y) for x,y in zip(pa,pb)) or
            (len(pa)==2 and all(same(x,y) for x,y in zip(pa,reversed(pb)))))


def compact(release,p,v,path,line):
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
                secondary_ra_deg=p.get('secondary_ra'),secondary_dec_deg=p.get('secondary_dec'),
                secondary_wave_air_A=p.get('secondary_wave'),morphology=p.get('morphology'),
                wave_air_A=p.get('wave'),trial_redshift=z,observed_EW_A=ew_obs,flux_cgs=p.get('flux'),flux_err_cgs=p.get('flux_err'),
                snr=p.get('snr'),rest_EW_A=ew_obs/(1+z) if z is not None and ew_obs is not None else None,
                spectral_rest_EW_lower_A=low_obs/(1+z) if z is not None and low_obs is not None else None,
                imaging_rest_EW_lower_A=(v.get('imaging_EW') or {}).get('ew0_lower_A'),
                catalog_class=p.get('source_type'),catalog_redshift=p.get('z_hetdex'),
                disposition=v['disposition'],independent_shot_max_snr=max(
                    (r['snr'] for r in v.get('repeats',[]) if r.get('snr') is not None),default=None),
                checks_passed=' | '.join(v.get('passed',[])),reasons=' | '.join(v.get('reasons',[])),
                evidence_ledger=f'{release}/{path}',jsonl_line=line)


def external_yield(directory):
    """Count eligible hypothesis rows and distinct successfully retrieved spectrum arrays."""
    with gzip.open(directory/'receipts.json.gz','rt') as file:
        receipts=json.load(file)
    successful_batches=Counter()
    batches=Counter()
    for row in receipts['catalog']:
        service=row['service']
        batches[service]+=1
        successful_batches[service]+=any(a.get('status')=='complete' for a in row['attempts'])
    arrays={r['array_sha256'] for r in receipts['spectra'].values()
            if r.get('status')=='measured' and r.get('array_sha256')}
    rows=Counter()
    with gzip.open(directory/'catalog.csv.gz','rt') as file:
        for r in csv.DictReader(file):
            rows['eligible']+=1
            rows['spectral_fit']+=int(r['measured_spectral_entries'])>0
            rows['agree']+=r['external_redshift_agrees'].lower()=='true'
            rows['prior_agree']+=r['prior_VAC_agrees'].lower()=='true'
    return batches,successful_batches,len(arrays),rows


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
    if os.environ.get("STAGE")=="vet" and not all((EVIDENCE/t/'focus'/'vet-summary.json').exists() for t in RELEASES):
        print("The first-pass vetting is checkpointed. The independent-spectroscopy and repeated-geometry jobs will publish the combined claim.")
        return
    if os.environ.get("STAGE")=="focus" and not all((EVIDENCE/t/'focus'/'saved-pair-summary.json').exists() for t in RELEASES):
        print("Both supplements are checkpointed. The finalization stage will reconcile individual visits and independent spectroscopy before cataloging the result.")
        return
    assert all((EVIDENCE/t/'focus'/'saved-pair-summary.json').exists() for t in RELEASES)
    allstats={}
    new=[]
    followed=[]
    counts=Counter()
    calib={}
    evidence_details={}
    unknown_atlases={}
    for tag in RELEASES:
        d=EVIDENCE/tag
        scan=read(d/'scan-summary.json')
        vet=read(d/'vet-summary.json')
        blanks=read(d/'offsource-summary.json')
        repeat=read(d/'repeat-summary.json')
        assert scan['observed']['spectra']==scan['input_rows'] and scan['shuffled']['spectra']>0
        assert scan['injections'].get('unblended',0)>0 and sum(vet['shortlist_counts'].values())>0
        calib[tag]=summarize(d)
        evidence_details[tag]=summarize_evidence(d)
        (d/'evidence-summary.json').write_text(json.dumps(evidence_details[tag],indent=2,allow_nan=False)+'\n')
        allstats[tag]=(scan,vet,blanks,repeat)
        tallied=Counter()
        raw_counts=Counter()
        effective_remote=Counter()
        unknown=Counter()
        for line,entry,ledger,initial in adjudications(d):
            p,v=entry['peak'],entry['vet']
            raw_counts[(p['channel'],initial)]+=1
            tallied[(p['channel'],v['disposition'])]+=1
            if 'NED' in v or 'native_ew_photometry' in v:
                followed.append(compact(tag,p,v,ledger,line))
                effective_remote[p['channel']]+=1
            if v['disposition']=='candidate':
                new.append(compact(tag,p,v,ledger,line))
            audit=v.get('unassigned_line_audit') or {}
            if audit:
                unknown['audited']+=1
                unknown['anchored']+=bool(audit.get('anchor_tests'))
                if audit.get('status')=='unmatched at independent redshift':
                    unknown['unmatched']+=1
                    atlases={a.get('repo') for a in audit.get('atlas_receipts',[]) if a.get('parsed_entries',0)>0}
                    unknown['unmatched_with_both']+= {'sdss/mangadap','linetools/linetools'}<=atlases
        unknown_atlases[tag]=unknown
        for ch,t in vet['status_counts'].items():
            for status,n in t.items():
                assert raw_counts[ch,status]==n
        for line,entry,path in additional_ew(d):
            p,v=entry['peak'],entry['vet']
            tallied[('ew',v['disposition'])]+=1
            effective_remote['ew']+=1
            row=compact(tag,p,v,path,line)
            followed.append(row)
            if v['disposition']=='candidate':
                new.append(row)
        vet['remote_assessed']=dict(effective_remote)
        counts.update(tallied)
    write_table(EVIDENCE/'candidates.csv',new,FIELDS)
    write_table(EVIDENCE/'followup.csv',followed,FIELDS)
    groups=[]
    for row in new:
        if not any(same_source(row,other) for other in groups):
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
      f'The follow-up ledger records external-check attempts for {len(followed):,} line or spatial hypotheses. Across all stages, {nunres:,} hypotheses remained unresolved. '
      'The releases share targets and exposures. Counts describe observations or feature hypotheses and can include repeated measurements of a galaxy.',
      '', '## Search coverage','',
      '| Release | Source observations scanned | Distinct release source IDs | Fitted peaks | Raw shuffled spectra | Off-source apertures | Same-position stacks from distinct shots |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for tag,(s,v,b,r) in allstats.items():
        sections.append(f"| {tag.upper()} | {s['input_rows']:,} | {s['distinct_release_source_ids']:,} | {s['observed'].get('peaks',0):,} | {s['shuffled']['spectra']:,} | {b.get('spectra',0):,} | {r['counts'].get('physical_groups',0):,} |")
    if new:
        sections+=['','### Sources passing the checks','',
                   '| Channel | Release / detection ID | RA, Dec (deg) | Trial z | Flux ± 1σ (cgs) | Spectral rest EW (Å) | Spectral EW lower bound (Å) | Imaging EW lower bound (Å) |',
                   '|---|---|---|---:|---:|---:|---:|---:|']
        for q in new:
            vals=[q[k] for k in ('rest_EW_A','spectral_rest_EW_lower_A','imaging_rest_EW_lower_A')]
            ew=' | '.join(f'{v:.1f}' if v is not None else '—' for v in vals)
            z=f"{q['trial_redshift']:.5f}" if q['trial_redshift'] is not None else 'unassigned'
            sections.append(f"| {NAMES[q['channel']]} | {q['release'].upper()} / {q['detectid']} | {q['ra_deg']:.6f}, {q['dec_deg']:.6f} | {z} | {q['flux_cgs']:.3g} ± {q['flux_err_cgs']:.2g} | {ew} |")
    sections+=['','The source search uses published PSF-weighted extractions and catalog emission knots. Diffuse sources absent from these inputs are outside the candidate search. Off-source cube apertures serve as controls. Each public 1D spectrum combines fibres and dithers from one shot. A map from that shot tests spatial morphology.',
              '', '## Leads and dispositions','',
              '| Channel | Single-visit screen, SC1 / SC2 | Excluded | Already published | Unresolved | Passed | External follow-up, SC1 / SC2 |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for ch in CHANNELS:
        ss=' / '.join(str(allstats[t][0]['observed'].get(ch,0)) if ch!='blob' else 'spatial' for t in RELEASES)
        rem=' / '.join(str(allstats[t][1]['remote_assessed'].get(ch,0)) for t in RELEASES)
        sections.append(f"| {NAMES[ch]} | {ss} | {counts[ch,'excluded']:,} | {counts[ch,'published']:,} | {counts[ch,'unresolved']:,} | {counts[ch,'candidate']:,} | {rem} |")
    sections+=['','Spatial proposals, subthreshold helium–hydrogen line pairs, imaging-led EW measurements and independent-visit stacks enter after the single-visit screen. A hypothesis can appear in more than one channel. Ranked follow-up gives priority to independent visits or DESI counterparts. Rows outside external follow-up retain an unresolved disposition.',
      '', 'An isolated line anywhere on this wavelength grid admits a Lyα redshift. The unassigned-line channel retains strong features displaced from common transitions at the catalog redshift. An unmatched line at an independent host redshift means no coincidence in the retrieved lists within the adopted velocity window. Wavelength coverage and completeness vary across those lists. A projected source can also contribute a line. The 5540 Å air red edge cuts off He II at about z = 2.378. Both C III] components fall within VIRUS only below about z = 1.903 in the requested interval. At higher redshift the metal test needs a redder spectrum. The evidence distinguishes an out-of-band component from a measured upper limit. A candidate galaxy spectrum cannot on its own determine a Population III or direct-collapse origin.',
      '', '| Unassigned-line follow-up | Audited hypotheses | With external redshift | Unmatched at that redshift in retrieved lists | Of unmatched, both external atlases returned transitions |',
      '|---|---:|---:|---:|---:|']
    for tag,u in unknown_atlases.items():
        sections.append(f"| {tag.upper()} | {u['audited']} | {u['anchored']} | {u['unmatched']} | {u['unmatched_with_both']} |")
    sections+=[
      '', '[All external follow-ups](evidence/followup.csv) · [Rows passing every gate](evidence/candidates.csv). Each release’s compressed vetting ledger contains IDs, comparisons, rejected alternatives, and pending checks. Fluxes are in erg s⁻¹ cm⁻², spectral densities per Å, wavelengths in observed air Å and EWs in rest-frame Å. Each trial redshift depends on the named identification.',
      '', '## Controls and recovery','',
      f'{nrec:,} of {ninj:,} successfully injected, unblended primary lines were recovered ({100*nrec/ninj:.1f}%) over the specified mix of amplitudes and widths. This is a spectral recovery fraction conditional on an extracted source. It does not include the probability of entering the public catalog, sky coverage, spatial completeness, or success in external follow-up.',
      '', '| Stream | Spectra | Fitted peaks | He II leads | Extreme EW leads | Off-list leads |','|---|---:|---:|---:|---:|---:|']
    for tag,(s,v,b,r) in allstats.items():
        for key,a in [('observed',s['observed']),('shuffled',s['shuffled']),('off-source',b)]:
            sections.append('| '+tag.upper()+' '+key+' | '+' | '.join(f'{a.get(k,0):,}' for k in ['spectra','peaks','heii','ew','unknown'])+' |')
    sections+=['','| Release | Matched real spectra / mock spectra | Real / mock leads per 1,000 spectra | Matched mock passes / leads | Raw-shuffle sample spectra | Raw mock passes / leads |',
               '|---|---:|---:|---:|---:|---:|']
    for tag in RELEASES:
        d=EVIDENCE/tag
        matched=read(d/'nuisance-summary.json')
        raw=read(d/'raw-shuffle-vet-sample.json')
        def vet_yield(kind):
            v=read(d/kind/'vet-summary.json')
            return sum(c.get('candidate',0) for c in v['status_counts'].values()),sum(v['shortlist_counts'].values())
        mk,mn=vet_yield('matched-shuffle-vet')
        rk,rn=vet_yield('shuffle-vet')
        sample=[matched[key] for key in ('observed','mock')]
        den=[arm.get('spectra',0) for arm in sample]
        rates=[]
        for n,arm in zip(den,sample):
            h=sum(arm.get(ch,0) for ch in ('heii','ew','unknown'))
            rates.append(f'{1000*h/n:,.1f}' if n else '—')
        sections.append(f"| {tag.upper()} | {den[0]:,} / {den[1]:,} | {' / '.join(rates)} | {mk:,} / {mn:,} | {raw['selected_spectra']:,} | {rk:,} / {rn:,} |")
    sections+=['','Each lead is a channel hypothesis. A peak can qualify for multiple channels. Passes count hypotheses that completed confirmation. External follow-up of control leads used the same ranked first-pass rules. The control rate is conditional on its construction and the matched population.',
               '', '| Published He II positive controls | Matched of eight | Single-line screen recovered | Weaker paired-line screen recovered | Recovered by both | Recovered by either |',
               '|---|---:|---:|---:|---:|---:|']
    additional=[]
    for tag in RELEASES:
        d=EVIDENCE/tag
        known=read(d/'positive-controls.json')
        weaker=read(d/'weak-heii-known-eight.json')
        weak=read(d/'weak-heii-summary.json')
        first={r['published_detectid'] for r in known['controls'] if r['recovered_heii_screen']}
        second={r['prior_detectid'] for r in weaker if any(o['weak_pair_blind_recovery'] for o in r['observations'])}
        matched=({r['published_detectid'] for r in known['controls'] if r['matched']} |
                 {r['prior_detectid'] for r in weaker if r['observations']})
        sections.append(f"| {tag.upper()} | {len(matched)} / {known['targets']} | {len(first)} | {len(second)} | {len(first & second)} | {len(first | second)} |")
        blank=read(d/'halo-injections-summary.json')
        ctr=allstats[tag][3]['counts']
        additional.append(f"{tag.upper()} produced {weak['observed_clean_pairs']:,} weaker spectral pairs against "
                          f"{weak['expected_displaced_pairs_at_true_exposure']:.1f} expected from displaced-wavelength controls before vetting. "
                          f"It also tested {ctr.get('control_groups',0):,} wrong-position, distinct-shot stacks and "
                          f"{blank['usable_position_wavelength_pairs']} usable blank position/wavelength pairs for halo injection.")
    sections+=['','Recovery counts use unique published objects within each release. Each screen used its own S/N and pair cuts. '+ ' '.join(additional),'',
               '| Release | Nominal injected S/N | Clean trials | Primary line recovered |',
               '|---|---:|---:|---:|']
    for t in RELEASES:
        for sn in (5,7,10,15,25):
            bins=[a for a in calib[t] if a['dimension']=='nominal_snr' and float(a['level'])==sn and a['selection']=='peak']
            n=sum(a['trials'] for a in bins)
            yes=sum(a['recovered'] for a in bins)
            if n:
                sections.append(f"| {t.upper()} | {sn} | {n} | {yes} ({100*yes/n:.1f}%) |")
    sections+=['','The raw wavelength-by-wavelength shuffle mixes bright and faint continua and destroys adjacent-pixel covariance. Its very large positive counts make it a stress test. After observing this excess, we added a same-shot block shuffle that preserves source continua. Its rates apply to the matched population. Off-source apertures preserve cube interpolation and spatial noise. Empty sky has no catalog redshift, so its off-list test is undefined. The per-release tables break recovery down by input flux, S/N, wavelength, profile width and input rest EW. Wilson intervals are nominal per-trial binomial intervals. Six trials share each injected background, and the intervals omit that clustering. The paired-visit controls substitute unrelated sources observed in the same shots.',
               '', '## Inputs, evidence and reproduction','']
    native_dirs=[(t,EVIDENCE/t/'ew-imaging') for t in RELEASES if (EVIDENCE/t/'ew-imaging'/'summary.json').exists()]
    if native_dirs:
        sections+=['','### Native Legacy Surveys equivalent widths','',
                   'Independent line identification and corroboration select the native-imaging cohort. Tractor fluxes and PSF-matched coadds provide redward continuum measurements. Limits require calibrated band coverage and a valid flux uncertainty. The two-band upper bound refers to a compact source with continuum-dominated r and z photometry, under the stated continuum and aperture model. It does not constrain diffuse emission beyond those apertures.',
                   '', '| Release / stream | EW rows in first register | Independently supported rows rechecked | Further supported rows outside the imaging scout | Rows with a lower bound | Lower bound >240 Å | Two-band model upper bound <240 Å | Confirmed |',
                   '|---|---:|---:|---:|---:|---:|---:|---:|']
        for tag,directory in native_dirs:
            summary=read(directory/'summary.json')
            for name,row in summary['counts'].items():
                vals=[row.get(k,0) for k in ('ew_rows_in_register','independently_corroborated_eligible',
                      'outside_scout_with_independent_support','with_EW_lower','lower_above_240',
                      'two_band_model_upper_below_240','disposition_candidate')]
                sections.append(f"| {tag.upper()} {name} | "+' | '.join(f'{x:,}' for x in vals)+' |')
        for tag,_ in native_dirs:
            base=f'evidence/{tag}/ew-imaging'
            sections.append(f'\n{tag.upper()} · [Every measured or inconclusive EW]({base}/outcomes.csv) · [Tightest bounds and confirmed rows]({base}/best.json) · [Native photometry and checksums]({base}/native.jsonl.gz).')
        sections.append('')
    sections+=['Every unresolved narrow He II / Lyα pair meeting the flux-ratio and spectral cuts, and every unresolved EW lead with a qualifying lower bound in the frozen first-pass register, entered the supplementary positional search of independent spectral catalogs. The table distinguishes successful catalog queries from attempts. Full rechecks preserve the original ledger and each revised decision.',
               '', '| Release | Eligible source / line hypotheses | Distinct positions | Rechecked hypotheses |',
               '|---|---:|---:|---:|']
    for t in RELEASES:
        f=read(EVIDENCE/t/'focus'/'summary.json')['counts']
        sections.append(f"| {t.upper()} | {f.get('eligible',0)} | {f.get('distinct_query_positions',0)} | {f.get('followup_records',0)} |")
    sections+=['','| Supplementary spectroscopy | SPARCL batches with valid response / attempted | CDS DESI batches with valid response / attempted | Distinct retrieved spectrum files by checksum | Eligible hypotheses with an external line fit | With an agreeing warning-free redshift | With an earlier agreeing DESI visual redshift |',
               '|---|---:|---:|---:|---:|---:|---:|']
    for t in RELEASES:
        batches,good,narrays,rows=external_yield(EVIDENCE/t/'focus')
        sections.append(f"| {t.upper()} | {good['SPARCL']} / {batches['SPARCL']} | {good['CDS_DESI_DR1']} / {batches['CDS_DESI_DR1']} | {narrays} | {rows['spectral_fit']} | {rows['agree']} | {rows['prior_agree']} |")
    sections.append('')
    sections.append('For supplementary redshifts, agreement requires |Δz| < 0.006. A single spectrum or position can contribute to several eligible hypotheses. Independent line fits and catalog redshift agreement are recorded separately.')
    sections+=['','| Native DESI DR1 files and single-cone catalog retry | Positions rechecked | Matched | No DESI DR1 counterpart | Partial or unavailable | Source hypotheses with new evidence re-evaluated |',
               '|---|---:|---:|---:|---:|---:|']
    for t in RELEASES:
        retry=read(EVIDENCE/t/'focus'/'static-summary.json')
        status=retry['position_status']
        sections.append(f"| {t.upper()} | {retry['searched_positions']} | {status.get('matched',0)} | {status.get('no_desi_dr1_counterpart',0)} | {status.get('partial_or_unavailable',0)} | {retry['new_adjudications']} |")
    sections+=['','The retry includes failed spectroscopic catalog queries and unresolved spectrum retrievals. The receipts distinguish files searched, sky pixels with no published product, empty target cones and failed requests. A completed DESI DR1 lookup constrains that survey only.',
               '', '[SC1 static catalog audit](evidence/sc1/focus/static-summary.json) · [SC2 static catalog audit](evidence/sc2/focus/static-summary.json).']
    sections.append('')
    for t in RELEASES:
        spatial=read(EVIDENCE/t/'repeat-blobs'/'summary.json')['counts']
        sections.append(f"{t.upper()} had {spatial.get('geometry_components',0):,} two-knot systems whose geometry and wavelengths recurred in distinct shots. Of these, {spatial.get('systems_requiring_two_shot_cube_check',0):,} met the local and prior-catalog selection for independent cube checks, {spatial.get('systems_with_measurable_bridge_in_distinct_shots',0):,} had measurable bridges in two shots, and {spatial.get('systems_with_both_knots_and_bridge_confirmed_twice',0):,} passed both knot-centered spatial tests twice.")
    sections+=['','The endpoint audits give the positions, trial identifications and dispositions of every system that passed both spatial tests twice. [SC1 endpoint measurements and source checks](evidence/sc1/blob-followup/systems.csv) · [SC2 endpoint measurements and source checks](evidence/sc2/blob-followup/systems.csv).']
    systems=[]
    for t in RELEASES:
        for item in read(EVIDENCE/t/'blob-followup'/'systems.json'):
            a,b=item['knots']
            reasons='; '.join(item.get('reasons') or item.get('passed') or [])
            systems.append(f"| {t.upper()} / {item['group']} | {a['ra']:.6f}, {a['dec']:.6f} | {b['ra']:.6f}, {b['dec']:.6f} | {a['wave_air_A']:.1f} | {a['trial_redshifts']['Ly_alpha']:.5f} | {item['disposition']} | {reasons.replace('|', '/')} |")
    if systems:
        sections+=['','| Release / system | First knot RA, Dec (deg) | Second knot RA, Dec (deg) | Observed line (Å) | Trial Lyα z | Disposition | Evidence or remaining ambiguity |','|---|---|---|---:|---:|---|---|',*systems]
    sections.append('')
    sections+=['| Check on ranked features | SC1 | SC2 |',
               '|---|---:|---:|']
    for key,name in [('NED','NED cone searches with a valid response'),('SIMBAD','SIMBAD cone searches with a valid response'),
                     ('Legacy','Legacy Surveys cone searches with a valid response'),('HSC','Direct HSC cutouts with measured photometry'),
                     ('cube','Cube opened and target narrowband window tested')]:
        values=[]
        for t in RELEASES:
            a=evidence_details[t]['archive_checks'][key]
            values.append(str(a.get('matched',0)+a.get('empty_cone',0)+a.get('measured',0)))
        sections.append('| '+name+' | '+' | '.join(values)+' |')
    for t in RELEASES:
        missing=allstats[t][1].get('catalog_errors') or {}
        if missing:
            labels={'LAN':'the published HETDEX Lyα-nebula catalog','DESI':'the HETDEX-specific DESI follow-up catalog'}
            detail=' and '.join(labels.get(k,k) for k in sorted(missing))
            sections+=['',f'{t.upper()} could not retrieve {detail} during the initial source and control vetting. That limited novelty checks and the imaging shortlist for this release. The independent spectral and repeated-geometry results above report their own coverage.']
    sections+=['','For each tested cube, the extension statistic subtracts the median from disjoint offset-wavelength windows and divides by 1.4826 times their median absolute deviation (MAD). The numerator counts source-row hypotheses with an excess ≥5 times this scaled MAD. The denominator counts source rows with a measurable extension statistic in their primary cube. Multiple source rows can refer to the same physical system.','',
               '| Spatial test on source rows | SC1 excess / measurable source rows | SC2 excess / measurable source rows |',
               '|---|---:|---:|']
    for key,name in [('two_knots','Emission between two knots'),('single_center','Diffuse emission around one center')]:
        values=[]
        for t in RELEASES:
            a=evidence_details[t]['spatial_arms'][key]['excess_rate_when_measurable']
            values.append(f"{a['numerator']} / {a['denominator']}")
        sections.append('| '+name+' | '+' | '.join(values)+' |')
    sections+=['','The release summaries retain each unresolved reason, paired-metal upper limits, and exact ledger rows for the strongest remaining leads.','',
               '[SC1 evidence breakdown](evidence/sc1/evidence-summary.json) · [SC2 evidence breakdown](evidence/sc2/evidence-summary.json).','']
    sections+=['| Release input | Bytes | SHA-256 |','|---|---:|---|']
    for t in RELEASES:
        for item in read(EVIDENCE/t/f'inputs-{t}.json'):
            sections.append(f"| {item['file']} | {item['bytes']:,} | `{item['sha256']}` |")
    for t in RELEASES:
        a=f'evidence/{t}'
        sections+=['',f'{t.upper()} · [Input checksums and mirrors]({a}/inputs-{t}.json) · [Ancillary and cube receipts]({a}/vet-inputs.json) · [Vetting ledger]({a}/vetted.jsonl.gz) · [Injections]({a}/injections.jsonl.gz) · [Completeness]({a}/completeness.csv) · [He II positive controls]({a}/positive-controls.json) · [Weaker paired lines]({a}/weak-heii-summary.json) · [Imaging-led EW measurements]({a}/ew-rescue-summary.json) · [Supplementary spectroscopy]({a}/focus/summary.json) · [Query receipts]({a}/focus/receipts.json.gz) · [Eligible spectroscopy ledger]({a}/focus/catalog.csv.gz) · [Repeated spatial systems]({a}/repeat-blobs/summary.json) · [Raw-shuffle follow-up]({a}/shuffle-vet/vet-summary.json) · [Continuum-preserving control]({a}/nuisance-summary.json) · [Matched control vetting]({a}/matched-shuffle-vet/vet-summary.json) · [Halo recovery]({a}/halo-injections-summary.json) · [Repeat controls]({a}/repeat-summary.json) · [Off-source spectra selection]({a}/offsource.jsonl.gz).']
    sections+=['','The workflow downloads both releases, verifies the publisher MD5 values and computes an independent SHA-256. A full run starts with empty release output directories and records its Python environment. The final job commits the evidence, builds this report and runs both repository checks. [Exact workflow commands](../../.github/workflows/hetdex-blind-emitters.yml).',
               '', 'Before the workflow is on the default branch, the following workflow-only commit launches a full reproduction on the research branch.',
               '', '```sh','git switch science/hetdex-blind-emitter-search','git pull --ff-only',
               'flow=.github/workflows/hetdex-blind-emitters.yml',
               'echo "# Reproduce the full HETDEX search." >> "$flow"','git add "$flow"',
               'git commit -m "[full-hetdex] Reproduce both public releases"','git push origin HEAD','```',
               '', 'Once GitHub lists the workflow on the default branch, its full dispatch is also available.',
               '', '```sh','gh workflow run hetdex-blind-emitters.yml --ref science/hetdex-blind-emitter-search -f stage=full','```',
               '', '[Spectral method](SPECTRAL_METHODS.md) · [Spatial, imaging and literature checks](VETTING.md) · [HETDEX source catalog 1 and public data](https://doi.org/10.5281/zenodo.7448504) · [HETDEX source catalog 2 and public data](https://doi.org/10.5281/zenodo.19581262) · [HETDEX data products and associated catalogs](https://hetdex.org/data-results/) · [Prior HETDEX He II search](https://arxiv.org/abs/2609.16134).',
               '', 'These measurements use the public source catalogs of Mentuch Cooper et al. (2023, ApJ 943, 177) and Mentuch Cooper et al. (2026, ApJS 284, 67), obtained with VIRUS on the Hobby–Eberly Telescope. The release records list the HETDEX collaboration, telescope partners, funding and requested acknowledgements.']
    (DIR/'README.md').write_text('\n'.join(sections)+'\n')
    if new:
        entry=f'**Rare ultraviolet emitters in HETDEX.** {len(groups)} source/channel associations passed the search criteria in {nobs:,} source observations.'
    else:
        entry=f'**A search for rare emitters in HETDEX.** Across {nobs:,} source observations, no new source passed all checks for unusually strong helium or hydrogen emission, extended hydrogen clouds or an unidentified line.'
    root_catalog(f'- {entry} *{label}* · [Details](results/{DIR.name}/README.md)', label)


if __name__=='__main__':
    main()
