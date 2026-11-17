"""Cross-match image-selected sky positions against the field exclusions and EDR sources.

Usage: python vet_positions.py positions.csv CACHE_DIR COMPARISON_DIR [OUTPUT_DIR]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

positions_file, cache, out = map(Path, sys.argv[1:4])
peaks = pd.read_csv(positions_file)
comparisons = pd.read_csv(out / 'field_comparisons.csv')
if len(sys.argv) > 4:
    out = Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
stone = fits.getdata(cache / 'nexus_wide01_deep01_stacked_sources_allband.fits.gz')
edr = fits.getdata(cache / 'wide_ep1_01_src_catalog_v1.fits')
spec_targets = []
for epoch in range(1, 9):
    data = fits.getdata(cache / f'nexus_msa_deep_ep{epoch:02}_targets.fits')
    for obj in data[(data['z_VI_conf'] == 2) & (data['z_VI'] > 0)]:
        spec_targets.append((obj['RA'], obj['DEC'], int(obj['MSAID']), obj['z_VI'], epoch))
spec_targets = pd.DataFrame(spec_targets, columns=['ra','dec','msaid','z','epoch'])


def coords(ra, dec):
    return SkyCoord(np.asarray(ra)*u.deg, np.asarray(dec)*u.deg)


def prop(x):
    return round(float(x), 5) if np.isfinite(x) and x > 0 else None


pc = coords(peaks.ra_deg, peaks.dec_deg)
ec = coords(edr['RA'], edr['DEC'])
sc = coords(stone['RA'], stone['DEC'])
qc = coords(spec_targets.ra, spec_targets.dec)
cc = coords(comparisons.ra_deg, comparisons.dec_deg)
rows = []
links = []
for pos, p in zip(peaks.itertuples(), pc):
    ds = p.separation(cc).arcsec
    for j in np.flatnonzero(ds <= 2):
        obj = comparisons.iloc[j]
        links.append((pos.id, obj.catalog, obj.designation, round(ds[j],4), obj.status))
    def nearest_class(cats):
        ids = np.flatnonzero(comparisons.catalog.isin(cats))
        j = ids[np.argmin(ds[ids])]
        return comparisons.iloc[j].designation, np.round(ds[j],4)
    designation, known_sep = nearest_class(['Zhuang2026_Table2','TNS_TransientsFYI','TNS_Fink'])
    stone_name, vs = nearest_class(['Stone2026'])
    msa_name, ms = nearest_class(['NEXUS_MSA_priority0'])
    ed_sep = p.separation(ec).arcsec
    ei = np.argsort(ed_sep)[:2]
    e = edr[ei[0]]
    si = np.argmin(p.separation(sc).arcsec)
    # Same EDR identifier is mandatory when transferring redshift to an EDR host.
    # Use the positional check as a second guard against inconsistent identifiers.
    linked = np.flatnonzero(stone['NID'] == e['ID'])
    s = None
    if linked.size:
        j = linked[np.argmin(ec[ei[0]].separation(sc[linked]).arcsec)]
        if ec[ei[0]].separation(sc[j]).arcsec < .2:
            s = stone[j]
    qdist = ec[ei[0]].separation(qc).arcsec
    qi = np.flatnonzero(qdist < .15)
    zs = sorted(set(round(float(n),3) for n in spec_targets.iloc[qi].z))
    # A conflicting solution is retained in the audit cell with all robust values.
    ztext = ';'.join(str(z) for z in zs)
    sep = float(ed_sep[ei[0]])
    relation = ('centroid_coincident' if sep <= .15 else
                'nearby_offset' if sep <= 1 else 'distant_nearest')
    rows.append(dict(id=pos.id, ra_deg=pos.ra_deg, dec_deg=pos.dec_deg,
                     nearest_report=designation, report_sep_arcsec=known_sep,
                     nearest_published_variable=stone_name, variable_sep_arcsec=vs,
                     nearest_msa_priority0=msa_name, msa_priority0_sep_arcsec=ms,
                     nearest_edr_id=int(e['ID']), edr_sep_arcsec=round(sep,4),
                     second_edr_id=int(edr[ei[1]]['ID']), second_edr_sep_arcsec=round(ed_sep[ei[1]],4),
                     centroid_relation=relation, edr_kron_radius_arcsec=round(float(e['R_KRON']),4),
                     sep_over_edr_kron=round(sep/float(e['R_KRON']),3) if e['R_KRON'] > 0 else None,
                     stone_id=int(s['ID']) if s is not None else None,
                     edr_zphot=prop(s['Z_PHOT']) if s is not None else None,
                     edr_zspec_stone=prop(s['Z_SPEC']) if s is not None else None,
                     edr_zspec_qdr_conf2=ztext, qdr_nearest_to_edr_arcsec=round(float(qdist.min()),4),
                     nearest_stacked_id=int(stone[si]['ID']), stacked_sep_arcsec=round(p.separation(sc[si]).arcsec,4),
                     stacked_zphot=prop(stone[si]['Z_PHOT'])))
pd.DataFrame(rows).to_csv(out/'position_context.csv', index=False)
pd.DataFrame(links, columns=['id','catalog','designation','separation_arcsec','status']).sort_values(['id','separation_arcsec']).to_csv(out/'matches_within_2arcsec.csv', index=False)
