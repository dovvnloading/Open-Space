"""Combine independent image, catalog and targeted checks for all 42 groups."""
import pandas as pd
from campaign import HERE


def main():
    groups=pd.read_csv(HERE/'coincidences.csv')
    pix=pd.read_csv(HERE/'pixel_apertures.csv').merge(groups[['candidate','sign']],on='candidate')
    pix['passes']=pix.snr*pix.sign>=5
    twice=pix.groupby(['candidate','band']).passes.all().groupby('candidate').sum()
    hosts=pd.read_csv(HERE/'host_neighbors.csv')
    explicit={
        'a370_long-003':'known broad-absorption-line quasar SMM J02399-0136 L1; z_spec 2.798',
        'ceers_june-002':'early dither 3 clips the galaxy core at a detector edge in both bands',
        'cosmos_medium-007':'individual exposure ranges overlap in both bands; medians change by at most 1.3 percent',
        'cosmos_medium-009':'foreground stellar counterpart with long-baseline motion; individual exposures reproduce no catalog-scale fade',
        'cosmos_medium-010':'individual exposures overlap across epochs in both bands; catalog amplitude fails independent check',
        'cosmos_medium-012':'one usable early aperture per tested band; overlapping or small exposure changes fail to reproduce catalog amplitude',
        'cosmos_medium-014':'individual exposures agree across epochs; large catalog fade fails',
        'cosmos_medium-022':'exposure-confirmed brightening at galaxy center; AGN or projected transient unresolved'
    }
    rows=[]
    for g in groups.itertuples():
        n=int(twice.get(g.candidate,0))
        neighbor=hosts[(hosts.candidate==g.candidate)&(hosts.catalog=='cosmos2020')&(hosts.sep_arcsec<=.5)].sort_values('sep_arcsec')
        if g.candidate in explicit:reason=explicit[g.candidate]
        elif n<2:reason='does not reproduce in two filters at both aperture radii'
        elif not neighbor.empty and neighbor.iloc[0].lp_type==1:reason='stellar SED preferred in COSMOS2020; no independent extragalactic identification'
        else:reason='mosaic variability only; no independent transient classification'
        rows.append(dict(candidate=g.candidate,ra=g.ra,dec=g.dec,catalog_bands=g.bands,
                         bands_passing_both_apertures=n,disposition=reason))
    pd.DataFrame(rows).to_csv(HERE/'dispositions.csv',index=False)


if __name__=='__main__':main()
