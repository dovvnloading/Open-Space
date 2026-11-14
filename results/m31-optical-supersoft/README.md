# A candidate optical eruption for XMMU J004347.8+413310 in Andromeda

Two optical observers recorded a 16.7-mag transient at the position of a later supersoft X-ray source on 2014 May 6–7. The Kitt Peak Hα survey independently measured a transient 91 days later, catalogued as M31N 2014-08a. Four XMM-Newton visits then detected a soft X-ray source between 270 and 625 days after May 6. After registering the two deepest XMM fields against Chandra positions, their centroids land 0.44″ and 0.88″ from the May 7 measurement. These observations make the 2014 optical transient a candidate counterpart to XMMU J004347.8+413310, catalogued as M31N 2015-06c. They also extend the observed X-ray timeline 147 days earlier than the June pointing cited in its discovery report.

## Who measured what

Nishiyama and Kabashima reported the May transient as PNV J00434777+4133129 on seven unfiltered exposures. They quoted a 19.2-mag image limit and a prior March image limit of 18.9 mag. T. Yusa confirmed it the following night with a different telescope. Their positions differ by 0.92″. Rector et al. later listed the August Hα detection among the newly reported novae in their survey. The open Hα photometry contains one detection bracketed by faint limits. The discovery observations used different filters and the May–August gap leaves open whether they sampled one sustained event. The May identifier matters here because the current catalogue's M31N 2014-05a is a different object. [Original optical reports](http://www.cbat.eps.harvard.edu/unconf/followups/J00434777+4133129.html) · [Rector et al. (2022)](https://doi.org/10.3847/1538-4357/ac87ad) · [Living Catalogue of Novae in M31](https://ashafter.sdsu.edu/m31_nova_table.html)

Henze et al. reported the 2015 June 28 discovery and fitted a 34 ± 1 eV blackbody. They measured 30 ± 2 eV in 2016 January, when they still lacked an optical counterpart and called for an archive search. The February 2015 detection and the independent 2015 June 27 pointing below were identified by cross-matching public EPIC source lists. The original analyses establish the supersoft spectrum. This investigation measures positions and epochs using released detection products. [Discovery, ATel 8227](https://www.astronomerstelegram.org/?read=8227) · [Follow-up, ATel 8826](https://www.astronomerstelegram.org/?read=8826) · [Northern disk X-ray survey](https://arxiv.org/abs/1809.08020)

| UTC date | Measurement | Days since May 6 optical detection |
| --- | --- | ---: |
| 2014 May 6 | 16.7 unfiltered, Nishiyama & Kabashima | 0.0 |
| 2014 May 7 | 16.7 unfiltered, Yusa | 0.7 |
| 2014 Aug 5 | 17.51 ± 0.07 Hα, Kitt Peak | 91.2 |
| 2015 Feb 1 | XMM 0744350901, source 53, HR1 = −0.998 | 270.2 |
| 2015 Jun 27 | XMM 0761970101, source 79, HR1 = −0.805 | 416.3 |
| 2015 Jun 28 | XMM 0763120101, source 3, HR1 = −0.941 | 417.5 |
| 2015 Aug 30 | Fainter than 20.72 in Hα, Kitt Peak | 481.1 |
| 2016 Jan 21 | XMM 0763120201, source 3, HR1 = −0.903 | 624.8 |

HR1 compares counts in 0.5–1.0 and 0.2–0.5 keV, using their difference divided by their sum. All four source records are point-like and clear of EPIC flags. Their detection likelihoods are 131, 170, 11,546 and 5,127, respectively. The first and last detections span 355 days. The sampling leaves the turn-on, turn-off and continuity between visits unresolved.

## Astrometry and competing sources

The raw June 28 position lies 3.07″ from the August Hα position. The original 2015 telegram reported a 1.5″ three-sigma radius at its own centroid. To test whether frame offsets matter, we registered each XMM field to the Chandra Source Catalog with a two-coordinate translation. We excluded the central 12″ before fitting. The primary 12′ fit used bright, unresolved XMM sources with detection likelihood ≥25, catalogued statistical position error ≤2″, and no second XMM detection within 9″. Each reference had a Chandra match within 5″ and no second Chandra source within 8″. We clipped separations beyond 2.75″ of the median translation in three passes. The closest Chandra source is 7.3″ from the August optical position, so Chandra contributes field references only.

| Deep XMM pointing | References | East, north shift | Scatter per coordinate | Separation from May 6 / May 7 / August positions |
| --- | ---: | --- | --- | --- |
| 2015 Jun 28 | 51 | −0.10″, +2.13″ | 0.77″, 0.78″ | 0.64″ / 0.44″ / 1.27″ |
| 2016 Jan 21 | 45 | −0.90″, −2.33″ | 0.67″, 0.63″ | 1.03″ / 0.88″ / 1.58″ |

Sasaki et al. (2018) independently registered these XMM data using USNO-B1.0 and 2MASS. Their table A.1 lists source 53 as a supersoft source, with position 00h43m47.788s, +41°33′12.38″ and quoted positional error 0.51″. That centroid lies 0.56″ from the May 6 detection, 0.77″ from the May 7 detection and 1.12″ from the August Hα position. This cross-check uses the same XMM observations with a different astrometric reference. [Sasaki et al., section 2.2 and table A.1](https://arxiv.org/pdf/1809.08020)

The deep centroids agree to 0.54″ after registration. Using an 8′ or 18′ aperture in place of 12′ moves either deep centroid by less than 0.4″. The scatter in the table measures dispersion among the individual reference sources. Optical astrometric errors also matter. Rector et al. give a global WCS accuracy better than 0.6″. The May observers do not publish a positional covariance.

A persistent Pan-STARRS source of mean r = 21.16 mag (42 detections, object 157860109495944763) lies 0.94″ and 1.16″ from the two registered X-ray centroids. The May discovery report already identifies a faint LGGS object in the field. Either a quiescent counterpart or an unrelated neighboring star could account for that source. The nearest separate Chandra source, 2CXO J004347.1+413317, is also close enough to enter an XMM extraction aperture. Blending in the XMM PSF can affect a position or a hardness ratio. In the optical catalogue there are zero other nova entries within 30″ of the August position and 21 within 5′, excluding the optical and X-ray designations under study. Those counts describe the local environment. Selection across the whole galaxy precludes using them as a discovery p-value. The data support a candidate association and a pre-discovery X-ray detection. The two catalogue designations are insufficient evidence of two separate eruptions.

## Selection and reproduction

The search compared 1,435 eruption catalogue positions with 17,493 XMM-Newton pipeline detections from 146 observations in a rectangle covering M31. Post-2011 optical entries were screened for a point-like source within 4″ observed between 0 and 1,200 days later, with EPIC likelihood >12, HR1 <0.15 and HR2 <0. A wider 7″ radius was used to inspect competing detections. The search recovered published counterparts including M31N 2012-06a and the annual recurrent M31N 2008-12a. The target here stood out because the X-ray observers explicitly sought an optical counterpart. The published Hα source, the two May optical positions, and the current catalogue's separate X-ray entry can each be traced to specific original records. A search of those records and identifier-based literature queries did not locate a prior publication making this association. Publication priority remains open. [XMM Science Archive](https://nxsa.esac.esa.int/) · [Chandra Source Catalog](https://cxc.harvard.edu/csc/) · [MPE optical catalogue](https://www.mpe.mpg.de/~m31novae/opt/m31/M31_table.html)

Run from this directory with Python 3.10 or later. The replay uses only the standard library, 625 frozen XMM detections from the four visits, 417 Chandra positions, and the frozen optical and nearby-source rows. It calculates the offsets, repeat-visit separations, matched reference counts and local nova counts.

```sh
python replay.py > replay-output.json
diff -u measurements.json replay-output.json
```

`optical.csv` preserves the individual passbands and faint limits. `xmm-field.csv` includes the statistical position errors, likelihoods, hardness uncertainties, estimated counts and pipeline flags. `sources.json` records the archive queries, input sizes and original response hashes. The hashes describe the retrieved snapshots and live services may update. The nova catalogue is by A. W. Shafter, K. Hornoch and J. Zhao (accessed 2026 November 14, CC BY 4.0), incorporating pre-May-2022 data from W. Pietsch's MPE catalogue. CBAT observers, Rector et al., Henze et al., the XMM-Newton project, the Chandra X-ray Center and the Pan-STARRS1 surveys provided the underlying observations and products.
