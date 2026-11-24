# Research-note preparation review

Checked on **23 November 2026**, alongside the scientific follow-up packets.
The two existing manuscripts remain **drafts for author review**. This pass
corrected concrete scientific and presentation issues; it does not certify
all references or scientific interpretations for submission.

## Changes made

The OGLE draft now distinguishes 5,641 parsed alerts, 5,519 fitted events and
125 ranked candidates; uses the reproduced baseline range I = 18.112–18.129;
and replaces the July astrometric heuristic with the conditional spring
forecast. Its abstract identifies the unmeasured mass and the sensitivity of
the quoted black-hole probability to geometric weighting. The existing
photometric forecast is unchanged.

The infrared draft now treats the radio increase as tentative and describes
the SPHEREx color through the comparison sample. It includes the interacting
supernova counterexample to the proposed coronal-line diagnostic, with a
primary citation. The nearby-target message also dates the 2024 residual
flux and the November catalog searches explicitly. Both manuscripts record
ChatGPT's role in this revision and leave final author verification pending.

## Checks against primary documentation

The current [RNAAS preparation guidance](https://journals.aas.org/research-note-preparation-guidelines/)
specifies a 1,500-word limit including references and captions, a required
abstract within 150 words, and one figure or one table. RNAAS is editorially
moderated rather than peer reviewed. The [AAS AI guidelines](https://journals.aas.org/author-llm-guidelines/)
require disclosure and author responsibility for the submitted content.
These policies were retrieved on this review date; the manuscript has not
been submitted.

| Check | Result and scope |
| --- | --- |
| Figure/table structure | OGLE has one figure; NEOWISE has one table. |
| Rough source-level word screen | Both drafts are below 1,500 words and both abstracts below 150 with a simple TeX-command-stripping count; the exact current counts and source hashes are in [the receipt](publication-checks.json). This is not `texcount` or an editorial count. |
| Six UAT keyword identifiers | All names and identifiers agree with the official UAT source file retrieved on this date. |
| Masterson et al. 2024 | NASA's Swift bibliography confirms ApJ 961, 211 and the listed first authors; DOI 10.3847/1538-4357/ad18bb. |
| Necker et al. 2025 | Publisher/CDS records confirm A&A 695, A228; expanded the manuscript entry to Necker, Graikou, Kowalski, et al. |
| Smith et al. 2009 | The author manuscript confirms the coronal-line/infrared SN 2005ip counterexample and ApJ 695, 1334–1350. |
| Author metadata and archive | ORCID, affiliation, email and Zenodo DOI placeholders remain deliberately unfilled. |

Primary reference records:

- [Masterson entry in NASA's 2024 Swift bibliography](https://swift.gsfc.nasa.gov/results/refereed/2024.html)
- [Necker catalog ReadMe and author list](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/695/A228?format=html&tex=true)
- [Smith et al. author manuscript](https://arxiv.org/abs/0809.5079)
- [Official UAT repository](https://github.com/astrothesaurus/UAT), source `UAT_list.json`, SHA-256 `d1e2731ad825236ee11b7587a9be978e377980d429452cf9050f36060a9c8289`.

The remaining manuscript references were not all independently checked in
this pass. The broader priority question also remains open: a bounded web
search is not an official TNS, ADS or ATel clearance.

## Rendering attempt and reproducible build

The official CTAN distribution supplied **AASTeX 7.0.1**, whose class is
`aastex701.cls`. Both drafts now name that class. The downloaded package
SHA-256 is `8b76872176f499607d63c477629ffd06b90f1d765b3910e1d2fa02d2306e0ab9`.
[CTAN package](https://ctan.org/pkg/aastex) · [AAS package guidance](https://journals.aas.org/aastex-package-for-manuscript-preparation/)

**Compilation and visual inspection could not be completed in this
environment.** The available pdfLaTeX executable lacks `pdflatex.fmt` and
the installed TeX support tree. Tectonic's required support-bundle index
returned HTTP 403 through both the official relay and its documented data
endpoint. The AAS 7.0.2 ZIP also returned HTTP 403; CTAN's 7.0.1 download
succeeded. No rendered paper PDF is claimed or delivered.

With a complete TeX installation and the official class available, compile
each source in its own directory and inspect every page:

```sh
cd papers/rnaas/ogle-2026-blg-0058
latexmk -pdf -interaction=nonstopmode -halt-on-error note.tex
cd ../neowise-nuclear-flares
latexmk -pdf -interaction=nonstopmode -halt-on-error note.tex
```

The OGLE figure is the previously deposited `figure1.pdf`; its photometric
predictions were not changed. The layout, references, table width, abstract
count and final word count still need checking on rendered output.

## Remaining publication decisions

The author must supply their actual metadata, confirm the analysis and
disclosure, complete reference/priority checks, and choose a permanent data
release. No identity or DOI is inferred or invented here. The existing
[paper index](../../papers/rnaas/README.md) records the planned OGLE courtesy
contact and follow-up sequence; the updated message files are drafts, with
recipients and personal statements still to be verified before sending.
