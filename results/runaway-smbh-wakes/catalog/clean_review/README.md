# Visual review of the clean extragalactic MER selection

None of the 198 catalogue-selected VIS entries shows a convincing straight stellar wake terminating at a distinct galaxy in the inspected images. Of the 198 entries, 102 show a long crossing streak, 35 show a direct stellar diffraction spike, 29 show a detached or segmented streak, 20 fall in a bright-star pattern of uncertain origin, and 11 show multiple long streaks. One further track crosses a resolved bright galaxy. In a 120 arcsec context image it extends well beyond both sides of that galaxy. Further 360 arcsec context traced that line to Gaia DR3 2160979470332604672, a G=6.75 star at 269 arcsec from the catalogue entry. The parent catalogue report includes the full stellar-spike image. A single image cannot determine the exposure origin of a detached track.

The 14 pages show the full set of 42 arcsec VIS and H image pairs in numeric object-ID order, with cyan ticks marking the catalogued position. The independent display stretch in each band describes morphology only. A solid black or gray H region can be unobserved or masked. In particular, the gallery does not establish a color for every source. Observations are per catalogue entry. Several entries trace parts of the same streak in nearby or overlapping cutouts.

`observations.csv` records the visual label, evidence, page and panel, and comparison with the separate Gaia axis alignment. A visually identified bright-star ray requires an apparent radial connection in the displayed field. An uncertain bright-star field records nearby diffraction structure with ambiguous parentage. Some other long tracks may have bright parents outside the displayed field. Apparent projection near a galaxy is recorded separately from the visual class. `h_primary_axis_measured` records whether the parent's specified H aperture and parallel-control measurement succeeded. It does not assert an H detection.

The parent selection table was filtered to `mer_det_quality_flag == 0` in the three deep fields. The Local Dark Nebula 1641 pointings were excluded. Input FITS URLs and SHA-256 values for the 198 pairs are in the parent catalogue measurement table. The per-pair thumbnail SHA-256s are in the observation table. The extra VIS/H context query has its own JSON receipt.

After reconstructing the paired thumbnails with the parent catalogue script, the pages can be regenerated with the following commands. Directory variables refer to the selected catalogue CSV, frozen tile CSV, and thumbnail cache.

```bash
python make_sheets.py --selected "$SELECTION" --tiles "$TILES" --stamps "$STAMPS" --out .
python wide_context.py --stamps "$STAMPS" --out .
```
