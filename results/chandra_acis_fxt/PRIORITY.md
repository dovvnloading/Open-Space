# Two localized flares with independent X-ray checks

**Chandra Source Catalog 2.2 already flagged the 78-photon flare in ObsID 23324 as variable in that observation. The 37-photon flare in ObsID 32380 came from a field pointed at a ROSAT source for a stellar-activity study, and the source was detected again the next day.** Both are measured flares. Their distances and physical origins remain unresolved. [Measured triggers and exact inputs](evidence/priority/44985879443/leads.json) · [Independent Chandra visits](evidence/priority/44985879443/history/history.json)

| ObsID | ICRS RA, Dec (deg) | 90% radius | MJD UTC | Aperture t95−t05 | 60 s peak aperture flux, 0.5–7 keV |
|---|---|---:|---:|---:|---:|
| 23324 | 127.769030, +9.705384 | 1.43″ | 59569.55085 | 2,881 s | (1.48 ± 0.66) × 10⁻¹² erg cm⁻² s⁻¹ |
| 32380 | 243.548595, −14.359688 | 1.24″ | 61199.79665 | 162 s | (4.04 ± 1.12) × 10⁻¹² erg cm⁻² s⁻¹ |

The flux errors cover counting statistics in the selected peak interval. The reported durations span the 5th to 95th percentiles of aperture photon times within the selected flare window. These follow-up measurements were chosen after the blind search and have no independent discovery-trial interpretation. [Photometry for 23324](evidence/priority/44985879443/23324/photometry/flux.json) · [Photometry for 32380](evidence/priority/44985879443/32380/photometry/flux.json)

## ObsID 23324

The flare was about nine arcminutes from the targeted star LHS 2021. The 3,600 s window contained 78 photons in 75 frames and 78 detector pixels. The aperture had six photons during 10.09 ks before the flare and 19 in 6.38 ks afterward. Within half the extraction radius these counts were 0, 44 and 12 before, during and after. A separate X-ray source 31.4″ away had 168, 58 and 98 photons in the same intervals. Its count rate remained similar while the flare source brightened and decayed. The detector positions followed the telescope dither. [Image and detector traces](evidence/priority/44985879443/23324/photometry/diagnostics.png) · [Independent aperture counts](evidence/priority/44985879443/23324/photon_check.json)

CSC 2.2 identifies 2CXO J083104.6+094219 at a separation of 0.64″, with broad-band variability index 10 and 115.36 net source counts in ObsID 23324. Six other public ACIS pointings passed a broad pointing-coordinate search; none supplied a verified independent exposure at the flare aperture. They give no prior or subsequent count-rate constraint. [CSC measurements](evidence/priority/44985879443/history/target-00/csc.json) · [Exposure checks](evidence/priority/44985879443/history/history.json)

Gaia DR3 599875677054390784 lies 1.54″ away with G = 20.64 and a two-parameter astrometric solution. The associated Pan-STARRS stacked source has r−i = 1.24 ± 0.04 and i−z = 0.96 ± 0.02. Its faintness, limited morphology measurements and nearby bright star leave its class uncertain. A different Gaia source 5.00″ away has G = 14.24 and parallax 0.781 ± 0.020 mas. The AllWISE source 4.90″ from the flare position is closer to that bright star. The evidence admits a cool stellar counterpart and supplies no secure galaxy or redshift. [Optical and infrared catalog responses](evidence/priority/44985879443/optical/review.json)

The event input is `acisf23324N003_evt2.fits.gz`, SHA-256 `507a5aaf26b492ef2b5120242d9d33e281803667c8c1edc125fa905e6a09bbd0`. [URL, byte count and decompressed checksum](evidence/priority/44985879443/23324/input.json)

## ObsID 32380

The 180 s window contained 37 source photons and one photon in the twelve-times-larger annulus. The aperture contained 15 photons in 4.18 ks before the flare and 122 in 5.72 ks afterward. A further observation, ObsID 32390 on 2026-06-09, contained 28 aperture and 13 annular photons in 9.95 ks of supported exposure, giving a net mean aperture rate of 0.00271 counts s⁻¹. Persistent or recurrent X-ray emission is plausible. [Flare diagnostic](evidence/priority/44985879443/32380/photometry/diagnostics.png) · [Repeat-visit photometry](evidence/priority/44985879443/history/history.json)

The spacecraft targeted 2RXSJ161409.7-142143 under proposal 20201077, *Stellar activity with TESS and Chandra*. The pointing coordinates were within 5″ of the flare centroid. A bright red Gaia object at 1.89″ has G = 14.04 and BP−RP = 3.45, with a two-parameter solution. Gaia DR3 marks this source (4329327244825074688) as photometrically variable. Optical, infrared and X-ray source positions differ by a few arcseconds, making motion or blending relevant. The targeting, colors and subsequent detection favor stellar activity. A geometric match to the older ROSAT position alone is weak given its 35″ uncertainty. [Proposal and public pointing metadata](evidence/inventory/inventory.json) · [Counterpart measurements](evidence/priority/44985879443/optical/review.json)

The event input is `acisf32380N001_evt2.fits.gz`, SHA-256 `3c5f46567187114d7deb741d2a709bbc4e3566801a90a9cdb3f23a7c66a944cf`. [URL, byte count and decompressed checksum](evidence/priority/44985879443/32380/input.json)

The archived GCN search parsed 46,206 circulars. It returned 17 and 23 documents on nearby dates respectively. The reviewed shortlists supplied no positional counterpart. Position parsing covered 720 of 2,735 documents associated with Einstein Probe. Non-GCN reports and notices without usable coordinates remain outside this check. [GCN search and circular receipts](evidence/priority/44985879443/alerts/review.json) · [Actions run](https://github.com/dovvnloading/Open-Space/actions/runs/44998051147)
