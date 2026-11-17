# Two additional ALMA claims

**ALMA source U is inconsistent with a bound planet at 500–1,000 au.** We remeasured the publicly released 2016 December Band 6 and 2017 March Band 7 continuum images and found a 0.123″ centroid displacement in 85.8 days. Earth moved 0.929 au transverse to the first sightline over that interval. Keeping the object at the measured sky position would require over nine times the maximum displacement possible for a bound planet in this distance range. This is a direct test of the source identification. The astrophysical nature of U remains open.

**Gna also fails the bound 500–1,000 au test** if its two reported positions correspond to one body. Its authors withdrew the preprint in January 2016 pending more detections. Their original bound interpretation was a 220–880 km body at roughly 12–25 au. The alternative planet-sized interpretation required an unbound object. [Original preprint, Table 1 and §4](https://arxiv.org/pdf/1512.02650v1) · [Author's withdrawal](https://arxiv.org/abs/1512.02650v2)

## U measurements

Liseau et al. reported one secure 405 GHz detection in the final 2016 paper, 4.24 ± 0.49 mJy at JD 2457144.632219328, FK5 J2000 14h39m28.491s, −60°49′51.83″. The paper entertained a low-albedo extreme trans-Neptunian object with unknown distance and size. Its figure considered diameters up to 6,000 km. Its later astrometric reanalysis shifted the 2015 position by 0.92″ and measured 3.2 ± 0.5 mJy. The subsequent visit detections at nearly the same sky position are tabulated in the 2019 follow-up. Dates for our propagation use the MJD column of its Table 4 for the same field and band, and coordinates use Table 6. [2016 journal article, Fig. 5 and Table 5](https://www.aanda.org/articles/aa/pdf/2016/10/aa29135-16.pdf) · [2019 follow-up, §§5.1–5.2 and Tables 4 and 6](https://arxiv.org/pdf/1904.03043)

| Public ALMA continuum product | Date | Band | Fitted peak (mJy/beam) | Peak / local rms | Fixed-position peak / local rms |
|---|---|---:|---:|---:|---:|
| uid://A001/X11e/X38 | 2015-05-02 | 8 | 3.242 | 12.62 | 3.35 |
| uid://A001/X87d/X97 | 2016-12-30 | 6 | 0.706 | 5.25 | 5.00 |
| uid://A001/X87d/X9b | 2017-03-26 | 7 | 2.048 | 5.69 | 5.26 |
| uid://A001/X87d/X9f | 2017-03-26 | 8 | 2.228 | 3.46 | 1.40 |

The fixed point is the published 2017 Band 7 position. We fitted the header's elliptical restoring beam and a local sloping background in the primary-beam-corrected Stokes I continuum image, allowing up to 1.5″ per axis for the table's fitted peak. The noise estimate is the scaled median absolute deviation of pixels at the same primary-beam radius, excluding U. The ratios describe local contrast with correlated pixels. Band 8 in 2017 gives little independent support at the fixed position. Band 6 and Band 7 give similar positions when fitted separately. In the discovery map the fitted peak lies 0.67″ from the 2017 position. These are image-domain measurements of the released reductions; we did not independently calibrate the visibilities. This project used only public ALMA material for U. [ALMA Science Archive](https://almascience.eso.org/aq/) · [Measured positions, fluxes and beam sizes](alma_pixels.csv) · [Product IDs, exact FITS member names and SHA-256](receipts.json)

## Dynamics and 2023–2026 counterfactuals

For each pair, we projected the heliocentric Earth displacement perpendicular to the first measured sightline. At a geocentric range up to 1,002 au, the angular difference accounts for at most 1,002 sin(θ) au of this displacement. We added 2″ to the measured separation, allowed up to three extra days for changing light time and used the solar escape speed at 495 au as an upper speed throughout each short interval. With starting radius at least 500 au, a bound body cannot cross the 5 au buffer over either interval. Gna's dates have no time of day, so we allowed movement of either observation within the full UTC date using a combined one-day time allowance and a 0.018 au/day upper bound on Earth's speed.

| Input pair | Measured separation | Minimum required body displacement | Maximum bound displacement |
|---|---:|---:|---:|
| U, 2016-12-30 / 2017-03-26, published astrometry | 0.088″ | 0.9190 au | 0.0971 au |
| Gna, 2014-03-20 / 2014-04-14, published astrometry | 6.249″ | 0.0686 au | 0.0317 au |

These inequalities apply to all orbit headings and eccentricities across the stated distance interval, conditional on associating each pair of reported positions. They do not exclude objects elsewhere on the sky. For Gna we use the numerical Table 1 positions, which give a slightly larger displacement than the paper's approximate 5.9″ sentence. They are 19h15m23.124s, −07°02′48.765″ and 19h15m23.534s, −07°02′50.106″, with respective peaks of 2.8 ± 0.2 and 3.0 ± 0.1 mJy/beam. The authors give a third epoch upper limit below 1.05 mJy/beam on 2014-05-26. We have not independently remeasured the Gna images. [Gna original Table 1](https://arxiv.org/pdf/1512.02650v1) · [Numerical tests](result.json) · [Transcribed inputs](astrometry.csv)

The accompanying 144 monthly rows summarize 103,680 counterfactual U positions, beginning from the original 2015 discovery direction, at 500, 700 and 1,000 au and 720 evenly spaced circular headings per radius. Propagation uses a solar two-body circular solution, ERFA's geocentre, and iterated light time. Each bounding box encloses a ring at its stated distance. The later ALMA observations rule out these tracks as forecasts for U. Consequently no 2023–2026 image exclusion fraction is assigned. [Monthly ephemerides](counterfactual_monthly.csv)

## Replay

Run `python test_candidates.py` with NumPy and pyerfa for the dynamical test and monthly envelopes. Add `--full` to generate the four individual yearly compressed track tables. Run `python alma_pixels.py /path/to/cache` with Astropy and SciPy as well to retrieve the four calibrated ALMA products, verify archive and exact member hashes and repeat the four image fits. The product links are pinned in the receipts from the two project queries.
