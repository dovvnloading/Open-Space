# Recovery of nearby known Y dwarfs

The primary CatWISE selection recovered **31 of 41** UltracoolSheet entries with infrared spectral types beginning with Y and adopted catalogue parallaxes ≥50 mas. Thirty passed the image-fit quality and W2 ≤16.5 cuts. Each object counts once, including objects with several CatWISE fragments.

| Distance from adopted catalogue parallax | Known controls | Primary selection | Bright tier |
|---|---:|---:|---:|
| ≤5 pc | 2 | 2 | 2 |
| >5 to 10 pc | 17 | 15 | 14 |
| >10 to 20 pc | 22 | 14 | 14 |

The sample uses the catalogue's reported parallaxes without deconvolution or a cut on parallax significance. Its membership reflects published discoveries and spectral classifications. It measures recovery on these known objects, not detection probability for the underlying population.

Each search downloads all CatWISE rows within 64.8″ of the predicted WISE-epoch position. To identify the control, it requires a corrected motion vector within the larger of 0.2″/yr and 30% of the published motion, plus a position within 2.5″. The initial reference propagated the UltracoolSheet J2000 coordinates and proper motions. That gave 25 primary and 24 bright recoveries. WISE 0855's propagated coordinate differed by 3.19″ from a prediction using measured census astrometry close to the WISE reference epoch. Both computations searched the same already-downloaded cone.

The epoch audit identifies a unique 20 pc census anchor for 39 of 41 controls using position, motion and parallax agreement. The two others retain their original coordinates. At 2.5″ the audit recovers three CatWISE fragments of WISE 0855 and four of WISE J163940.86−684744.6. Increasing the radius to 5″ or 10″ gives the same 31 unique primary recoveries. The CSV files expose the selection, alternate reference position, motion tolerance and per-object multiplicity. The dedicated WISE 0855 NEOWISE and SPHEREx tests and the independent W2 epoch linker provide additional measurements of this control.

Run these from the result directory after restoring the pilot inputs.

~~~bash
python controls/recover.py
python controls/check_epochs.py
~~~
