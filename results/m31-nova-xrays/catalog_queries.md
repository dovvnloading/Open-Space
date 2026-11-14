# Rebuilding the catalog inputs

Get the optical table from [Shafter, Hornoch & Zhao](https://ashafter.sdsu.edu/m31_nova_table.csv). Its full downloaded SHA-256 was e154ab8e9dd72cba3556849fd1a51b263a575a2d8916a7a0ca70c1b5edb21616. Select rows whose JD is between 2456000 and 2459505.78 inclusive, and keep name, jd, ra_deg and dec_deg.

Use ADQL with the [Chandra CSC 2.2 TAP service](https://cda.cfa.harvard.edu/csc22tap). The following queries give the regional catalog sources, source-observation measurements and associations. For the sky control and pair search, join the first two October observations' wide-band measurements to the master position. The interval identifiers in the compound join key matter for a few older observations.

~~~sql
SELECT name, ra, dec, err_ellipse_r0, significance
FROM csc22.master_source
WHERE ra BETWEEN 9 AND 12 AND dec BETWEEN 40.1 AND 42.4

SELECT obsid, obi, region_id, gti_mjd_obs, instrument,
       flux_significance_w, flux_significance_b, flux_significance_u,
       flux_significance_s, flux_significance_m, flux_significance_h
FROM csc22.observation_source
WHERE ra_nom BETWEEN 9 AND 12 AND dec_nom BETWEEN 40.1 AND 42.4

SELECT ma.name, sa.obsid, sa.obi, sa.region_id
FROM csc22.master_stack_assoc AS ma
JOIN csc22.stack_observation_assoc AS sa
ON ma.detect_stack_id=sa.detect_stack_id AND ma.region_id=sa.region_id
WHERE sa.obsid IN (23871,26163)
~~~

Join associations to measurements by (obsid, obi, region_id), then to the master table by name. Keep every row with flux_significance_w >= 3 in either of the two visits. The replay takes the intersection of both visits. For longitudinal histories, repeat the association query with the last condition replaced by a restriction on ma.name to the seven CSC names in source_history.csv.gz. Keep one record per master source and observation interval. The additional history is for the weaker M31N 2019-07c sightline.

The [Chandra observatory TAP endpoint](https://cda.cfa.harvard.edu/cxctap) supplied the post-CSC observations. Deduplicate rows by obs_id before downloading.

~~~sql
SELECT obs_id,t_min,s_ra,s_dec,instrument_name FROM ivoa.ObsCore
WHERE s_ra BETWEEN 10.4 AND 10.95 AND s_dec BETWEEN 41.04 AND 41.5
AND t_min > 59100
~~~

The selected 2010 and 2011 HRC control observations sample the earlier history, including 11809, the strongest pre-eruption catalog measurement for the 2017-06c candidate. See observations.csv for the exact products. Chandra FITS DATE-OBS and the optical JD rounded to 0.1 day are used for the delays. Their last decimals specify arithmetic on the archived times, not an uncertainty on the explosion epoch.
