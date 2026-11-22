readme of the light curve data

format:
MJDcenter 2-20keV[ph/s/cm2] err 2-4keV err 4-10keV err 10-20keV err 

The MJDcenter is the center time of the day (1 day LC) and of the center 
time of the orbit (1 orbit LC).
The MJD column is equally separated by 24h for daily-data
6 hour for 6-hour-data or 1.5 hour for orbit-data.
The flux data is the average of the one or several scans during the period.
It is the time on ISS. It is neither geocenter, nor barycenter.

Note that duration of one MAXI scan is about 60s which comes in every 92 minutes.

If you want to know the exact scan time (center of the scan dwell time)
please download "Scan data" in v6L, v6m, and v7L pages.
Or please use MAXI ondemand (click on-demand tab in http://maxi.riken.jp )
Specify, for example, 0.02 (day) (i.e. <90 minutes) for the light curve bin.
Then you will get the center time of the scan in the MJD column.

In both cases the MJD is in ISS time. You will need a further 
correction to get the barycentric time.
