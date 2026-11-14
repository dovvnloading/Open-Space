# Independent harmonic measurement of OGLE-LMC-CEP-1429

The LMC east extension reported a roughly 1.196-day image signal near this known first-overtone Cepheid. We downloaded new MAST cutouts and the OGLE light curves, then fitted the catalog 3.599551-day pulsation and its first four Fourier orders simultaneously. The third order predicts **1.199850 days**. A 1.2-pixel aperture centered on the Gaia Cepheid gives the following semiamplitudes. [Prior provisional grid investigation](LMC_EAST_PRELIMINARY.md) · [OGLE Cepheid collection](https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/cep/) · [Measurements, exact inputs and aperture sensitivities](../results/lmc-east-cepheid-harmonics.json)

| Source | Accepted epochs | Fundamental amplitude | Third harmonic | Third / fundamental |
| --- | ---: | ---: | ---: | ---: |
| TESS Sector 1 | 1,267 | 14.49 e⁻/s | 0.895 e⁻/s | 6.17% |
| TESS Sector 2 | 1,196 | 15.84 e⁻/s | 0.986 e⁻/s | 6.23% |
| OGLE-III I | 792 | 0.0871 relative flux | 0.0045 relative flux | 5.16% |
| OGLE-IV I | 736 | 0.0904 relative flux | 0.0051 relative flux | 5.66% |

Four additional TESS visits (S3, S5, S8, S12) gave ratios of 8.99%, 6.92%, 10.53% and 6.81%. The range illustrates the sensitivity of this weak harmonic to observing geometry, contaminants and the photometric extraction. The OGLE-III and OGLE-IV V-band ratios were 5.16% and 5.44%. These measurements support the published Cepheid waveform as an explanation for the short-period signal. A distinct Gaia source 1.06″ from the Cepheid remains unresolved by TESS. The separate unpublished field-wide grid counts retain their provisional status. OGLE follow-up of this star has also reported other frequencies in its light curve. [Published Cepheid frequency analysis, online table A3](https://users.camk.edu.pl/smolec/cep1O_MNRAS/appendix_online.pdf)

The TESS fits use calibrated fluxes, a per-frame spatial nuisance fit outside four pixels and segment-specific quadratic time trends. OGLE magnitudes are converted to flux relative to each series' median. Four aperture radii and all fitted amplitudes and phases are in the evidence JSON. From a checkout with the archive dependencies installed, run the following command. [Analysis code](../examples/lmc_east_harmonic.py)

```sh
OPENBLAS_NUM_THREADS=1 python -m examples.lmc_east_harmonic runs/lmc-east-harmonic/cache runs/lmc-east-harmonic/result.json
```
