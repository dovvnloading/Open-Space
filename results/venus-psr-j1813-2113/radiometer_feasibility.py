#!/usr/bin/env python3
"""Falsifiable 2026 hardware detectability scenarios for the 2027 Venus pulsar occultation.

Research ADDITION: independent of orbital computations. This is radiometer
sensitivity, NOT a measurement of J1813-2113 at 2.25 GHz. 2026 NASA IVS
Warkworth 2023/24 report pp94-97 supplies S-band receiver ranges and SEFD.
Earlier 2021/22 report mentions a separate L-band feed on the 12-m, whose
current operation and L-band sensitivity the 2026 report does not establish.

No file in PR #21 is modified by this executable supplement.
"""
import json
import math

# Keith et al. MeerKAT L-band data, **NOT** measured S-band flux:
L_BAND_MEDIAN_MJY=0.479
L_BAND_MAX_MJY=0.570
L_BAND_MIN_MJY=0.208

# Warkworth: NASA/TP-20260002317, p. 95 specification table and p. 97 update.
# These are S-band SEFDs, NOT L-band SEFDs. Receiver extent: 2.2-2.3 GHz.
S_BAND_30M_SEFD_JY=2100.0
S_BAND_12M_SEFD_JY=4800.0  # improved interim LNA from 2024 report
BANDWIDTH_HZ=100e6  # entire documented 2.2-2.3 GHz: optimistic
NPOL=2               # assumes both independent circular pol channels usable
INTEGRATION_S=360.0  # one six-minute continuous off-occultation baseline
DUTY=0.10            # unmeasured 10% effective duty; not a source assertion

def folded_snr(flux_mjy,sefd_jy,bandwidth_hz=BANDWIDTH_HZ,
               duration_s=INTEGRATION_S,n_pol=NPOL,duty=DUTY):
    assert 0<duty<1 and bandwidth_hz>0 and duration_s>0 and sefd_jy>0
    # Effective periodic radiometer equation, ideal calibration, full bandwidth.
    return flux_mjy*1e-3/sefd_jy * math.sqrt(n_pol*bandwidth_hz*duration_s) * math.sqrt((1-duty)/duty)

def flux_needed_mjy(nsigma,sefd_jy,duty=DUTY):
    return nsigma / folded_snr(1.,sefd_jy,duty=duty)

def sefd_needed_jy(nsigma,flux_mjy,duty=DUTY):
    return folded_snr(flux_mjy,1.,duty=duty)/nsigma

def main():
    med30=folded_snr(L_BAND_MEDIAN_MJY,S_BAND_30M_SEFD_JY)
    med12=folded_snr(L_BAND_MEDIAN_MJY,S_BAND_12M_SEFD_JY)
    extreme=folded_snr(L_BAND_MAX_MJY,S_BAND_30M_SEFD_JY,duty=.01)
    result={
      "provenance":"NASA-TP-20260002317 Table 1; Keith et al. L-band flux (not S-band measurement)",
      "observed_frequency":"S band 2200-2300 MHz; source flux extrapolated hypothetically flat from MeerKAT L band",
      "assumptions":{"bandwidth_mhz":100,"n_polarizations":2,"duration_seconds":360,
                     "baseline_effective_duty_fraction":.10,
                     "optimistic_source_flux_mjy_at_S":L_BAND_MEDIAN_MJY},
      "Warkworth_30m_S": {
        "documented_sefd_jy":S_BAND_30M_SEFD_JY,"ideal_snr_6min":round(med30,4),
        "flux_for_5sigma_mjy":round(flux_needed_mjy(5,S_BAND_30M_SEFD_JY),3),
        "sefd_for_5sigma_jy":round(sefd_needed_jy(5,L_BAND_MEDIAN_MJY),2)
      },
      "Warkworth_12m_S":{"documented_interim_sefd_jy":S_BAND_12M_SEFD_JY,
                         "ideal_snr_6min":round(med12,4)},
      "favorable_extreme_conditional":{
        "L_band_brightest_value_assumed_still_valid_at_S_mjy":L_BAND_MAX_MJY,
        "duty_fraction":.01,"30m_ideal_snr_6min":round(extreme,4),
        "max_ideal_snr_if_S_flux_le_L_band_max_and_duty_ge_1pct":round(extreme,4)
      },
      "scope":"An S-band hypothetical upper scenario, not evidence that no Warkworth L-band receiver exists."
    }
    # Frozen numerical contracts; these do not infer the actual source spectral index.
    assert .18<med30<.19 and .08<med12<.09
    assert extreme<.75
    assert 13<flux_needed_mjy(5,S_BAND_30M_SEFD_JY)<13.1
    assert sefd_needed_jy(5,L_BAND_MEDIAN_MJY)<80
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
