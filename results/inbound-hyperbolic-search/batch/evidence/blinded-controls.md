# Blinded interstellar controls

With the three names hidden, the published Find_Orb initial score failed to recover the accepted speed of every control in this experiment. On the same measurements, removing the initial hyperbola penalty recovered all three speeds. This test used a 120-second CPU limit per fit. The 18-second catalogue budget timed out on all nine control variants. [Long-budget fits, input hashes and retained arc lengths](blinded_controls.json) · [Short-budget outcomes](blinded_controls_18s.json)

| Object | Blinded, default initial score | Blinded, no initial hyperbola penalty |
| --- | --- | --- |
| 1I/‘Oumuamua | Bound solution using 9/180 positions on a 0.18-day arc | 26.34 km/s, 164/180 positions |
| 2I/Borisov | Bound solution using 19/180 positions | 32.28 km/s, 180/180 positions |
| 3I/ATLAS | 14.40 km/s, perihelion 0.00370 au, no eccentricity error reported | 57.97 km/s, perihelion 1.362 au, 160/180 positions |

We uniformly sampled 180 sorted optical records in each frozen MPC arc. All variants for an object used identical times, station codes and positions. The first twelve bytes carried either the known permanent designation or the same anonymous label. The unrestricted variant additionally prepended the native `#interstellar` directive. This directive removes a penalty on preliminary eccentricity above 1.01. Elliptical solutions remain available. The known-ID, anonymized and unpenalized headers and astrometry hashes are retained for replay. [Script](../blinded.py) · [Upstream scoring rule](https://github.com/Bill-Gray/find_orb/blob/1de55cca12e7a585c9c7384acbb80538ea2afc3c/orb_func.cpp#L3778)

For several anomalous solutions Find_Orb reported zero RMS and no eccentricity uncertainty. These outputs cannot establish a precise orbit. The longer-budget experiment demonstrates a failure mode of this automatic search. Three controls do not calibrate a catalogue-wide false-negative rate. The independent [event-based integrations](../../orbits/README.md) use a separate fitter and many more observations.
