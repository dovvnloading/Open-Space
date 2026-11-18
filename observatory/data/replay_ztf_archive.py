"""Independently recover J0826+7002 nightly photometry from public ZTF DR24 data.

Run `python observatory/data/replay_ztf_archive.py` with numpy, pandas, astropy
and pyarrow (including S3 support). Requires public HTTPS and anonymous S3
access. Reads one sky partition and the IRSA science-exposure metadata. The
routine observatory build does not run this network audit.
"""

import hashlib
import io
from decimal import Decimal
from pathlib import Path
from urllib.request import urlopen

import astropy.units as u
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.fs as fs
import pyarrow.parquet as pq
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers


RESULTS = Path(__file__).resolve().parents[2] / "results/j0826-lens-delay"
TARGET = "J0826+7002"
PARTITION = (
    "ipac-irsa-ztf/ztf/enhanced/dr24/lc/hats/ztf_dr24_lc-hats/dataset/"
    "Norder=5/Dir=0/Npix=1949/part0.snappy.parquet"
)
EXPOSURES = (
    "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?"
    "POS=126.5333358%2C70.0464092&COLUMNS=field%2Cccdid%2Cqid%2Cfiltercode"
    "%2Cobsjd%2Cexptime%2Cmaglimit%2Cfilefracday&RESPONSEFORMAT=CSV"
)
CENTER = SkyCoord(126.5333358 * u.deg, 70.0464092 * u.deg)
RADIUS = 0.0020917 * u.deg
BANDS = {1: "zg", 2: "zr", 3: "zi"}


def cone_lightcurves():
    """Read position columns first so only intersecting row groups need photometry."""
    archive = fs.S3FileSystem(region="us-east-1", anonymous=True,
                              request_timeout=90, connect_timeout=15)
    partition = pq.ParquetFile(PARTITION, filesystem=archive)
    for group in range(partition.num_row_groups):
        coords = partition.read_row_group(group, columns=["objra", "objdec"])
        sky = SkyCoord(coords["objra"].to_numpy().astype(float) * u.deg,
                       coords["objdec"].to_numpy().astype(float) * u.deg)
        inside = CENTER.separation(sky) < RADIUS
        if inside.any():
            table = partition.read_row_group(group)
            yield from table.filter(pa.array(inside)).to_pylist()


def exposure_metadata():
    with urlopen(EXPOSURES, timeout=90) as response:
        data = response.read()
    print(f"IBE response SHA-256 {hashlib.sha256(data).hexdigest()}")
    metadata = pd.read_csv(io.BytesIO(data), dtype={"obsjd": str})
    metadata["mjd"] = [float(Decimal(jd) - Decimal("2400000.5")) for jd in metadata.obsjd]
    # The IBE time is the start of the exposure. The bulk lightcurve reports
    # heliocentric MJD at mid-exposure to five decimal places.
    iers.conf.auto_download = False
    palomar = EarthLocation.from_geodetic(-116.862 * u.deg, 33.356 * u.deg, 1712 * u.m)
    start = Time(metadata.mjd.to_numpy(), format="mjd", scale="utc", location=palomar)
    correction = start.light_travel_time(CENTER, kind="heliocentric").to_value(u.day)
    metadata["hmjd"] = metadata.mjd + correction + metadata.exptime / 172800
    return metadata


def match_exposures(source, metadata):
    band = BANDS[source["filterid"]]
    channel = source["rcid"]
    field = metadata[(metadata.field == source["fieldid"])
                     & (metadata.ccdid == channel // 4 + 1)
                     & (metadata.qid == channel % 4 + 1)
                     & (metadata.filtercode == band)].sort_values("hmjd")
    epochs = pd.DataFrame(source["lightcurve"])
    times = field.hmjd.to_numpy()
    lo = np.searchsorted(times, epochs.hmjd - 1 / 86400)
    hi = np.searchsorted(times, epochs.hmjd + 1 / 86400, side="right")
    if np.any(hi - lo != 1) or len(set(lo)) != len(lo):
        raise ValueError(f"Ambiguous or missing exposure match for {source['objectid']}")
    matches = field.iloc[lo]
    epochs["mjd"] = matches.mjd.to_numpy()
    epochs["dt_sec"] = 86400 * np.abs(epochs.hmjd.to_numpy() - matches.hmjd.to_numpy())
    epochs["band"] = band
    return epochs


def nightly(epochs):
    medians = epochs.groupby(["image", "band", "day"])[["mjd", "mag", "magerr"]].median()
    return medians.rename(columns={"mag": "magnitude", "magerr": "median_catalog_error"})


def main():
    components = pd.read_csv(RESULTS / "targets.csv").query("system == @TARGET")
    positions = SkyCoord(components.ra_deg.to_numpy() * u.deg, components.dec_deg.to_numpy() * u.deg)
    metadata = exposure_metadata()
    detections = []
    for source in cone_lightcurves():
        reference = SkyCoord(source["objra"] * u.deg, source["objdec"] * u.deg)
        hits = components[reference.separation(positions).arcsec < 0.85]
        if len(hits) > 1:
            raise ValueError(f"Ambiguous image for object {source['objectid']}")
        if len(hits) == 1:
            epochs = match_exposures(source, metadata)
            epochs["image"] = hits.image.iloc[0]
            detections.append(epochs)
    epochs = pd.concat(detections, ignore_index=True)
    print(f"{len(detections)} light curves, {len(epochs)} detections, "
          f"maximum absolute exposure-time residual {epochs.dt_sec.max():.4f} s")
    # DR24 HATS omits the CGI per-detection limitmag. The comparison tests
    # whether the available cuts independently recover the published night set.
    clean = epochs[(epochs.catflags == 0) & (epochs.magerr > 0) & (epochs.magerr < 0.28)].copy()
    clean["day"] = np.floor(clean.mjd).astype(int)

    frozen = pd.read_csv(RESULTS / "ztf-nights.csv").query("system == @TARGET").copy()
    frozen["day"] = np.floor(frozen.mjd).astype(int)
    frozen = frozen.set_index(["image", "band", "day"]).sort_index()
    recovered = nightly(clean)
    if not recovered.index.equals(frozen.index):
        raise ValueError("Recovered image, band and night keys differ from the deposited sample")
    expected = frozen[recovered.columns]
    discrepancy = (expected - recovered).abs().max()
    print(f"{len(clean)} accepted detections, {len(expected)} distinct image-band nights")
    print("Largest absolute differences before decimal display rounding:")
    print(discrepancy.to_string(float_format=lambda value: f"{value:.9g}"))
    if (discrepancy > [5.1e-8, 5.6e-7, 5.1e-7]).any():
        raise ValueError("Native bulk photometry or exposure timestamps exceed audit tolerances")

    # Nine significant decimal digits preserve float32 photometry. Converting
    # to decimal before the median resolves ties at the published rounding unit.
    for column in ["mag", "magerr"]:
        clean[column] = clean[column].map(lambda value: float(format(value, ".9g")))
    displayed = nightly(clean).round({"mjd": 7, "magnitude": 6, "median_catalog_error": 6})
    pd.testing.assert_frame_equal(displayed, expected, check_exact=True)
    print("Every published MJD, magnitude and median catalog error matches exactly at printed precision.")


if __name__ == "__main__":
    main()
