"""Offline checks for the preregistered real-FITS byte-range reader.

Mocks actual FITS byte offsets/WCS; tests do not masquerade as ACT observations.
"""
from __future__ import annotations
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_preregistered_field as s


class RealFitsPilotSyntheticTests(unittest.TestCase):
    def test_registered_subset_catalog_has_no_snr_selection(self):
        tests, refs = s.selection()
        self.assertEqual(len({str(r["map_id"]) for r in tests+refs}), len(tests)+len(refs))
        self.assertTrue(all(r["tube_slot"] == "pa4" and r["frequency"] == "f150" for r in tests+refs))
        self.assertTrue(all(s.in_interval(r["start_time"], s.DATE_START, s.DATE_END) for r in tests))
        self.assertTrue(all(not s.in_interval(r["start_time"], s.DATE_START, s.DATE_END) for r in refs))
        self.assertTrue(all(s.in_interval(r["start_time"], s.REF_START, s.REF_END) for r in refs))
        self.assertTrue(any(r["start_time"].startswith("2018-09-10") for r in tests))
        self.assertFalse(any(r["start_time"].startswith("2018-09-10") for r in refs))

    def test_space_separated_catalog_timestamps_obey_calendar_day_boundaries(self):
        self.assertTrue(s.in_interval("2018-09-10 23:06:59", s.DATE_START, s.DATE_END))
        self.assertFalse(s.in_interval("2018-09-10 23:06:59", s.REF_START, s.DATE_START))
        self.assertFalse(s.in_interval("2018-09-13 00:00:00", s.DATE_START, s.DATE_END))
        self.assertTrue(s.in_interval("2018-09-13 00:00:00", s.REF_START, s.REF_END))

    def test_gaussian_beam_fwhm_uncertainty_is_propagated(self):
        psfs, receipt = s.beam_family()
        self.assertEqual(sorted(psfs), [1.2,1.4,1.6])
        self.assertEqual(receipt["fwhm_sensitivity_envelope_arcmin"], [1.2,1.6])
        self.assertFalse(receipt["validated_for_single_observation"])
        self.assertTrue(psfs[1.6].sum() > psfs[1.2].sum())

    def test_mapcat_ctime_string_is_not_unix_float(self):
        from datetime import datetime, timezone
        dt = datetime.fromisoformat("2018-09-11 23:06:30.000000").replace(tzinfo=timezone.utc)
        self.assertGreater(dt.timestamp(), 1_500_000_000)

    def test_real_fits_data_byte_offsets_for_t_plane(self):
        w = WCS(naxis=2)
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.crval = [s.CENTER_RA, s.CENTER_DEC]
        w.wcs.crpix = [101., 101.]
        w.wcs.cdelt = [-s.PIX_DEG, s.PIX_DEG]
        ary = np.zeros((3,201,201), dtype="float32")
        ary[0,:,:] = 9.
        ary[1,:,:] = 1000.
        ary[2,:,:] = -2000.
        buff = io.BytesIO()
        fits.PrimaryHDU(data=ary, header=w.to_header()).writeto(buff)
        payload = buff.getvalue()
        def fake(url, start, stop, header=False):
            self.assertLess(stop, len(payload))
            return payload[start:stop+1], {
                "status":206, "range":f"bytes {start}-{stop}/{len(payload)}"}
        ra, dec, aperture = s.tile_wcs()
        with patch.object(s, "fetch_range", side_effect=fake):
            data, receipt = s.access_tile("path/fixture.fits", ra, dec)
        self.assertTrue(receipt["roi_overlap"])
        self.assertGreater(receipt["science_bytes"], 0)
        self.assertEqual(receipt["status"], 206)
        self.assertAlmostEqual(data[80,80], 9., places=4)
        self.assertTrue(np.allclose(data[aperture], 9., atol=1e-4, equal_nan=False))
        self.assertIn("ordered_range_digests_sha256", receipt)

    def test_reference_requires_three_disjoint_nights(self):
        def pair(value, ivar=1.):
            return np.full((3,3),value),np.full((3,3),ivar)
        recs = [
            ({"start_time":"2018-09-04T01:00:00"}, pair(3)),
            ({"start_time":"2018-09-04T02:00:00"}, pair(5)),
            ({"start_time":"2018-09-05T01:00:00"}, pair(0)),
        ]
        _, iv, nights = s.ref_coadd(recs)
        self.assertEqual(len(nights), 2)
        self.assertEqual(int(iv.sum()), 0)
        recs.append(({"start_time":"2018-09-06T01:00:00"}, pair(0)))
        ref, iv, nights = s.ref_coadd(recs)
        self.assertEqual(len(nights), 3)
        self.assertTrue(np.all(iv == 3))  # same-night maps add no independent ivar
        self.assertTrue(np.allclose(ref, 4/3, atol=1e-12))

    def test_aperture_area_not_full_fits_bounding_box(self):
        _, _, mask = s.tile_wcs()
        area = mask.sum()*(s.PIX_DEG**2)
        self.assertAlmostEqual(area, np.pi*s.RADIUS_DEG**2, delta=.007)
        self.assertLess(area, s.N**2*s.PIX_DEG**2)


if __name__ == "__main__":
    unittest.main()
