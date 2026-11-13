import tempfile
import unittest
from unittest.mock import patch

from open_space import wise, wise_archive

FIELD = {"id": "test", "ra_deg": 359.9999, "dec_deg": 30.0, "radius_arcmin": 3.0, "role": "exploration"}


def rows(pm_x=0.5, pm_y=-0.32, seasons=20):
    values = []
    for i in range(seasons):
        for j in range(6):
            mjd = wise.REF_MJD + (i - 5) * 182.625 + j / 12
            t = (mjd - wise.REF_MJD) / wise.DAYS_YEAR
            ra, dec = wise.sky(pm_x * t + (j - 2.5) * 0.004, pm_y * t, FIELD)
            values.append(
                {
                    "source_id": f"{i:02d}-{j}",
                    "scan_id": f"{i:02d}",
                    "frame_num": str(j),
                    "ra": ra,
                    "dec": dec,
                    "mjd": mjd,
                    "w1mpro": 12,
                    "w1sigmpro": 0.03,
                    "w1snr": 25,
                    "w1rchi2": 1,
                    "w1sat": 0,
                    "cc_flags": "0000",
                    "qual_frame": 10,
                    "qi_fact": 1,
                    "moon_masked": "00",
                    "saa_sep": 18,
                    "sso_flg": 0,
                    "nb": 1,
                    "na": 0,
                    "sigra": 0.1,
                    "sigdec": 0.1,
                    "allwise_cntr": "null",
                }
            )
    return values


def response(rows, count, timestamp="one"):
    names = wise_archive.COLUMNS["gaia"]
    return (
        f"\\DATETIME= '{timestamp}'\n\\RowsRetrieved = {count}\n"
        + "|"
        + "|".join(k.ljust(16) for k in names)
        + "|\n"
        + "".join(" " + " ".join(str(v).ljust(16) for v in row) + "\n" for row in rows)
    ).encode()


class WiseTests(unittest.TestCase):
    def test_moving_source_crosses_ra_zero_and_matches_gaia(self):
        gaia = [
            {
                "source_id": "7241",
                "ra": FIELD["ra_deg"],
                "dec": 30.0,
                "pmra": 501,
                "pmdec": -318,
                "ref_epoch": 2016.0,
                "ruwe": 1.09,
            }
        ]
        sources, counts = wise.detect(rows(), gaia, FIELD)
        s = sources[0]
        self.assertEqual(counts["qualified_tracks"], 1)
        self.assertAlmostEqual(s["pmra_mas_yr"], 500, delta=1)
        self.assertAlmostEqual(s["pmdec_mas_yr"], -320, delta=1)
        self.assertEqual(s["disposition"], "catalog_motion_recovered")
        self.assertEqual(s["quality"]["epochs"], 20)
        self.assertTrue(s["quality"]["compatible_phases"])

    def test_static_and_low_coverage(self):
        s, _ = wise.detect(rows(pm_x=0, pm_y=0), [], FIELD)
        self.assertEqual(s[0]["disposition"], "below_motion_threshold")
        s, _ = wise.detect(rows(seasons=6), [], FIELD)
        self.assertEqual(s, [])
        s, _ = wise.detect(rows(), [], FIELD)
        self.assertEqual(s[0]["disposition"], "requires_catalog_and_image_vetting")

    def test_quality_and_distinct_frames(self):
        rs = rows()
        for r in rs:
            r["frame_num"] = "0"
        sources, counts = wise.detect(rs, [], FIELD)
        self.assertEqual(sources, [])
        self.assertEqual(counts["visit_sources"], 0)
        rs[0]["cc_flags"] = "D000"
        self.assertEqual(wise.quality_reason(rs[0]), "exposure_or_artifact_flag")
        rs[0]["cc_flags"] = "0000"
        rs[0]["sso_flg"] = 1
        self.assertEqual(wise.quality_reason(rs[0]), "saturation_blend_or_fit")

    def test_confused_catalog_match_is_flagged(self):
        pair = [
            {
                "ra": FIELD["ra_deg"],
                "dec": 30,
                "pmra": 500,
                "pmdec": -320,
                "ref_epoch": 2016,
                "ruwe": 1,
                "source_id": str(i),
            }
            for i in (100, 101)
        ]
        source = wise.detect(rows(), pair, FIELD)[0][0]
        self.assertEqual(source["gaia"]["matches_within_2arcsec"], 2)
        self.assertEqual(source["disposition"], "catalog_comparison_followup")

    def test_canonical_archive_preserves_ids_and_discards_request_noise(self):
        a = ("23456789012345", "359.99", "30", "null", "null", "2016", "null")
        b = ("13456789012345", "359.89", "30", "0.03", "1", "2016", "1.7")
        x = wise_archive.canonical_table(response([a, b], 2), "gaia")
        assert x == wise_archive.canonical_table(response([b, a], 2, "two"), "gaia")
        self.assertIn(b"13456789012345,359.89", x)
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            wise_archive.canonical_table(response([a], 2), "gaia")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            wise_archive.canonical_table(response([b, b], 2), "gaia")
        with self.assertRaisesRegex(ValueError, "row limit"):
            wise_archive.canonical_table(response([], wise_archive.MAX_ROWS + 1), "gaia")

    def test_cache_mismatch_and_input_limits(self):
        import io

        b = ("13456789012345", "359.89", "30", "0.03", "1", "2016", "1.7")
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(wise_archive, "urlopen", return_value=io.BytesIO(response([b], 1))):
                rows, metadata = wise_archive.acquire(FIELD, tmp, "gaia")
            self.assertEqual(rows[0]["source_id"], b[0])
            with self.assertRaisesRegex(ValueError, "verification"):
                wise_archive.acquire({**FIELD, "input_sha256": {"gaia": "a" * 64}}, tmp, "gaia")
            self.assertEqual(wise_archive.acquire(FIELD, tmp, "gaia")[1], metadata)
        with self.assertRaises(ValueError):
            wise_archive.field_params({**FIELD, "id": "../bad"})
        with self.assertRaises(ValueError):
            wise_archive.field_params({**FIELD, "radius_arcmin": float("nan")})
