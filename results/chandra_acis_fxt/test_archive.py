"""Offline archive regression fixtures for the manual survey workflow."""
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

import archive


def make_row(obsid, **changes):
    return {
        "obsid": str(obsid), "instrument": "ACIS-S", "exposure_mode": "TE",
        "grating": "NONE", "status": "archived", "public_avail_date": "2021-01-01T00:00:00",
        "start_date": "2020-03-17T01:02:03", "exposure_time": "4.125", "event_count": "700",
        "target_name": "fixture target", "ra": "0", "dec": "-25.5",
        "proposal_number": "TEST01", "proposal_type": "TOO", "proposal_title": "follow-up",
        "proposal_category": "fixture field", "proposal_pi": "Example", **changes,
    }


def fixture_download(content):
    def download(url, dest, data=None):
        dest.write_bytes(content)
        return {"url": url, "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
    return download


class TestArchive(unittest.TestCase):
    def test_release_window_imaging_and_observing_intent(self):
        rows = [
            make_row(10),  # Observed in 2020, released on first inclusive boundary.
            make_row(11, public_avail_date="2026-11-21T23:59:59"),
            make_row(12, public_avail_date="2026-11-22", start_date="2026-01-01"),
            make_row(13, public_avail_date="2020-12-31"),
            make_row(14, exposure_mode="CC"), make_row(15, grating="HETG"),
            make_row(16, public_avail_date=""), make_row(17, status="observed"),
            make_row(18, exposure_mode="UNSPECIFIED", proposal_number=""),
        ]
        pages = [([{"n": str(len(rows))}], {}), (rows[:5], {}), (rows[5:], {}), ([], {})]
        with tempfile.TemporaryDirectory() as path, patch.object(archive, "tap", side_effect=pages) as tap:
            directory = Path(path)
            selected = archive.inventory(directory, "2021-01-01", "2026-11-21")
            self.assertEqual([r["obsid"] for r in selected], [10, 11, 18])
            self.assertEqual(json.loads((directory / "inventory.json").read_text()), selected)
            first = selected[0]
            self.assertEqual(first["target_name"], "fixture target")
            self.assertEqual((first["catalog_ra_deg"], first["catalog_dec_deg"]), (0.0, -25.5))
            self.assertEqual(first["proposal_type"], "TOO")
            self.assertEqual(first["proposal_title"], "follow-up")
            self.assertEqual(first["proposal_pi"], "Example")
            self.assertEqual(first["proposal_category"], "fixture field")
            self.assertEqual(first["exposure_s"], 4125)
            self.assertEqual(parse_qs(urlsplit(first["observation_url"]).query)["obsid"], ["10"])
            self.assertIn("TEST01", parse_qs(urlsplit(first["proposal_url"]).query)["QUERY"][0])
            self.assertIsNone(selected[-1]["proposal_url"])
            exceptions = json.loads((directory / "inventory_exceptions.json").read_text())
            self.assertEqual({r["obsid"]: r["reason"] for r in exceptions}, {
                "16": "missing_public_date", "17": "not_marked_archived", "18": "check_mode_from_FITS",
            })
            self.assertIn("o.obsid > 14", tap.call_args_list[2].args[2])

    def test_latest_interleaved_streams_preserve_both_cycles(self):
        names = [
            "acisf00010N002_e1_evt2.fits.gz", "acisf00010N003_e1_evt2.fits",
            "acisf00010N003_e1_evt2.fits.gz", "acisf00010N003_e2_evt2.fits.gz",
            "acisf00011N009_evt2.fits.gz", "acisf00010N003_evt1.fits.gz",
            "hrcf00010N004_evt2.fits.gz",
        ]
        html = "".join(f'<a href="{name}">download</a>' for name in reversed(names))
        self.assertEqual(archive.product_names(html, 10), names[2:4])

    def test_partial_catalog_fails(self):
        pages = [([{"n": "3"}], {}), ([make_row(10)], {}), ([], {})]
        with tempfile.TemporaryDirectory() as path, patch.object(archive, "tap", side_effect=pages):
            with self.assertRaisesRegex(ValueError, "1 of 3"):
                archive.inventory(Path(path), "2021-01-01", "2026-11-21")
            self.assertFalse((Path(path) / "inventory.json").exists())

    def test_tap_trailing_error_disallows_partial_table(self):
        # A TAP server may emit an error after the streamed rows.
        raw = b'''<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE>
          <INFO name="QUERY_STATUS" value="OK"/><TABLE><FIELD name="obsid"/>
          <DATA><TABLEDATA><TR><TD>10</TD></TR></TABLEDATA></DATA></TABLE>
          <INFO name="QUERY_STATUS" value="ERROR">truncated</INFO></RESOURCE></VOTABLE>'''
        with tempfile.TemporaryDirectory() as path, patch.object(archive, "download", fixture_download(raw)):
            with self.assertRaisesRegex(ValueError, "truncated"):
                archive.tap(Path(path), "cut", "SELECT TOP 1 obsid FROM cxc.observation")
            self.assertEqual(gzip.decompress((Path(path) / "cut.xml.gz").read_bytes()), raw)

    def test_fits_identity_pointing_and_imaging_mode(self):
        from astropy.io import fits

        events = fits.BinTableHDU.from_columns(
            [fits.Column(name="TIME", format="D", array=[10.0, 11.0])], name="EVENTS"
        )
        events.header.update({
            "OBS_ID": 10, "INSTRUME": "ACIS", "GRATING": "NONE", "READMODE": "TIMED",
            "RA_NOM": 1.0, "DEC_NOM": -25.0, "RA_TARG": 1.1, "DEC_TARG": -25.1,
        })
        with tempfile.TemporaryDirectory() as path:
            product = Path(path) / "fixture.fits.gz"
            hdus = fits.HDUList([fits.PrimaryHDU(), events])
            hdus.writeto(product)
            receipt = archive.inspect_product(product, 10)
            self.assertEqual(receipt["rows"], 2)
            self.assertEqual(receipt["fits_sha256"], hashlib.sha256(gzip.decompress(product.read_bytes())).hexdigest())
            self.assertEqual([receipt[k] for k in ("ra_nom_deg", "dec_nom_deg", "ra_targ_deg", "dec_targ_deg")],
                             [1.0, -25.0, 1.1, -25.1])
            with self.assertRaisesRegex(ValueError, "OBS_ID"):
                archive.inspect_product(product, 20)
            events.header["GRATING"] = "LETG"
            hdus.writeto(product, overwrite=True)
            with self.assertRaisesRegex(ValueError, "imaging configuration"):
                archive.inspect_product(product, 10)

    def test_exact_replay_and_same_length_revision_change(self):
        # Same length makes this exercise the checksum, not just truncation.
        original = b"fixture A"
        name = "acisf00010N003_evt2.fits.gz"
        with tempfile.TemporaryDirectory() as path:
            directory = Path(path)
            manifest = directory / "manifest.json"
            archive.dump(manifest, [{"url": archive.ROOTS[0] + "/0/10/primary/" + name,
                                     "filename": name, "bytes": len(original),
                                     "sha256": hashlib.sha256(original).hexdigest()}])
            output = directory / "products"
            with patch.object(archive, "download", fixture_download(original)):
                archive.replay(manifest, output)
            with patch.object(archive, "download", fixture_download(b"fixture B")):
                with self.assertRaisesRegex(ValueError, "differs from manifest"):
                    archive.replay(manifest, output)
            self.assertEqual((output / name).read_bytes(), original)
            self.assertFalse((output / (name + ".verified")).exists())


if __name__ == "__main__":
    unittest.main()
