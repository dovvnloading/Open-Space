"""PDS4 NASA primary archive decoder, independent label-offset unit tests."""
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    from examples.maven_pds_kp_primary_audit import CHANNELS, parse_pds, pds_url, OFFICIAL_MD5
    HAS_SCIENCE = True
except ImportError:
    HAS_SCIENCE = False


@unittest.skipUnless(HAS_SCIENCE, "needs the optional numpy+pandas validation environment")
class PDSPrimaryFormatTests(unittest.TestCase):
    def test_fixed_width_pds_field_locations(self):
        # Construct a very small genuine PDS4-style table with all 17 fields.
        names = ["Time (UTC/SCET)"] + list(CHANNELS.values())
        widths = [19] + [16] * len(CHANNELS)
        rec_length = sum(widths) + 2
        xml_fields = []
        start = 1
        for i, (name, w) in enumerate(zip(names, widths), 1):
            xml_fields.append(f"<Field_Character><name>{name}</name><field_number>{i}</field_number><field_location>{start}</field_location><field_length>{w}</field_length></Field_Character>")
            start += w
        label = f"<Product_Observational><File_Area_Observational><Table_Character><offset>0</offset><records>2</records><Record_Character><fields>{len(names)}</fields><record_length>{rec_length}</record_length>{''.join(xml_fields)}</Record_Character></Table_Character></File_Area_Observational></Product_Observational>"
        def record(t, n):
            return (t + "".join(f"{float(n + i):16.4E}" for i in range(16)) + "\r\n").encode("ascii")
        with tempfile.TemporaryDirectory() as tmp:
            tab = Path(tmp) / "unit.tab"
            xml = Path(tmp) / "unit.xml"
            tab.write_bytes(record("2025-08-24T05:00:00", 1000) + record("2025-08-24T05:00:08", 2000))
            xml.write_text(label)
            df, schema = parse_pds(tab, xml)
            self.assertEqual(len(df), 2)
            self.assertEqual(schema["record_length"], rec_length)
            self.assertEqual(df.iloc[0]["sep_ion_fov1f"], 1000)
            self.assertEqual(df.iloc[0]["sep_ion_fov1f_quality"], 1001)
            self.assertEqual(df.iloc[1]["sep_electron_fov2r_quality"], 2015)

    def test_primary_archive_is_specific_pds4_release_not_mirror(self):
        self.assertIn("mvn_kp_insitu_20250824_v26_r01", pds_url("tab"))
        self.assertEqual(OFFICIAL_MD5, "377f757a180b34c7af46b2937831c181")


if __name__ == "__main__":
    unittest.main()
