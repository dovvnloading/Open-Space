"""Synthetic Gaia association, integer precision, catalog contract and HTTP-cache regressions."""
from pathlib import Path
import tempfile
from unittest import TestCase
from unittest.mock import patch
from urllib.parse import urlencode
from xml.etree import ElementTree as ET

import campaign
import gaia_properties as gaia


A = "4100000000000000000"
B = "4100000000000000001"
C = "4100000000000000002"


def votable(fields, rows, status="OK", idtype="long"):
    root = ET.Element("VOTABLE", version="1.3", xmlns="http://www.ivoa.net/xml/VOTable/v1.3")
    resource = ET.SubElement(root, "RESOURCE", type="results")
    ET.SubElement(resource, "INFO", name="QUERY_STATUS", value=status)
    table = ET.SubElement(resource, "TABLE")
    for field in fields:
        props = {"name": field, "datatype": idtype} if field == "source_id" else {
            "name": field, "datatype": "char", "arraysize": "*"}
        ET.SubElement(table, "FIELD", **props)
    data = ET.SubElement(ET.SubElement(table, "DATA"), "TABLEDATA")
    for values in rows:
        tr = ET.SubElement(data, "TR")
        for value in values:
            ET.SubElement(tr, "TD").text = value
    return ET.tostring(root)


def optical_row(s, arcsec):
    ra = arcsec/3600
    return {"catalog_values": {"source_id": s, "ra": str(ra), "dec": "0", "ref_epoch": "2016.0",
                "pmra": None, "pmdec": None, "phot_g_mean_mag": None},
            "offset_arcsec": arcsec, "within_xray_r90": arcsec <= 1,
            "derived": {"motion_propagated": False, "offset_event_epoch_arcsec": arcsec}}


def save_optical(root):
    lead = {"obsid": 12, "input_sha256": "a"*64, "trigger_key": "real/1", "event_mjd": 59000,
            "ra_deg": 0, "dec_deg": 0, "position_error_arcsec": 1}
    prior, inside, outside = optical_row(A, 3), optical_row(B, 0.5), optical_row(C, 4)
    p = root/"lead_optical"/"observations"/"12"/lead["input_sha256"]/"real"/"1.json.gz"
    p.parent.mkdir(parents=True)
    campaign.save_obs(p, {"schema": 1, "lead": lead, "prior_gaia_ids": [A], "complete": True, "catalogs": {
        "gaia_cone": {"status": "rows_returned", "rows": [inside, outside]},
        "gaia_id_"+A: {"status": "rows_returned", "rows": [prior]}}})
    return p


class GaiaPropertiesTest(TestCase):
    def test_exact_identifier_schema_and_catalog_coverage(self):
        for bad in (float(A), False, "4.1e18", "041", str(2**63)):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                gaia.sid(bad)
        fields = ["source_id", "phot_g_mean_mag"]
        result = gaia.read_table(votable(fields, [(A, "18.0"), (B, None)]), fields, [A, B, C], gaia.SOURCE)
        self.assertEqual(result[A][0]["phot_g_mean_mag"], "18.0")
        self.assertIsNone(result[B][0]["phot_g_mean_mag"])
        self.assertEqual(result[C], [])
        for raw in (votable(fields, [(A, "18")], idtype="double"), votable(fields, [], status="OVERFLOW"),
                    votable(fields, [(C, "2")]), votable(fields, [(A, "2"), (A, "2")])):
            with self.assertRaises(ValueError):
                gaia.read_table(raw, fields, [A], gaia.SOURCE)
        schema = [dict(table_name=gaia.SOURCE, column_name=f, datatype=t) for f, t in
                  (("source_id", "adql:BIGINT"), ("phot_g_mean_mag", "REAL"))]
        self.assertEqual(gaia.schema_columns(schema, gaia.SOURCE), fields)
        with self.assertRaisesRegex(ValueError, "unavailable"):
            gaia.schema_columns(schema, gaia.VAR)

    def test_association_routes_through_original_optical_file(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path = save_optical(root)
            lead, = gaia.associations(root)
            self.assertEqual(set(lead["associations"]), {A, B})
            self.assertTrue(lead["associations"][A]["primary_match_records"])
            self.assertFalse(lead["associations"][A]["positions"][0]["within_xray_r90"])
            self.assertTrue(lead["associations"][B]["positions"][0]["within_xray_r90"])
            self.assertEqual(lead["associations"][B]["positions"][0]["motion_propagated"], False)
            self.assertAlmostEqual(lead["associations"][B]["positions"][0]["epoch_ra_deg"], 0.5/3600)
            self.assertEqual(lead["optical_records"][0]["sha256"], gaia.digest(path.read_bytes()))

    def test_incremental_requests_and_offline_replay_keep_absence_and_failure_distinct(self):
        calls = []
        def transport(label, endpoint, src, replay, memos, params):
            self.assertFalse(replay)
            q = params["QUERY"]
            calls.append(q)
            url = endpoint+"?"+urlencode(params)
            if q == gaia.SCHEMA_QUERY:
                fields = gaia.SCHEMA_FIELDS
                schema = {gaia.SOURCE: ("source_id", "phot_g_mean_mag", gaia.DSC[0]),
                          gaia.VAR: ("source_id", "best_class_name", "best_class_score")}
                raw = votable(fields, [(tab, col, "long" if col == "source_id" else "char", None, None)
                                       for tab, names in schema.items() for col in names])
            else:
                fields = q.split(" FROM ")[0].removeprefix("SELECT ").split(",")
                ids = q.split(" IN (")[1][:-1].split(",")
                if gaia.VAR in q and C in ids:
                    memos[label] = {"url": url, "ok": False, "attempts": []}
                    raise ValueError("Synthetic DNS failure")
                rows = ([(x, "18.0", None) for x in ids if x != B] if gaia.SOURCE in q else
                        [(x, "AGN", ".54") for x in ids if x == A])
                raw = votable(fields, rows)
            h = gaia.digest(raw)
            name = "fixture-"+h+".bin"
            (src/name).write_bytes(raw)
            receipt = {"body_file": name, "sha256": h, "size_bytes": len(raw)}
            memos[label] = {**receipt, "url": url, "ok": True, "attempts": [receipt]}
            return raw

        with tempfile.TemporaryDirectory() as d, patch.object(gaia, "read_input", side_effect=transport):
            out = Path(d)
            sources, _ = gaia.lookup(out, [A, B], False, "10-1")
            self.assertEqual(sources[B][gaia.SOURCE]["status"], "no_returned_source_row")
            self.assertEqual(sources[B][gaia.VAR]["status"], "no_published_row")
            self.assertIsNone(sources[A][gaia.SOURCE]["rows"][0][gaia.DSC[0]])
            calls.clear()
            updated = gaia.lookup(out, [A, B, C], False, "11-1")
            sources, schemas = updated
            self.assertTrue(all(q == gaia.SCHEMA_QUERY or f" IN ({C})" in q for q in calls))
            self.assertEqual(sources[C][gaia.VAR]["status"], "unavailable")
            self.assertTrue(schemas[-1]["tables"][gaia.VAR]["failures"][0]["receipt"])
            with patch.object(gaia, "read_input", side_effect=AssertionError("offline")):
                self.assertEqual(gaia.lookup(out, [A, B, C], True, "12-1"), updated)
                success = gaia.campaign.read(out/sources[A][gaia.SOURCE]["receipt"])
                raw = out/"source"/success["run_attempt"]/success["http"]["body_file"]
                raw.write_bytes(b"tampered")
                with self.assertRaisesRegex(ValueError, "Altered Gaia response"):
                    gaia.lookup(out, [A, B, C], True, "12-1")
