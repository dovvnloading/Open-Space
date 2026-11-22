"""Candidate-only preservation, byte-exact reconstruction, and damage detection."""
import gzip
import hashlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from astropy.io import fits

import candidate_inputs as ci


class TestCandidateInputs(unittest.TestCase):
    def test_empty_review_emits_no_event_bytes_or_downloads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            decisions=root/"decisions.json"
            decisions.write_text('{"decisions": []}')
            (root/"publication.json").write_text("{}")
            with (patch.object(ci,"BASE",root/"candidate_inputs"),patch.object(ci.campaign,"EVIDENCE",root),
                  patch.object(ci,"reviewed",return_value=({},{"reviewed_inputs_sha256":"0"*64})),
                  patch.object(ci.archive,"replay",side_effect=AssertionError("empty review fetched data")) as http):
                first=ci.preserve(decisions,root,[])
                self.assertEqual(first,ci.preserve(decisions,root,[]))
                self.assertEqual(ci.manifest(ci.BASE)["products"],{})
                self.assertFalse((ci.BASE/"objects").exists())
                http.assert_not_called()

    def test_split_and_restore_verify_transport_and_fits_bytes(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(ci,"PART",73):
            root=Path(tmp)
            artifact=root/"artifact"
            artifact.mkdir()
            path=artifact/"acisf00042N003_evt2.fits.gz"
            table=fits.BinTableHDU.from_columns([fits.Column(name="TIME",format="D",array=[1.,2.])],name="EVENTS")
            table.header.update({"OBS_ID":"42","INSTRUME":"ACIS","READMODE":"TIMED","GRATING":"NONE","REVISION":3})
            raw=io.BytesIO()
            fits.HDUList([fits.PrimaryHDU(),table]).writeto(raw)
            path.write_bytes(gzip.compress(raw.getvalue(),mtime=0))
            src={"obsid":42,"filename":path.name,"bytes":path.stat().st_size,"sha256":ci.sha(path),
                 "fits_sha256":hashlib.sha256(raw.getvalue()).hexdigest(),"rows":2,"revision":3,
                 "readmode":"TIMED","grating":"NONE","url":"https://cxc.cfa.harvard.edu/cdaftp/byobsid/2/42/primary/"+path.name}
            with patch.object(ci.archive,"replay",side_effect=AssertionError("retained source was ignored")):
                chosen,method,_=ci.capture(src,[artifact],None,root)
            self.assertEqual(chosen,path)
            self.assertEqual(method["kind"],"retained_actions_artifact")
            store=root/"store"
            parts=ci.split(chosen,src,store)
            self.assertGreater(len(parts),1)
            entry={"archive_receipt":src,"parts":parts}
            restored=ci.restore(entry,store,root/"out")
            self.assertEqual(path.read_bytes(),restored.read_bytes())
            with self.assertRaisesRegex(ValueError,"revision"):
                ci.verify(path,{**src,"revision":4})
            blob=store/parts[0]["file"]
            blob.write_bytes(b"x"*parts[0]["bytes"])
            with self.assertRaisesRegex(ValueError,"original bytes"):
                ci.restore(entry,store,root/"second")
