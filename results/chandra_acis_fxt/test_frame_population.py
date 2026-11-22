from collections import Counter
import gzip
from pathlib import Path
import tempfile
import unittest

import frame_population as survey


class RecordedS4Population(unittest.TestCase):
    def test_identical_stream_selection_and_keyed_real_rows(self):
        e={"ccd_id":8,"instrument_pass":True,"max_counts_per_frame":3,"on_counts":4,"peak_counts":3,
           "distinct_frames":2,"ra_deg":10.,"dec_deg":30.,"trial_duration_s":600,"seed_score":9.,
           "trigger_key":"real/0","vet":{"status":"no_counterpart"}}
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            path=root/"photometry/1/event/real/0/event_rows.csv.gz"
            path.parent.mkdir(parents=True)
            table=("region,phase,chipx,chipy,time_met_s,in_peak\n"
                   "aperture,trigger,13,40,18.0,1\naperture,trigger,15,40,18.0,1\n"
                   "aperture,trigger,16,41,18.0,1\naperture,trigger,16,41,180,0\n"
                   "aperture,before,13,40,0,0\nannulus,trigger,22,40,18.0,1\n")
            path.write_bytes(gzip.compress(table.encode(),mtime=0))
            e["diagnostics"]={"files":{"event_rows.csv.gz":survey.ingest.sha(path)}}
            controls=[{"replica":j,"triggers":[{**e,"trigger_key":f"control-{j}/0"}]} for j in range(24)]
            p={"input":{"filename":"event","sha256":"a"*64},"search":{"triggers":[e,{**e,"ccd_id":7},
               {**e,"instrument_pass":False},{**e,"max_counts_per_frame":2}],"controls":controls}}
            counts={s:Counter() for s in survey.STREAMS}
            rows=survey.collect_product(root,1,p,counts,"b"*64)
            self.assertEqual(len(rows),25)
            self.assertEqual(counts["real"],{"raw_s4":3,"instrument_pass_s4":2,"selected_s4":1,"selected_after_catalog":1})
            self.assertEqual((rows[0]["max_same_node_row_and_time"],rows[0]["time_tags_with_repeated_node_row"]),(2,1))
            self.assertEqual(rows[-1]["row_audit"],"no archived control photon table")
            path.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError,"hash differs"):
                survey.collect_product(root,1,p,counts,"b"*64)
