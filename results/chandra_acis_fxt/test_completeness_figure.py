import tempfile
from pathlib import Path
import unittest

from completeness_figure import render, table


class MeasuredCompletenessFigure(unittest.TestCase):
    def test_pools_widths_by_actual_attempts_and_preserves_empty_cells(self):
        with tempfile.TemporaryDirectory() as work:
            folder=Path(work)
            grid=folder/"grid.csv"
            grid.write_text("duration_s,input_counts,psf_sigma_sky_px,scheduled,attempted,recovered\n"
                            "60,8,1.2,2,2,2\n60,8,6,10,8,3\n180,16,3,1,0,0\n")
            good,total=table(grid)
            self.assertEqual((good[0,0],total[0,0]),(5,10))
            self.assertEqual(total[1,1],0)
            output=folder/"recovery.png"
            render(grid,output,10,5)
            self.assertEqual(output.read_bytes()[:8],b"\x89PNG\r\n\x1a\n")
            with self.assertRaisesRegex(ValueError,"survey totals"):
                render(grid,output,11,5)
