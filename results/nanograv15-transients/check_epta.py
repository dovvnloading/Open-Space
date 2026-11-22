"""Reader regression checks for publisher-deleted and active FORMAT 1 rows."""
from pathlib import Path
from tempfile import TemporaryDirectory

import astropy.units as u
from pint.toa import read_toa_file

from epta import stage_tim


def main():
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        original = work / "source.tim"
        original.write_text(
            "FORMAT 1\n"
            "Cactive 1400 55752.0 1 @\n"
            "I???? original_filename 1400 55753.0 1 @\n"
            "C???? original_filename 1400 55754.0 1 @\n"
            "edited by hand\n"
            "TIME 3\n"
            "EFAC 2\n"
            "EQUAD 4\n"
            "INCLUDE other.tim\n"
            "c 1400 55755.0 3 @ -junk -sys telescope -addsat -1\n"
            "SKIP\n"
            "ignored 1400 55754.0 1 @\n"
            "NOSKIP\n"
            "END\n"
            "ignored 1400 55754.0 1 @\n")
        child = work / "other.tim"
        child.write_text("FORMAT 1\nother 700 55753.1 1 @ -sys child\n"
                         "TIME 10\nFMIN 800\ncut 700 55753.2 1 @\n"
                         "other 1400 55754.0 1 @ -sys later\n")
        receipt = {"normalized_comments":0, "clock_valid_mjd_by_site":{}, "clock_cuts":{},
                   "tim_file_counts":{}, "ignored_line_locations":{}}
        staged = work / "staged.tim"
        staged.write_text(stage_tim(original,work,receipt,set()))
        child.write_text(stage_tim(child,work,receipt,set()))
        toas, _ = read_toa_file(staged)
        by_sys = {t.flags["sys"]:t for t in toas}
        assert len(toas) == len(by_sys) == 3
        assert abs((by_sys["child"].mjd.mjd - 55753.1)*86400) < 1e-6
        assert abs((by_sys["later"].mjd.mjd - 55754.)*86400 - 10.) < 1e-6
        assert abs((by_sys["telescope"].mjd.mjd - 55755.)*86400 - 2.) < 1e-6
        assert by_sys["telescope"].error.to_value(u.us) == 10
        assert by_sys["child"].error.to_value(u.us) == by_sys["later"].error.to_value(u.us) == 1
        counts = receipt["tim_file_counts"]["source.tim"]
        assert counts["tempo2_I_deleted_toas"] == 1
        assert counts["incomplete_nondirective_lines"] == 1
        assert counts["valueless_flags_ignored"] == counts["addsat_flags_applied"] == 1
        assert receipt["tim_file_counts"]["other.tim"]["local_window_cuts"] == 1


if __name__ == "__main__":
    main()
