"""Reader regression checks for publisher-deleted and active FORMAT 1 rows."""
from pathlib import Path
from tempfile import TemporaryDirectory

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
            "INCLUDE other.tim\n"
            "c 1400 55755.0 1 @ -junk -sys telescope -addsat -1\n"
            "SKIP\n"
            "ignored 1400 55754.0 1 @\n"
            "NOSKIP\n"
            "END\n"
            "ignored 1400 55754.0 1 @\n")
        (work / "other.tim").write_text("FORMAT 1\nother 700 55753.1 1 @\n")
        receipt = {"normalized_comments":0, "clock_valid_mjd_by_site":{}, "clock_cuts":{},
                   "tim_file_counts":{}, "ignored_line_locations":{}}
        staged = work / "staged.tim"
        staged.write_text(stage_tim(original,work,receipt,set()))
        toas, _ = read_toa_file(staged)
        assert len(toas) == 2
        moved = [t for t in toas if t.flags.get("sys") == "telescope"]
        assert len(moved) == 1
        assert abs((moved[0].mjd.mjd - 55755.)*86400 + 1.) < 1e-6
        counts = receipt["tim_file_counts"]["source.tim"]
        assert counts["tempo2_I_deleted_toas"] == 1
        assert counts["incomplete_nondirective_lines"] == 1
        assert counts["valueless_flags_ignored"] == counts["addsat_flags_applied"] == 1


if __name__ == "__main__":
    main()
