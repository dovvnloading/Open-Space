"""Verify the committed PR changes and its one additional counted catalog entry."""
import re
import subprocess

import report


def git(*args):
    return subprocess.check_output(["git", *args])


def verify():
    base=git("merge-base", "origin/main", "HEAD").strip().decode("ascii")
    changed=git("diff", "--name-only", "-z", base, "HEAD").decode("utf-8").split("\0")
    allowed={"README.md", ".github/workflows/chandra-acis-fxt.yml"}
    unexpected=[p for p in changed if p and p not in allowed and not p.startswith("results/chandra_acis_fxt/")]
    if unexpected:
        raise ValueError(f"Unrelated paths in PR: {unexpected}")
    original=git("show", base+":README.md").decode("utf-8").splitlines(keepends=True)
    revised=git("show", "HEAD:README.md").decode("utf-8").splitlines(keepends=True)
    if len(original)<9 or len(revised)<10:
        raise ValueError("Unexpected root README layout")
    markers=r"(?<!\*)\*(Candidate|Measured|Prediction|Correction|Null result)\*(?!\*)"
    rows=[i for i,line in enumerate(revised) if line.startswith("- ") and report.TAG in line]
    if len(rows)!=1 or len(re.findall(markers,revised[rows[0]]))!=1:
        raise ValueError("Expected exactly one counted Chandra catalog entry")
    label=re.findall(markers,revised[rows[0]])[0]
    titles={"Measured":"new measurements","Candidate":"candidate discoveries","Prediction":"forecasts",
            "Correction":"corrections","Null result":"carefully recorded null"}
    totals=[re.findall(markers,line) for line in revised if line.startswith("- ")]
    counts=[name for names in totals for name in names]
    expected=original[8]
    for phrase,total in [("entries",len(counts)),*[(title,counts.count(k)) for k,title in titles.items()]]:
        expected,n=re.subn(r"\b\d+(?= "+re.escape(phrase)+r"\b)",str(total),expected)
        if n!=1:
            raise ValueError(f"Unrecognized root summary field: {phrase}")
    if revised[8]!=expected:
        raise ValueError("The ninth line differs from the expected label counts")
    revised[8]=original[8]
    index=rows[0]
    del revised[index]
    if revised!=original and index<len(revised) and revised[index].strip()=="":
        del revised[index]
    if revised!=original:
        raise ValueError("The root README differs beyond the count summary and Chandra entry")
    print(f"Verified the PR paths and the single {label} catalog entry relative to {base}.")


if __name__=="__main__":
    verify()
