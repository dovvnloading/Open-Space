"""Run all archive queries for KUG 1256+241 in the dedicated GitHub workflow."""
import importlib
import os
import traceback
from concurrent.futures import ThreadPoolExecutor
from common import OUT, PREVIOUS, input_file, save, status


def main():
    if os.environ.get("GITHUB_ACTIONS") != "true":
        raise SystemExit("Run the dedicated GitHub Actions workflow for live collection.")
    # Include inherited material used in interpreting the three survey data points.
    for path in ["outputs/candidates.csv", "outputs/rass_hardness/source_diagnostics.csv",
                 "outputs/seyfert_literature.json", "outputs/broker_cones.json"]:
        input_file(PREVIOUS / path)
    def collect(name):
        try:
            importlib.import_module(name).run()
        except Exception as exc:
            status(name, "unavailable", f"Module stopped early: {exc}")
            (OUT / f"{name}-error.txt").write_text(traceback.format_exc())
        print(f"Completed {name}", flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(collect, ["spectra", "photometry", "xray", "literature"]))
    save()


if __name__ == "__main__":
    main()
