"""Emit one JSON document for the pinned pulsar timing benchmark."""
import argparse
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path

from verify import replay


class ArgumentError(ValueError):
    """Invalid arguments to the machine-readable entry point."""


class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentError(message)


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-html", type=Path, help="Pinned Jodrell Bank HTML file")
    parser.add_argument("--atnf-html", type=Path, help="Pinned ATNF HTML file")
    parser.add_argument("--paper-pdf", type=Path, help="Verify original Table 2 using pdftotext")
    try:
        args = parser.parse_args(argv)
        with redirect_stdout(sys.stderr):
            measured = replay(args.catalog_html, args.atnf_html, args.paper_pdf)
        response = {
            "schema_version": 1, "status": "verified", **measured,
            "claim_scope": {
                "measured": "Retrospective forecast error for time to the next reported glitch.",
                "prior_work": "The J1341 size–wait correlation and the seven Parkes glitches were previously published.",
                "discovery_status": "No new glitch or physical mechanism established.",
                "external_permutation": "Conditional on exchangeability of waits at fixed predictions; serial dependence can invalidate it.",
                "gain_bounds": "External measurement sensitivity with the historical fit held fixed; future-event uncertainty is unmeasured.",
            },
        }
        # Serialize before emitting any output. Non-finite results are a failure.
        payload = json.dumps(response, ensure_ascii=True, allow_nan=False)
    except Exception as exc:
        code = 2 if isinstance(exc, ArgumentError) else 1
        print(json.dumps({"schema_version": 1, "status": "error",
                          "error": {"type": type(exc).__name__, "message": str(exc)}}))
        return code
    print(payload)
    return 0


if __name__ == "__main__":
    sys.exit(main())
